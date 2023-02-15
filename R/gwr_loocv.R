


#' Fits a GWR regression to tabular data contained in a data.frame
#'
#' This function fits a gwr model to a data.frame which includes geographic information
#'
#'
#' @param form a standard R formula in character format e.g. z~x+y
#' @param rdata a data.frame containing the variables in 'form', along with columns named 'lat' and 'lon'.
#' @param name a character vector used in any filenames created.
#' @param out.dir a filepath for outputs
#' @param abw adaptive bandwidth size.  Sets bandwidth as radius required to capture n neighbors.  Can be left blank if 'bw' is specified.
#' @param bw fixed bandwidth in km.  Can be left blank if 'abw' is set.  Error will be triggered if both 'bw' and 'abw' are set.
#' @param idx optional vector of row numbers at which to perform LOOCV.
#' @param n.samp optional max number of locations at which to perform LOOCV.
#' @return a gwr object with the original inputs, original data, fit, statistics.
#' @export
gwr_loocv <- function(form,data,name="",bw=NaN,adaptive=F,idx=NULL,n.samp=3000,verbose=T,buffer=3,fast=T){
  t0 <- unclass(Sys.time())
  if(!all(c("lon","lat") %in% colnames(data))) stop("need 'lon' and 'lat' colnames in data")
  if(is.null(idx)) {idx <- if(nrow(data)>n.samp) sample(nrow(data),n.samp) else 1:nrow(data)}
  '%dopar%' <- foreach::'%dopar%'
  '%do%' <- foreach::'%do%'

  if(!test_cluster(F)) stop("need cluster")
  if(foreach::getDoParWorkers()==1)  stop("need cluster")
  if(name=="") name <- paste0(gsub("\\+","_",sub("~","_vs_",form)),"_",ifelse(adaptive,"abw","bw"),bw)
  rn <- strsplit(form,"~")[[1]][1]
  vn <- all.vars(as.formula(form))[-1]
  data <- data[,c("id","lon","lat",vn,rn)]
  rownames(data) <- 1:nrow(data)
  ss <- complete.cases(data[,c(rn,vn)])
  data <- data[ss,]
  idx <- idx[!idx %in% which(!ss)]

  if(fast){
    X <- model.matrix(as.formula(form),data[,])
    np <- ncol(X)
    nc <- np+2
    XY <- as.matrix(terra::crds(terra::project(terra::vect(data[ss,],crs=geog),my_aea)))
  } else {
    np <-  ncol(model.matrix(as.formula(x@form),x@data[1:3,]))
    nc <- np*4+2
  }

  t1 <- unclass(Sys.time())
  if(verbose) cat("making LOO predictions for",form,"at",length(idx),"pts using cluster with",foreach::getDoParWorkers(),"cores\n")
  system.time(Out <- foreach(i=1:length(idx),.errorhandling='pass') %dopar% {
    if(fast){
      d <- sqrt((XY[,1]-XY[idx[i],1])^2+(XY[,2]-XY[idx[i],2])^2)/1000
      use_bw <- if(adaptive) sort(d)[bw] else bw
      md <- if(adaptive) sort(d)[bw*buffer] else bw*buffer

      # loo fit ~~~~~~~~~~~~~~~~~~
      ss2 <- d<md & d>0 & ss
      z <- if(np==1) try(lm.wfit(x=matrix(X[ss2,],sum(ss2)),y=data[ss2,rn],w=dnorm(d[ss2]/use_bw))$coefficients) else  try(lm.wfit(x=X[ss2,],y=data[ss2,rn],w=dnorm(d[ss2]/use_bw))$coefficients)
      y1 <- if(inherits(z,"try-error")) rep(NA,np) else z
      y1 <- if(length(y1)<np) rep(NA,np) else y1

      # full fit ~~~~~~~~~~~~~~~~~~
      ss3 <- d<md & ss
      z <-  if(np==1) try(lm.wfit(x=matrix(X[ss3,],sum(ss3)),y=data[ss3,rn],w=dnorm(d[ss3]/use_bw))$coefficients) else  try(lm.wfit(x=X[ss3,],y=data[ss3,rn],w=dnorm(d[ss3]/use_bw))$coefficients)
      y2 <- if(inherits(z,"try-error")) rep(NA,np) else z
      y2 <- if(length(y2)<np) rep(NA,np) else y2
      n <- sum(ss2)
      out <- c(y1,y2,bw=use_bw,n=n)
    } else {
      data$d <- as.vector(geodist::geodist(data[,c("lon","lat")],data[i,c("lon","lat")],measure="haversine"))/1000
      use_bw <- if(adaptive) sort(data$d)[bw+1] else bw
      data$wts <- dnorm(data$d/use_bw)
      ss2 <- data$d<use_bw*buffer & data$d>0 & ss
      if((n<-sum(ss2))>10){
        nu <- if(n<100) qr(data[ss2,vn])$rank else 9999
        if(nu>=np){
          fit <- predict(m0 <- lm(form,weights=wts,data=data[ss2,]),newdata=data[i,])
          out <- list(fit,summary(m0)$coefficients,n,nu,use_bw)
        } else out <- list(NA,matrix(NA,np,4),n,nu,use_bw)
      } else out <- list(NA,matrix(NA,np,4),n,NA,use_bw)
    }
    out
  }) # 10k pts 11.7s slow; 10s fast; 9.7s fast single
  t2 <- unclass(Sys.time())
  if(any(ss3<-sapply(Out,inherits,"error"))) stop(Out[[which(ss3)[1]]])
  data$loo.fit <- data$full.fit <- NA

  if(fast){
    out <- t(abind::abind(Out,along=2))
    yh1 <- if(np==1) X[idx,]*out[,1:np] else apply(X[idx,]*out[,1:np],1,sum)
    yh2 <- if(np==1) X[idx,]*out[,(np+1):(2*np)] else apply(X[idx,]*out[,(np+1):(2*np)],1,sum)
    # coefs <- abind::abind(loo=out[,1:np],full=out[,(np+1):(2*np)],along=3)
    data$loo.fit[idx] <- yh1
    data$full.fit[idx] <- yh2
    bw_used <- out[,"bw"]
  } else {
    data$loo.fit[idx] <- pmax(sapply(Out,function(x) x[[1]]),0)
    coefs <- aperm(abind::abind(lapply(Out,function(x) x[[2]]),along=3),c(3,1,2))
    bw_used <- sapply(Out,function(x) x[[5]])
  }
  data$loo.err <- pmax(data$loo.fit,0)-data[,rn]
  tt <- round(c(t.prep=t1-t0,t.calc=t2-t1,t.stats=unclass(Sys.time())-t2,t.tot=unclass(Sys.time())-t0),3)
  loostats <- c(rsq=rsq(pmax(data[idx,rn],0),data$loo.fit[idx]),mae=mean(abs(data$loo.err[idx]),na.rm=T),t=unclass(Sys.time())-t0,n.na=sum(is.na(data$loo.err[idx])),avg_bw=mean(bw_used),tt)
  if(verbose) {cat("done\n");print(loostats)}

  result <- new("gwr", form=form, name=name, rdata=data, loostats=loostats, bw=bw, adaptive=adaptive, bw_used=bw_used, fit_idx=idx, buffer=buffer)
  return(result)
}



