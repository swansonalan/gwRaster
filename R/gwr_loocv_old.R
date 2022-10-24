


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
gwr_loocv_old <- function(form,rdata,name="",out.dir="./",bw=NaN,adaptive=F,idx=NULL,n.samp=3000,verbose=T,buffer=3){
  t0 <- unclass(Sys.time())
  # require(foreach);require(abind);require(geodist);require(doParallel)
  if(!all(c("lon","lat") %in% colnames(rdata))) stop("need 'lon' and 'lat' colnames in rdata")
  # if((is.na(abw) & is.na(bw)) | (!is.na(abw) & !is.na(bw))) stop("need to specify one of abw or bw")
  if(is.null(idx)) {idx <- if(nrow(rdata)>n.samp) sample(nrow(rdata),n.samp) else 1:nrow(rdata)}
  if(!test_cluster(F)) stop("need cluster")
  if(foreach::getDoParWorkers()==1)  stop("need cluster")
  if(name=="") name <- paste0(gsub("\\+","_",sub("~","_vs_",form)),"_",ifelse(adaptive,"abw","bw"),bw)
  # fixbw <- !is.na(bw)
  # vn <- sapply(strsplit(strsplit(form,"~")[[1]][2],"\\+")[[1]],"\\*")[[1]]
  rn <- strsplit(form,"~")[[1]][1]
  vn <- all.vars(as.formula(form))
  ss <- complete.cases(rdata[,vn])
  np <-  ncol( model.matrix(as.formula(form),rdata[1:3,]))

  if(verbose) cat("making LOO predictions for",form,"using cluster with",foreach::getDoParWorkers(),"cores\n")
  system.time(Out <- foreach::foreach(i=idx,.packages=c('geodist'),.errorhandling='pass') %dopar% {
    rdata$d <- as.vector(geodist::geodist(rdata[,c("lon","lat")],rdata[i,c("lon","lat")],measure="haversine"))/1000
    use_bw <- if(adaptive) sort(rdata$d)[bw+1] else bw
    rdata$wts <- dnorm(rdata$d/use_bw)
    ss2 <- rdata$d<use_bw*buffer & rdata$d>0 & ss
    if((n<-sum(ss2))>10){
      nu <- if(n<100) qr(rdata[ss2,vn])$rank else 9999
      if(nu>=np){
        fit <- predict(m0 <- lm(form,weights=wts,data=rdata[ss2,]),newdata=rdata[i,])
        out <- list(fit,summary(m0)$coefficients,n,nu,use_bw)
      } else out <- list(NA,matrix(NA,np,4),n,nu,use_bw)
    } else out <- list(NA,matrix(NA,np,4),n,NA,use_bw)
    as.vector(out)
  })

  if(any(ss<-sapply(Out,inherits,"error"))) stop(Out[[which(ss)[1]]])
  rdata$loo.fit <- NA
  rdata$loo.fit[idx] <- pmax(sapply(Out,function(x) x[[1]]),0)
  rdata$loo.err <- rdata$loo.fit-rdata[,rn]
  coefs <- aperm(abind::abind(lapply(Out,function(x) x[[2]]),along=3),c(3,1,2))
  bw_used <- sapply(Out,function(x) x[[5]])
  loostats <- c(rsq=rsq(rdata[idx,rn],rdata$loo.fit[idx]),mae=mean(abs(rdata$loo.err[idx]),na.rm=T),t=unclass(Sys.time())-t0,n.na=sum(is.na(rdata$loo.err[idx])),avg_bw=mean(bw_used))
  if(verbose) {print(loostats);cat("\n")}

  result <- new("gwr", form=form, name=name, out.dir=out.dir, rdata=rdata, fit=rdata$loo.fit, coefs=coefs, loostats=loostats, bw=bw, adaptive=adaptive, bw_used=bw_used, fit_idx=idx, buffer=buffer)
  return(result)
}



