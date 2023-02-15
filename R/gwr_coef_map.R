#' Makes GWR coefficient estimates on a grid
#'
#' This function fits a gwr model to a data.frame which includes geographic information
#'
#'
#' @param object a gwr object, typically created using gwr_loocv()
#' @param tmplt.name filename of a raster to use as grid, or a raster object
#' @param verbose logical; whether or not to report progress
#' @return a gwr object with the original inputs, original data, fit, statistics.
#' @export
gwr_coef_map <- function(x,tmplt.name,name="",out.dir=NULL,buffer=3,verbose=T,overwrite=F,fast=T){
  t0 <- unclass(Sys.time())
  if(!test_cluster(F)) stop("need cluster")
  if(foreach::getDoParWorkers()==1)  stop("need cluster")
  if(is.null(out.dir)) out.dir <- paste0(dirname(tmplt.name),"/gwr_coefs/")
  if(!file.exists(out.dir)) dir.create(out.dir)
  coef.name <- paste0(out.dir,x@name,"_coefs.tif")
  '%dopar%' <- foreach::'%dopar%'
  '%do%' <- foreach::'%do%'

  # masked.coef.name <- paste0(x@out.dir,"coefs/",x@name,"_coefs_masked.tif")
  if(file.exists(coef.name) & overwrite) file.remove(coef.name)
  if(!file.exists(coef.name)){
    rn <- strsplit(x@form,"~")[[1]][1]
    vn <- all.vars(as.formula(x@form))[-1]
    ss <- complete.cases(x@rdata[,c(rn,vn)])
    if(fast){
      X <- model.matrix(as.formula(x@form),model.frame(x@rdata))
      np <- ncol(X)
      nc <- np+2
    } else {
      np <-  ncol( model.matrix(as.formula(x@form),x@rdata[1:3,]))
      nc <- np*4+2
    }
    Tmplt <- terra::rast(tmplt.name)
    tmplt <- cbind(x=terra::xFromCell(Tmplt,1:terra::ncell(Tmplt)),y=terra::yFromCell(Tmplt,1:terra::ncell(Tmplt)),dem=terra::values(Tmplt))
    out <- matrix(NA,nrow(tmplt),nc)
    length(idx <- which(!is.na(tmplt[,3])))
    if(verbose)  cat("estimating GWR coefficients on a grid\n")
    XY <- as.matrix(terra::crds(terra::project(terra::vect(x@rdata,crs=geog),crs(Tmplt,T))))
    x@rdata$x <- XY[,1];x@rdata$y <- XY[,2]

    if(verbose) cat(paste0("Calculating GWR coefficients on a ",nrow(Tmplt)," by ",ncol(Tmplt)," grid (",length(idx)," pts) using ",foreach::getDoParWorkers()," cores\n"))
    Out <- foreach(i=1:length(idx),.errorhandling='pass') %dopar% {
      xy <- tmplt[idx[i],]
      d <- sqrt((XY[,1]-xy[1])^2+(XY[,2]-xy[2])^2)/1000
      use_bw <- if(x@adaptive) sort(d)[x@bw] else x@bw
      md <- if(x@adaptive) sort(d)[x@bw*x@buffer] else x@bw*x@buffer

      if(fast){
        ss2 <- d<md & ss
        z <- if(np==1) try(lm.wfit(x=matrix(X[ss2,],ncol=1),y=x@rdata[ss2,rn],w=dnorm(d[ss2]/use_bw))$coefficients) else  try(lm.wfit(x=X[ss2,],y=x@rdata[ss2,rn],w=dnorm(d[ss2]/use_bw))$coefficients)
        y <- if(inherits(z,"try-error")) rep(NA,np) else z
        y <- if(length(y)<np) rep(NA,np) else y
        n <- sum(ss2)
      } else {
        tmp <- cbind(x@rdata[,c(rn,vn)],d=d)[d<md & ss,]
        tmp$wts <- dnorm(tmp$d/use_bw)
        if(length(unique(tmp[,rn]))>(np+1)){
          y <- summary(lm(x@form,weights=wts,data=tmp))$coefficients
          if(nrow(y)<np) y <- matrix(NA,np,4)
        } else y <- matrix(NA,np,4)
        n <- nrow(tmp)
      }
      c(as.vector(y),use_bw,n)
    }
    if(any(ss<-sapply(Out,inherits,"error"))) stop(Out[[which(ss)[1]]])
    if(verbose) cat("done, writing raster\n")
    out[idx,] <- t(abind::abind(Out,along=2))
    if(!fast) colnames(out) <- c(paste0(rep(c("int",tolower(vn)),4),"_",rep(c("est","se","t","p"),each=np)),"bwd","n")
    if(fast) colnames(out) <- c("int",tolower(vn),"bwd","n")
    Tmplt <- terra::rast(Tmplt,nlyrs=nc,vals=F)
    values(Tmplt) <- out
    writeRaster(Tmplt,coef.name,overwrite=T,NAflag=-9999,gdal=c("COMPRESS=LZW"))

    if(verbose) cat("done, time used",round(unclass(Sys.time())-t0,1),"seconds\n\n")
    out <- new("gwr.coef.fit",x,tmplt.name=tmplt.name,coef.name=coef.name,band.names=colnames(out))
  } else {
    if(verbose) cat(coef.name,"already exists, re-using\n")
    out <- new("gwr.coef.fit",x,tmplt.name="",coef.name=coef.name,band.names="")
  }
  return(out)
}





