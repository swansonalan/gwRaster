
#' plot a raster layer
#'
#' Makes a set of plots including:
#' map of fits
#' scatterplot of fits vs response
#' map of errors
#' qqnorm plot
#' variogram
#'
#' @param x a RasterLayer object
#' @return nothing
#' @export
setMethod("plot",signature(x="RasterLayer"),function(x,hill=NULL,sym=F,crop=1,zlab="",mar.fract=.15,...){
  v <- as.vector(sampleRegular(x,100000))
  if(sum(!is.na(v) & is.finite(v))<100) stop("no values in image")
  rng <-range(v,na.rm=T)
  prng <- if(crop<1) quantile(v,c((1-crop)/2,1-(1-crop)/2),na.rm=T) else rng
  if(sym) prng <- c(-1*max(abs(prng)),max(abs(prng)))

  brks <- seq(prng[1],prng[2],length.out=100)
  mids <- (brks[-1]+brks[-length(brks)])/2
  lmat <- matrix(mids,nrow=length(mids),ncol=5)
  cols <- if(sym) grDevices::colorRampPalette(c("blue","white","red"))(length(mids)) else fields::tim.colors(length(mids))
  # if(crop<1) brks[c(1,length(brks))]<-c(min(v,na.rm=T),max(v,na.rm=T))+c(-1,1)
  brks[c(1,length(brks))] <- c(rng[1]-.25*(rng[2]-rng[1]),rng[2]+.25*(rng[2]-rng[1]))

  if(!is.null(hill)){
    X <- col2rgb(cols)/255
    cols <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=.6)
  }

  # legend ~~~~
  op <- par(plt=c(1-mar.fract,1-.8*mar.fract,par()$plt[3:4]),mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9,new=F)
  image(list(z=t(lmat),y=mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=brks,col=cols,...)
  axis(side=4,at=pretty(mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
  mtext(zlab,side=4,line=1.5,xpd=NA,cex=.7)

  # actual plot ~~~~~
  par(plt=c(op$plt[1],1-1.05*mar.fract,par()$plt[3:4]),new=T)
  if(!is.null(hill)){
      raster::image(hill,col=gray(seq(0,1,length=500)),xaxt="n",yaxt="n",xlab="",ylab="",...)
      raster::image(x,breaks=brks,col=cols,add=T)
  } else {
    image(x,breaks=brks,col=cols,xaxt="n",yaxt="n",xlab="",ylab="")
  }
  if(isLonLat(x)) fields::US(add=T)
})

#' plot a raster layer
#'
#' Makes a set of plots including:
#' map of fits
#' scatterplot of fits vs response
#' map of errors
#' qqnorm plot
#' variogram
#'
#' @param x a RasterLayer object
#' @return nothing
#' @export
setMethod("plot",signature(x="SpatRaster"),function(x,hill=NULL,sym=F,crop=1,zlab="",mar.fract=.15,...){
  v <- as.vector(terra::spatSample(x,100000,"regular")[,1])
  if(sum(!is.na(v) & is.finite(v))<100) stop("no values in image")
  rng <- range(v,na.rm=T)
  prng <- if(crop<1) quantile(v,c((1-crop)/2,1-(1-crop)/2),na.rm=T) else rng
  if(sym) prng <- c(-1*max(abs(prng)),max(abs(prng)))
  if(sym) rng <- c(-1*max(abs(rng)),max(abs(rng)))

  brks <- seq(prng[1],prng[2],length.out=100)
  mids <- (brks[-1]+brks[-length(brks)])/2
  lmat <- matrix(mids,nrow=length(mids),ncol=5)
  cols <- if(sym) grDevices::colorRampPalette(c("blue","white","red"))(length(mids)) else fields::tim.colors(length(mids))
  # if(crop<1) brks[c(1,length(brks))]<-c(min(v,na.rm=T),max(v,na.rm=T))+c(-1,1)
  brks[c(1,length(brks))] <- c(rng[1]-.25*(rng[2]-rng[1]),rng[2]+.25*(rng[2]-rng[1]))

  if(!is.null(hill)){
    X <- col2rgb(cols)/255
    cols <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=.5)
  }

  # legend ~~~~
  op <- par(plt=c(1-mar.fract,1-.8*mar.fract,par()$plt[3:4]),mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9,new=F)
  terra::image(list(z=t(lmat),y=mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=brks,col=cols,...)
  axis(side=4,at=pretty(mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
  mtext(zlab,side=4,line=1.5,xpd=NA,cex=.7)

  # actual plot ~~~~~
  par(plt=c(op$plt[1],1-1.05*mar.fract,par()$plt[3:4]),new=T)
  if(!is.null(hill)){
    terra::image(hill,col=gray(seq(0,1,length=500)),xaxt="n",yaxt="n",xlab="",ylab="",...)
    terra::image(x,breaks=brks,col=cols,add=T)
  } else {
    terra::image(x,breaks=brks,col=cols,xaxt="n",yaxt="n",xlab="",ylab="")
  }
  if(isLonLat(x)) fields::US(add=T)
})

#' plot a SpatRaster layer
#'
#' Makes a set of plots including:
#' map of fits
#' scatterplot of fits vs response
#' map of errors
#' qqnorm plot
#' variogram
#'
#' @param x a RasterLayer object
#' @return nothing
#' @export
rplot <- function(x,hill=NULL,sym=F,crop=1,zlab="",mar.fract=.18,brks=NULL,...){
  if(is.null(brks)){
    v <- as.vector(terra::spatSample(x,100000,"regular")[,1])
    if(sum(!is.na(v) & is.finite(v))<100) stop("no values in image")
    brks <-  brkfun(v,crop=crop,sym=sym)
    } else brks <- brkfun(brks,crop=1,sym=sym)

  # legend ~~~~
  op <- par(plt=c(1-mar.fract,1-.8*mar.fract,par()$plt[3:4]),mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9,new=F)
  terra::image(list(z=t(brks$lmat),y=brks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=brks$brks,col=brks$cols)
  axis(side=4,at=pretty(brks$mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
  mtext(zlab,side=4,line=1.4,xpd=NA,cex=.8)

  # actual plot ~~~~~
  par(plt=c(op$plt[1],1-1.05*mar.fract,par()$plt[3:4]),new=T)
  if(!is.null(hill)){
    terra::image(hill,col=gray(seq(0,1,length=500)),xaxt="n",yaxt="n",xlab="",ylab="",...)
    terra::image(x,breaks=brks$brks,col=brks$alphacols,add=T)
  } else {
    terra::image(x,breaks=brks$brks,col=brks$cols,xaxt="n",yaxt="n",xlab="",ylab="",...)
  }
  if(isLonLat(x)) fields::US(add=T)
}


#' Coefficient plots for a GWR fit.
#'
#' Makes a set of plots for each cofficient of a gwr model
#'
#' @param x a gwr.coef.fit object created using gwr_coef_map()
#' @return nothing
#' @export
setMethod("plot",signature(x="gwr.coef.fit"),function(x,hill=NULL,crop=.99,title=NULL){
  print(class(x))
  img <- terra::rast(x@coef.name)
  vn <- c("intercept",all.vars(as.formula(x@form))[-1])
  np <-  ncol(model.matrix(as.formula(x@form),x@rdata[1:3,]))
  par(mfrow=c(ceiling(np/2),2),mar=c(1,1,1,1),oma=c(0,0,3,0))
  for(i in 1:np) rplot(terra::rast(x@coef.name,lyrs=i),hill=hill,crop=crop,sym=T,zlab=vn[i])
  if(is.null(title)) title <- paste0(x@form,", bandwidth=",round(x@loostats[5],1),"km")
  mtext(title,side=3,outer=T,line=1.25)
})

#' Plots for a GWR raster fit.
#'
#' Makes a set of plots for each cofficient of a gwr model
#'
#' @param x a gwr.coef.fit object created using gwr_coef_map()
#' @return nothing
#' @export
setMethod("plot",signature(x="gwr.raster.fit"),function(x,title=NULL,mar.fract=.18,show.errs=F,zlab=NULL,...){
  if(is.null(title)) title <- x@name
  rn <- strsplit(x@form,"~")[[1]][1]
  img <- terra::rast(x@fit.name)
  rbrks <- brkfun(as.vector(terra::spatSample(img,100000,"regular")[,1]),crop=.95)
  tmp <- x@rdata[!is.na(x@rdata$rast.fit),]
  tmp$err <- tmp$rast.fit-tmp[,rn]
  if(nrow(tmp)>5000) tmp <- tmp[sample(nrow(tmp),5000),]
  ebrks <- brkfun(tmp$err,sym=T,crop=.98)

  #dpar(mar=c(1,1,3,1))
  if(is.null(zlab)) zlab <- paste("predicted",rn)
  op <- par(plt=c(1-mar.fract*.95,1-.8*mar.fract,par()$plt[3:4]),mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9,new=F)
  graphics::image(list(z=t(rbrks$lmat),y=rbrks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=rbrks$brks,col=rbrks$cols)
  axis(side=4,at=pretty(rbrks$mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
  mtext(zlab,side=4,line=1.4,xpd=NA,cex=.8)

  if(show.errs){
    par(plt=c(.8*mar.fract,.95*mar.fract,par()$plt[3:4]),mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9,new=T)
    graphics::image(list(z=t(ebrks$lmat),y=ebrks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=ebrks$brks,col=ebrks$cols)
    axis(side=2,at=pretty(ebrks$mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
    mtext("error",side=2,line=1.5,xpd=NA,cex=.9)
    tmp <- project(terra::vect(tmp,crs=geog),crs(img))
    par(plt=c(mar.fract,1-mar.fract,par()$plt[3:4]),mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9,new=T)
  } else {
    par(plt=c(op$plt[1],1-mar.fract,par()$plt[3:4]),mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9,new=T)

    }

  if("hill" %in% names(x@img.fn)){
    hill <- terra::rast(x@img.fn["hill"])
    terra::image(hill,col=gray(seq(0,1,length=500)),xaxt="n",yaxt="n",xlab="",ylab="",main=title,...)
    terra::image(img,breaks=rbrks$brks,col=rbrks$alphacol,add=T)
  } else {
    terra::image(img,breaks=brks,col=cols,xaxt="n",yaxt="n",xlab="",ylab="")
  }
  if(show.errs) terra::points(tmp,ebrks$col.vec,pch=16,cex=.9)
})

#' Diagnostic plots of a LOOCV GWR fit.
#'
#' Makes a set of plots including:
#' map of fits
#' scatterplot of fits vs response
#' map of errors
#' qqnorm plot
#' variogram
#'
#' @param x a gwr object created using gwr_loocv()
#' @return nothing
#' @export
setMethod("plot",signature(x="gwr"),function(x,title=NULL,max.pts=7500){

  # data ~~~~~~~~~~~~~~~~
  tmp <- x@rdata[!is.na(x@rdata$loo.fit),]
  XY <- as.matrix(crds(project(vect(tmp,crs=geog),my_aea)))/1000
  tmp$x <- XY[,1];tmp$y <- XY[,2]
  tmp$loo.fit <- pmax(tmp$loo.fit,0)
  tmp$full.fit <- pmax(tmp$full.fit,0)
  ss <- if(nrow(x@rdata)<max.pts) 1:nrow(rdata) else sample(nrow(x@rdata),max.pts)
  rn <- strsplit(x@form,"~")[[1]][1]
  tmp$err <- tmp$full.fit - tmp[,rn]
  tmp$loo.err <- tmp$loo.fit - tmp[,rn]

  # variogram
  ss.vg <- if(nrow(tmp) > max.pts) sample(nrow(tmp),max.pts) else 1:nrow(tmp)
  md <- sqrt((max(tmp$x)-min(tmp$x))^2+(max(tmp$y)-min(tmp$y))^2)/2
  p=3
  vbrks <- c(0:10)^p/10^p*md
  vg0 <- geoR::variog(data=tmp$loo.err[ss.vg],coords=XY[ss.vg,],breaks=vbrks,max.dist=2000,messages=F,plot=F)
  tmp2 <- data.frame(v=vg0$v,u=vg0$u,n=vg0$n)

  # plotting colors ~~~
  err.brks <- brkfun(tmp$loo.err,sym=T,crop=.99)
  fit.brks <- brkfun(tmp$full.fit)

  # actual plot ~~~~~
  op <- par(mfrow=c(3,2),mar=c(5,4,3,1),oma=c(0,0,3,0))
  par(plt=c(.08,.85,.08,.88),new=F,mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9)
  plot(lat~lon,data=tmp,col=fit.brks$col.vec,pch=16,xlab="",ylab="")
  fields::US(add=T)
  par(plt=c(.87,.90,.08,.88),new=T)
  image(list(z=t(fit.brks$lmat),y=fit.brks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=fit.brks$brks,col=fit.brks$cols)
  axis(side=4,at=pretty(fit.brks$mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
  mtext("fit",side=4,line=1.5,xpd=NA,cex=.7)

  par(mar=c(5,4,3,1))
  xyr <- range(c(tmp$loo.fit,tmp[,rn]),na.rm=T)
  plot(tmp[ss.vg,rn],tmp$loo.fit[ss.vg],xlim=xyr,ylim=xyr,main=paste0("R^2=",round(x@loostats[1],3)),ylab="loocv fit",xlab=rn)
  abline(0,1,lty=2,col=2,lwd=2)

  par(plt=c(.08,.85,.08,.88),new=F,mgp=c(2,.75,0),cex.axis=.8,cex.lab=.9)
  plot(lat~lon,data=tmp,col=err.brks$col.vec,pch=16,xlab="",ylab="")
  fields::US(add=T)
  par(plt=c(.87,.90,.08,.88),new=T)
  image(list(z=t(err.brks$lmat),y=err.brks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=err.brks$brks,col=err.brks$cols)
  axis(side=4,at=pretty(err.brks$mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
  mtext("error",side=4,line=1.5,xpd=NA,cex=.7)

  par(mar=c(5,4,3,1))
  qqnorm(tmp$loo.err[ss.vg])
  qqline(tmp$loo.err[ss.vg])

  par(mar=c(5,4,3,1))
  plot(v~u,data=tmp2,type="b",ylab="semi-variance",xlab="distance (km)",main="Variogram",ylim=c(0,max(tmp2$v)))

  # plot(loo.fit~full.fit,data=tmp[ss.vg,]);abline(0,1,lty=2,col=2,lwd=2)
  qqnorm(tmp$loo.fit-tmp$full.fit,main="loo effect")
  if(is.null(title)) title <- paste0(x@form,", bandwidth=",round(x@loostats[5],1),"km")
  mtext(title,side=3,outer=T,line=1.25)
  par(op)
})

#' Utility function for creating breakpoints and colors for plotting
#'
#' Makes a set of plots including:
#' map of fits
#' scatterplot of fits vs response
#' map of errors
#' qqnorm plot
#' variogram
#'
#' @param v a numeric vector of values
#' @param crop sets the breaks and colors to cover this quantile of values.  useful if one or two outliers dominate the colorscale
#' @param sym logical.  if true, a blue/white/red colorscheme is used, with breaks centered on zero
#' @return a list with breaks, mids, colors, and matrix used for plotting the colorscale
#' @export
brkfun <- function(v,crop=1,sym=F,pad=T){
  if(sum(!is.na(v))==0){
    rng <- c(-1,1)
    crop <- 1
  } else rng <- range(v,na.rm=T)
  if(rng[1]==rng[2]) {rng <- c(rng[1]-1,rng[2]+1);crop<-1}
  prng <- if(crop<1) quantile(v,c((1-crop)/2,1-(1-crop)/2),na.rm=T) else rng
  if(sym) prng <- c(-1*max(abs(prng)),max(abs(prng)))
  if(sym) rng <- c(-1*max(abs(rng)),max(abs(rng)))

  brks <- seq(prng[1],prng[2],length.out=100)
  mids <- (brks[-1]+brks[-length(brks)])/2
  lmat <- matrix(mids,nrow=length(mids),ncol=5)
  cols <- if(sym) grDevices::colorRampPalette(c("blue","white","red"))(length(mids)) else fields::tim.colors(length(mids))
  # if(crop<1) brks[c(1,length(brks))]<-c(min(v,na.rm=T),max(v,na.rm=T))+c(-1,1)
  if(pad) brks[c(1,length(brks))] <- c(rng[1]-.25*(rng[2]-rng[1]),rng[2]+.25*(rng[2]-rng[1]))

  X <- col2rgb(cols)/255
  alphacols <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=.6)

  col.vec <- cols[cut(v,brks)]
  return(list(brks=brks,mids=mids,cols=cols,lmat=lmat,col.vec=col.vec,alphacols=alphacols))
}
