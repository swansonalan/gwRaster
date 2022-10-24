test.cluster <- function(verbose=T) {
  require(foreach);require(doParallel)
  x <- !inherits(nowarnings(try(foreach(i=1:3,.combine='c') %dopar% i)),"try-error")
  y <- getDoParRegistered()
  n <- getDoParWorkers()
  z <- x & y & n>1
  if(z) {if(verbose) cat("found",n,"core cluster\n")} else {if(verbose) cat("no valid cluster found\n")}
  z
}


gwr.cv <- function(form,abw=NULL,bw=NULL,idx=NULL,rdata,n.samp=1000,verbose=T){
  t0 <- unclass(Sys.time())
  require(foreach);require(abind);require(geodist);require(doParallel)
  if((is.null(abw) & is.null(bw)) | (!is.null(abw) & !is.null(bw))) stop("need to specify one of abw or bw")
  if(is.null(idx)) {idx <- if(nrow(rdata)>n.samp) sample(nrow(rdata),n.samp) else 1:nrow(rdata)}
  if(!test.cluster(F)) stop("need cluster")
  if(getDoParWorkers()==1)  stop("need cluster")
  fixbw <- !is.null(bw)
  # vn <- sapply(strsplit(strsplit(form,"~")[[1]][2],"\\+")[[1]],"\\*")[[1]]
  rn <- strsplit(form,"~")[[1]][1]
  vn <- all.vars(as.formula(form))
  ss <- complete.cases(rdata[,vn])
  np <-  ncol( model.matrix(as.formula(form),rdata[1:3,]))
  
  # cat("making LOO predictions for",form,"\n")
  system.time(Out <- foreach(i=idx,.packages=c('geodist'),.errorhandling='pass') %dopar% {
    rdata$d <- as.vector(geodist(rdata[,c("lon","lat")],rdata[i,c("lon","lat")],measure="haversine"))/1000
    if(!fixbw) bw <- sort(rdata$d)[abw+1]
    rdata$wts <- dnorm(rdata$d/bw)
    ss2 <- rdata$d<bw*3 & rdata$d>0 & ss
    if((n<-sum(ss2))>10){
      nu <- if(n<100) qr(rdata[ss2,vn])$rank else 9999
      if(nu>=np){
        fit <- predict(m0 <- lm(form,weights=wts,data=rdata[ss2,]),newdata=rdata[i,])
        out <- list(fit,summary(m0)$coefficients,n)
      } else out <- list(NA,matrix(NA,np,4),n,nu)
    } else out <- list(NA,matrix(NA,np,4),n,NA)
  })
  
  if(any(ss<-sapply(Out,inherits,"error"))) stop(Out[[which(ss)[1]]])
  rdata$loo.fit <- NA
  rdata$loo.fit[idx] <- pmax(sapply(Out,function(x) x[[1]]),0)
  rdata$loo.err <- rdata$loo.fit-rdata[,rn]
  coefs <- aperm(abind(lapply(Out,function(x) x[[2]]),along=3),c(3,1,2))
  loostats <- c(rsq=rsq(rdata[idx,rn],rdata$loo.fit[idx]),mae=mean(abs(rdata$loo.err[idx]),na.rm=T),t=unclass(Sys.time())-t0,n.na=sum(is.na(rdata$loo.err[idx])))
  if(verbose) {print(loostats);cat("\n")}
  return(list(fit=rdata$loo.fit,loostats,coefs))
}

make_gwr_30m <- function(form,rdata,name,img.fn,mask.name,out.dir,abw=NULL,bw=NULL,nmax=2500,do.loo=T,buff=3,idx=NULL,overwrite=F){
  # form is what you think it is
  # rdata is a data.frame with 
  
  t0 <- unclass(Sys.time())
  require(foreach);require(abind);require(doParallel);require(raster);require(rgdal)
  
  server <- strsplit(Sys.info()[[4]],"\\.")[[1]][1]
  tmpdir <- switch(server,vulcan="/mnt/SSD2TB/data/tmp/",bucephalus="/mnt/nvmeDrive/data/tmp/")
  
  if((is.null(abw) & is.null(bw)) | (!is.null(abw) & !is.null(bw))) stop("need to specify one of abw or bw")
  fixbw <- !is.null(bw)
  
  if(!test.cluster(F)) stop("need valid cluster")
  vn <- all.vars(as.formula(form))[-1]
  rn <- strsplit(form,"~")[[1]][1]
  np <-  ncol(model.matrix(as.formula(form),rdata[1:3,]))
  if(!all(vn %in% names(img.fn))) stop("image names don't match vars in formula")
  img.fn <- img.fn[match(vn,names(img.fn))];names(img.fn) <- vn
  
  # filenames ~~~~
  coef.name <- paste0(out.dir,"coefs/",name,"_coefs.tif")
  masked.coef.name <- paste0(out.dir,"coefs/",name,"_coefs_masked.tif")
  out.name <- paste0(out.dir,"fits/",name,"_preds.tif")
  dat.name <- paste0(out.dir,"tabdata/",name,"_loo_fits.Rdata")
  tmpnames <- paste0(out.dir,"coefs/",c("int",tolower(vn)),"_",name,"_coef_30m.tif")
  all.fn <- c(coef.name,masked.coef.name,out.name,tmpnames)
  if(overwrite & any(file.exists(all.fn))) cat("removing",sum(file.remove(all.fn[file.exists(all.fn)])),"previous files\n")
  
  if(do.loo){
      cat("making LOO predictions\n")
      ss <- complete.cases(rdata[,c(rn,vn)])
      if(is.null(idx)) idx <- if(sum(ss)<nmax) which(ss) else sample(which(ss),nmax)
      system.time(Out <- foreach(i=idx,.packages=c('geodist'),.errorhandling='pass') %dopar% {
        rdata$d <- as.vector(geodist(rdata[,c("lon","lat")],rdata[i,c("lon","lat")],measure="haversine"))/1000
        if(!fixbw) bw <- sort(rdata$d)[abw+1]
        md <- if(!fixbw) sort(rdata$d)[abw*buff] else bw*buff
        tmp <- rdata[rdata$d<md & ss,]
        tmp$wts <- dnorm(tmp$d/bw)
        fit <- if(length(unique(tmp[,rn]))>(np+1)) predict(lm(form,weights=wts,data=tmp[tmp$d>0,]),newdata=rdata[i,]) else NA
        fit
      })
      
      if(any(ss<-sapply(Out,inherits,"error"))) stop(Out[[which(ss)[1]]])
      rdata$loo.fit <- NA
      rdata$loo.fit[idx] <- pmax(unlist(Out),0)
      rdata$loo.err <- rdata$loo.fit-rdata[,rn]
    
      loostats <- c(rsq=rsq(rdata[,rn],rdata$loo.fit),mae=mean(abs(rdata$loo.err),na.rm=T),n=length(idx))
      print(loostats)
      cat("\n")
  } else loostats <- NA
  
  #################################################################################################
  # estimate GWR coef on coarse grid ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #################################################################################################
  t1 <- unclass(Sys.time())
  nc <- np*4+2
  Tmplt <- raster(mask.name)
  tmplt <- cbind(x=xFromCell(Tmplt,1:ncell(Tmplt)),y=yFromCell(Tmplt,1:ncell(Tmplt)),dem=values(Tmplt))
  out <- matrix(NA,nrow(tmplt),nc)
  length(idx <- which(!is.na(tmplt[,3])))
  if(!file.exists(coef.name)){
      cat("estimating GWR coefficients on a grid\n")
      ss <- complete.cases(rdata[,c(rn,vn)])
    
      system.time(Out <- foreach(i=1:length(idx),.errorhandling='pass') %dopar% {
        xy <- tmplt[idx[i],]
        rdata$d <- sqrt((rdata$x-xy[1])^2+(rdata$y-xy[2])^2)/1000
        if(!fixbw) bw <- sort(rdata$d)[abw]
        md <- if(!fixbw) sort(rdata$d)[abw*buff] else bw*buff
        tmp <- rdata[rdata$d<md & ss,]
        tmp$wts <- dnorm(tmp$d/bw)
        if(length(unique(tmp[,rn]))>(np+1)){ 
            out <- summary(lm(form,weights=wts,data=tmp),newdata=rdata[i,])$coefficients
            if(nrow(out)<np) out <- matrix(NA,np,4)
        } else out <- matrix(NA,np,4)
        c(as.vector(out),bw,nrow(tmp))
      })
      if(any(ss<-sapply(Out,inherits,"error"))) stop(Out[[which(ss)[1]]])
      out[idx,] <- t(abind(Out,along=2))
      colnames(out) <- c(paste0(rep(c("int",tolower(vn)),4),"_",rep(c("est","se","t","p"),each=np)),"bwd","n")
      Tmplt <- brick(Tmplt,nl=nc,values=F)
      values(Tmplt) <- out
      writeRaster(Tmplt,coef.name,overwrite=T,NAflag=-9999,options=c("COMPRESS=LZW","TIFW=NO"))
      cat("done\n\n")
  }
  
  # plot coefs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(!file.exists(masked.coef.name)){
      cat("making masked coef raster\n")
      img <- brick(coef.name)
      Cimg <- stack(foreach(i=2:np,.packages='raster') %dopar% {
          cimg <- raster(img,layer=i)
          pimg <- raster(img,layer=(np*3)+i)
          x <- values(cimg)
          x[values(pimg)>.01] <- NA
          values(cimg) <- x
          cimg
         })
      writeRaster(Cimg,masked.coef.name,overwrite=T,NAflag=-9999,options=c("COMPRESS=LZW","TIFW=NO"))
      cat("done\n\n")
  }

  # if(any(ss <- file.exists(tmpnames))) file.remove(tmpnames[ss])
  # make_gwr_raster(form,coef.name,img.fn,out.name,tmpnames,overwrite=F)
  
  # resample coefs to 30m ~~~~~~~~~~~~~~~~~~~~~~~~~
  t2 <- unclass(Sys.time())
  img <- raster(img.fn[1])
  res <- mean(res(img))
  ext <- as.vector(extent(img))
  dims <- ceiling(round(c((ext[2]-ext[1])/res,(ext[4]-ext[3])/res),3))
  proj <- paste0("'",projection(img),"'")
  
  if(any(ss <- file.exists(tmpnames))) file.remove(tmpnames[ss])
  cat("resampling GWR coefficients to 30m grid\n")
  system.time(Out2 <- foreach(i=1:np,.errorhandling='pass') %dopar% {
    sys.command <- paste("gdal_translate -q -a_nodata -9999 -r bilinear -b",i,"-projwin",ext[1],ext[4],ext[2],ext[3],"-outsize",paste(dims,collapse=" "),"-a_srs",proj,"-co 'COMPRESS=LZW' -co 'BIGTIFF=YES'",coef.name,tmpnames[i])
    system.time(system(sys.command)) # 3.7s
    1
  }) # 5.2min
  cat("done\n\n")
  
  # gdal_calc to make final 30m raster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  t3 <- unclass(Sys.time())
  fn <- c(tmpnames,img.fn)
  lets <- LETTERS[-8][1:length(fn)]
  bn <- paste(paste0("-",lets," ",fn),collapse=" ")
  form <- paste0("'maximum(A+",paste(paste0(lets[2:(np)],"*",lets[(np+1):(2*np-1)]),collapse="+"),",0)'")
  # form <- paste0("'maximum(A+",paste(paste0(lets[2:(np+1)],"*",lets[(np+2):(2*np)]),collapse="+"),",0)'")
  
  sys.command <- paste0("gdal_calc.py ",bn," --outfile=",out.name," --calc=",form," --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
  cat("calculating final raster using",form,"\n")
  system.time(system(sys.command)) # 18min
  cat("done\n\n")
  
  rdata$rast.fit <- raster::extract(raster(out.name),as.sp(rdata))
  rast.stats <- c(rsq=rsq(rdata[,rn],rdata$rast.fit),mae=mean(abs(rdata$rast.fit-rdata[,rn]),na.rm=T),mae_rast_vs_loo=mean(abs(rdata$loo.fit-rdata$rast.fit),na.rm=T))
  save(rdata,loostats,rast.stats,file=dat.name)
  cat(paste0("done, time used: ",round(( unclass(Sys.time())-t0)/60,1)," min\n"));cat("raster stats:\n")
  print(rast.stats)
  cat("\n\n")

  # file.remove(tmpnames)
  invisible(rdata)
}

make_gwr_raster <- function(form,coef.name,img.fn,out.name,tmpnames=NULL,overwrite=F){
  t0 <- unclass(Sys.time())
  vn <- all.vars(as.formula(form))[-1]
  if(any(!vn%in%names(img.fn))) {print(img.fn);stop(paste("img.fn names don't include variables in ",form))}
  np <- length(vn)+1
  img.fn <- img.fn[match(vn,names(img.fn))]
  
  if(is.null(tmpnames)) tmpnames <- paste0(sub(".tif","",out.name),"_",c("int",vn),"_30m_resamp.tif")
  all.fn <- c(tmpnames,out.name)
  if(overwrite & any(file.exists(all.fn))) cat("removing",sum(file.remove(all.fn[file.exists(all.fn)])),"previous files\n")
  if(!overwrite & any(file.exists(all.fn))) stop("file(s) exist, use overwrite=T")
  
  img <- raster(img.fn[1])
  res <- mean(res(img))
  ext <- as.vector(extent(img))
  dims <- ceiling(round(c((ext[2]-ext[1])/res,(ext[4]-ext[3])/res),3))
  proj <- paste0("'",projection(img),"'")
  
  cat("resampling GWR coefficients to 30m grid\n")
  system.time(Out2 <- foreach(i=1:(length(vn)+1),.errorhandling='pass') %dopar% {
    sys.command <- paste("gdal_translate -q -a_nodata -9999 -r bilinear -b",i,"-projwin",ext[1],ext[4],ext[2],ext[3],"-outsize",paste(dims,collapse=" "),"-a_srs",proj,"-co 'COMPRESS=LZW' -co 'BIGTIFF=YES'",coef.name,tmpnames[i])
    system.time(system(sys.command)) # 3.7s
    1
  }) # 5.2min
  cat("done\n\n")
  
  # gdal_calc to make final 30m raster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  t3 <- unclass(Sys.time())
  fn <- c(tmpnames,img.fn)
  lets <- LETTERS[-8][1:length(fn)]
  bn <- paste(paste0("-",lets," ",fn),collapse=" ")
  form <- paste0("'maximum(A+",paste(paste0(lets[2:(np)],"*",lets[(np+1):(2*np-1)]),collapse="+"),",0)'")
  # form <- paste0("'maximum(A+",paste(paste0(lets[2:(np+1)],"*",lets[(np+2):(2*np)]),collapse="+"),",0)'")
  
  sys.command <- paste0("gdal_calc.py ",bn," --outfile=",out.name," --calc=",form," --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
  cat("calculating final raster using",form,"\n")
  system.time(system(sys.command)) # 18min
  cat("done\n\n")
  
  # gdal_calc to make final 30m raster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  t3 <- unclass(Sys.time())
  fn <- c(tmpnames,img.fn)
  lets <- LETTERS[-8][1:length(fn)]
  bn <- paste(paste0("-",lets," ",fn),collapse=" ")
  form <- paste0("'maximum(A+",paste(paste0(lets[2:(np)],"*",lets[(np+1):(2*np-1)]),collapse="+"),",0)'")
  # form <- paste0("'maximum(A+",paste(paste0(lets[2:(np+1)],"*",lets[(np+2):(2*np)]),collapse="+"),",0)'")
  
  sys.command <- paste0("gdal_calc.py ",bn," --outfile=",out.name," --calc=",form," --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
  cat("calculating final raster using",form,"\n")
  system.time(system(sys.command)) # 18min
  cat("done\n\n")
  invisible(c(tmpnames,out.name))
}

plotfun <- function(form,name,out.dir,subsets=NULL,bnd=NULL,plot.name=NULL){
  #png(paste0(out.dir,"figs/",name,"_broot_fits2.png"),11,12,units="in",res=150)
  vn <- all.vars(as.formula(form))[-1]
  coef.name <- paste0(out.dir,"coefs/",name,"_coefs.tif")
  masked.coef.name <- paste0(out.dir,"coefs/",name,"_coefs_masked.tif")
  out.name <- paste0(out.dir,"fits/",name,"_preds.tif")
  tmpnames <- paste0(out.dir,"coefs/",c("int",tolower(vn)),"_",name,"_coef_30m.tif")
  dat.name <- paste0(out.dir,"tabdata/",name,"_loo_fits.Rdata")
  load(dat.name)
  rdata$bwd <- raster::extract(raster(coef.name,band=(length(vn)+1)*4+1),as.sp(rdata))
  
  ext <- extent(raster(out.name))
  tmp <- rdata[rdata$x>ext[1] & rdata$x<ext[2] & rdata$y>ext[3] & rdata$y<ext[4],]
  
  if(!is.null(plot.name)) png(plot.name,length(vn)*3,ifelse(length(vn)<4,4,3)*3,units="in",res=150)
  dpar(mfrow=c(ifelse(length(vn)<4,4,3),length(vn)),oma=c(0,0,3,0))
  # dpar(mfrow=c(3,pmax(length(vn),4)),oma=c(0,0,3,0))
  
  for(i in 1:length(vn)) mplot(raster(masked.coef.name,band=i),bnd=bnd,hill=hill,Sym=T,plot.names=vn[i],mar.fract=.2)
  for(i in 1:length(vn)) mplot(raster(tmpnames[i+1]),bnd=bnd,hill=hill_br,Sym=T,plot.names=vn[i],mar.fract=.2)
  if(!is.null(subsets)){
    for(i in 1:length(subsets)) mplot(crop(raster(out.name),extent(subsets[[i]])),bnd=bnd,hill=subsets[[i]],Sym=F,plot.names="fit subset")
  }
  dgplot(tmp[,c("dem30","gsr30","loo.err")],zlab="",title="loo error",Sym=T,crop=.95)
  mplot(tmp[,c("x","y","loo.err")],plot.names="loo error",Sym=T,crop=.95,hill=hill_br)
  mtext(paste0(form,", bandwidth=",round(mean(rdata$bwd,na.rm=T)),"km, loo R-squared=",round(loostats[1],3)),side=3,outer=T)
  if(!is.null(plot.name)) dev.off()
}

dgplot <- function(mdat,zlab=NULL,xlab=NULL,ylab=NULL,Sym=F,rev.cols=F,title=NULL,crop=1,close.last=T,close.first=T,col="tim",pch=21,mar.fract=.15,title.fract=.07){
  library(raster);library(fields)
  if(close.first) try(close.screen(all=TRUE))
  pt.cex <- if(nrow(mdat)>1500) .6 else 1.1
  if(is.null(xlab)) xlab <- colnames(mdat)[1]
  if(is.null(ylab)) ylab <- colnames(mdat)[2]
  if(is.null(zlab)) zlab <- colnames(mdat)[3]
  brks <- bfun(mdat[,3],sym=Sym,crop=crop,n=500)
  if(!Sym) brks$cols <- brks$tim.cols
  if(rev.cols) brks$cols <- rev(brks$cols)
  if(is.null(title)) par(plt=c(.12,1-mar.fract,.12,.97)) else par(plt=c(.12,1-mar.fract,.12,.88))
  plot(mdat[,1],mdat[,2],cex=1,xlab=xlab,ylab=ylab,type="n")
  if(pch==21) points(mdat[,1],mdat[,2],pch=21,bg=brks$cols[cut(mdat[,3],brks$brks)],cex=pt.cex) else points(mdat[,1],mdat[,2],pch=16,col=brks$cols[cut(mdat[,3],brks$brks)],cex=pt.cex)
  if(!is.null(title)) title(main=title,xpd=NA,line=0.4) # if(!is.null(main)) title(main=main,xpd=NA,line=0.5,cex.main=.8) else 
  if(is.null(title)) par(plt=c(1-mar.fract+.02,1-mar.fract+.07,.03,.97)) else par(plt=c(1-mar.fract+.02,1-mar.fract+.07,.03,.9))
  par(new=TRUE,cex.axis=.8,cex.lab=.8)
  image(list(z=t(brks$lmat),y=brks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=brks$brks,col=brks$cols)
  axis(side=4,at=pretty(brks$mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
  mtext(zlab,side=4,line=1.5,xpd=NA,cex=.9)
  if(close.last) close.screen(all=TRUE) else {par(plt=c(.03,.85,.03,.97))}
  }

mplot <- function(mdat,hill=NULL,zlab="",bnd=F,Sym=F,brks=NULL,rev.cols=F,common.scale=F,plot.names=NULL,makefile=F,
                  layout=NULL,title=NULL,crop=1,close.last=T,close.first=T,col="tim",n=NULL,alph=NULL,pts=NULL,pch=21,
                  mar.fract=.15,title.fract=.07,logscale=F){
  # plots points with colors based on values
  # main input is 'mdat', which should be a data.frame with the first two columns as coordinates, and the remaining columns the variables to plot.
  # 'makefile' is an optional plot name.  If set to T it generates a plot name automatically.  If F, it plots to the current graphics device.
  # 'title' is an overall title (for more than 1 variable)
  # 'plot.names' is a vector of names for the individual plots
  # 'layout' is a plots if more than one variable is plotted.  Function will figure this out autmatically if not supplied.
  # 'Sym' causes the color scheme to be blue-white-red centered on zero.
  # 'common.scale' causes each plot to have the same scale
  # 'zlab' adds a label to the scalebar
  # optional inputs include 'hill', for a raster hillshade to be used as the background, and 'bnd' which causes a US boundary to be plotted using 
  # the fields function 'US'.  The 'title' argument adds a plot title
  
  library(raster);library(fields)
  if(close.first) try(close.screen(all=TRUE))
  custom.brks <- !is.null(brks)
  # op <- par()
  # on.exit(close.screen(all=TRUE))
  raster_mode <- inherits(mdat,"RasterBrick") | inherits(mdat,"RasterStack") | inherits(mdat,"RasterLayer")
  if(is.null(alph)) alph <- if(is.null(hill)) 1 else .6
  if(custom.brks){
    brks <- bfun(brks,sym=Sym[1],alpha=alph,crop=crop,n=n)
    if(!Sym[1]) brks$cols <- brks$tim.cols
    if(col=="greens") brks$cols <- brks$greens
    if(rev.cols) brks$cols <- rev(brks$cols)
  }
  n.vars <- if(raster_mode) nlayers(mdat) else ncol(mdat)-2
  if(n.vars==1) common.scale <- F
  if(length(Sym)==1 & n.vars>1 & !custom.brks) Sym <- rep(Sym,n.vars)
  nc <- if(is.null(layout)) ceiling(sqrt(n.vars)) else layout[2]
  nr <- if(is.null(layout)) ceiling(n.vars/nc) else layout[1]
  if(is.character(makefile)){
    fname <- makefile
    makefile <- T
  } else if(makefile) fname <- paste0("/home/aswanson/tmp/rplot_",format(Sys.time(),"%Y%m%d_%H%m"),".png")
  if(makefile) png(fname,9,7,units="in",res=150)
  
  if(common.scale){
    s0 <- split.screen(rbind(c(0,1-mar.fract,0,1),c(1-mar.fract,1,0,1)))
    # screen(s0[1])
    if(!custom.brks) {
      if(logscale){
        
      } else {
        brks <- if(raster_mode) bfun(as.vector(sampleRegular(mdat,100000)),sym=Sym[1],alpha=alph,crop=crop,n=n) else bfun(as.vector(mdat[,3:ncol(mdat)]),sym=Sym[1],crop=crop,n=n)
        if(!Sym[1]) brks$cols <- brks$tim.cols
        if(col=="greens") brks$cols <- brks$greens
        if(col=="rwg") brks$cols <- brks$rwg
        if(rev.cols) brks$cols <- rev(brks$cols)
      }
    }
    screen(s0[2])
    # par(mar=rep(0,4))
    # plot(1~1)
    par(plt=c(0.02,.2,.05,.9))
    image(list(z=t(brks$lmat),y=brks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=brks$brks,col=brks$cols)
    axis(side=4,at=pretty(brks$mids),mgp=c(1.5,.35,0),cex.axis=1,tcl=-.1)
    mtext(zlab,side=4,line=1.4,xpd=NA,cex=1.1)
    screen(s0[1])
  }
  # print(n.vars)
  # if(is.null(plot.names)) plot.names <- if(raster_mode) names(mdat) else colnames(mdat)[3:ncol(mdat)] 
  if(!is.null(title)) {
    s2 <- split.screen(rbind(c(0,1,0,1-title.fract),c(0,1,1-title.fract,1)))
    screen(s2[2])
    mtext(title,side=3,outer=TRUE,line=-1.7,cex=1.25)
    screen(s2[1])
  }
  
  if(n.vars>1) s1 <- split.screen(c(nr,nc))
  pt.cex <- if(nrow(mdat)>1500) .6 else 1.1
  for(i in 1:n.vars){
    Zlab <- if(length(zlab)==1) zlab else zlab[i]
    if(raster_mode) mimg <- if(n.vars==1) mdat else raster(mdat,layer=i)
    if(n.vars>1) screen(s1[i])
    if(!common.scale & !custom.brks) {
      brks <- if(raster_mode){
        samp <- if(ncell(mimg)>30000) sampleRegular(mimg,30000) else values(mimg)
        bfun(samp,sym=Sym[i],alpha=alph,crop=crop,n=n) 
      } else bfun(mdat[,i+2],sym=Sym[i],crop=crop,n=n);
      if(!Sym[i]) brks$cols <- brks$tim.cols
      if(col=="greens") brks$cols <- brks$greens
      if(col=="rwg") brks$cols <- brks$rwg
      if(rev.cols) brks$cols <- rev(brks$cols)
    }
    if(!common.scale){if(is.null(plot.names)) par(plt=c(.03,1-mar.fract,.03,.97)) else par(plt=c(.03,1-mar.fract,.03,.88))}
    if(common.scale) {if(is.null(plot.names)) par(plt=c(.03,.97,.03,.97)) else par(plt=c(.03,.97,.03,.88))}
    if(is.null(hill)){
      if(raster_mode){
        image(mimg,breaks=brks$brks,col=brks$cols,xaxt="n",yaxt="n",xlab="",ylab="")
      } else {
        plot(mdat[,1],mdat[,2],cex=1,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
        if(pch==21) points(mdat[,1],mdat[,2],pch=21,bg=brks$cols[cut(mdat[,i+2],brks$brks)],cex=pt.cex) else points(mdat[,1],mdat[,2],pch=16,col=brks$cols[cut(mdat[,i+2],brks$brks)],cex=pt.cex)
      }
    } else {
      image(hill,col=gray(seq(0,1,length=500)),xaxt="n",yaxt="n",xlab="",ylab="")
      if(raster_mode){
        image(mimg,breaks=brks$brks,col=brks$cols,add=TRUE)
      } else 	{if(pch==21) points(mdat[,1],mdat[,2],pch=21,bg=brks$cols[cut(mdat[,i+2],brks$brks)],cex=pt.cex) else points(mdat[,1],mdat[,2],pch=16,col=brks$cols[cut(mdat[,i+2],brks$brks)],cex=pt.cex)}
    }
    if(!is.null(bnd)) if(is.logical(bnd))  US(add=TRUE) else plot(bnd,add=T)
    if(!is.null(pts)) points(pts[,1],pts[,2])
    if(!is.null(plot.names)) title(main=plot.names[i],xpd=NA,line=0.4,cex.main=ifelse(n.vars>9,.6,1)) # if(!is.null(main)) title(main=main,xpd=NA,line=0.5,cex.main=.8) else 
    if(!common.scale){
      if(is.null(plot.names)) par(plt=c(1-mar.fract+.02,1-mar.fract+.07,.03,.97)) else par(plt=c(1-mar.fract+.02,1-mar.fract+.07,.03,.9))
      par(new=TRUE,cex.axis=.8,cex.lab=.8)
      image(list(z=t(brks$lmat),y=brks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=brks$brks,col=brks$cols)
      axis(side=4,at=pretty(brks$mids),mgp=c(.3,.25,0),cex.axis=.8,tcl=-.1)
      mtext(Zlab,side=4,line=1.5,xpd=NA,cex=.9)
      
    }
  }
  # if(common.scale){
  #   screen(s0[2])
  #   # par(mar=rep(0,4))
  #   # plot(1~1)
  #   par(plt=c(0.02,.2,.05,.9))
  #   image(list(z=t(brks$lmat),y=brks$mids,x=1:5),xaxt="n",yaxt="n",ylab="",xlab="",breaks=brks$brks,col=brks$cols)
  #   axis(side=4,at=pretty(brks$mids),mgp=c(1.5,.35,0),cex.axis=1,tcl=-.1)
  #   mtext(Zlab,side=4,line=1.4,xpd=NA,cex=1.1)
  # }
  # if(!is.null(title)) {screen(s2[2]);mtext(title,side=3,outer=TRUE,line=-1.7,cex=1.25)}
  # if(n.vars>1 | !is.null(title)) close.screen(all=TRUE)
  if(close.last) close.screen(all=TRUE) else {par(plt=c(.03,.85,.03,.97))}
  if(makefile) dev.off()
  # if(makefile & sendfile) system(paste("~/sendmail.py",fname))
}
bfun <- function(x,sym=F,n=NULL,alpha=1,crop=1,log=F){
  require(fields);require(RColorBrewer)
  x<-x[!is.na(x) & is.finite(x)]
  mrng <-rng <- range(x,na.rm=TRUE) 
  if(crop<1) rng <- as.vector(quantile(x,c((1-crop)/2,1-(1-crop)/2)))
  if(crop>1) rng <- c(-crop,crop)
  if(rng[1]==rng[2]) rng <- rng+c(0,1)
  m <- 10^(round(log10(rng[2]-rng[1]))-2)
  rng <- c(m*floor(rng[1]/m),m*ceiling(rng[2]/m))
  if(rng[1]==min(x,na.rm=T)) rng[1]<-rng[1]-m
  if(rng[2]==max(x,na.rm=T)) rng[2]<-rng[2]+m
  if(sym==TRUE) {rng <- c(-1*max(abs(rng)),max(abs(rng)));mrng<-c(-1*max(abs(mrng)),max(abs(mrng)))}
  brks <- if(is.null(n)) seq(rng[1],rng[2],by=m) else seq(rng[1],rng[2],length=n)
  mids <- (brks[-1]+brks[-length(brks)])/2
  lmat <- matrix(mids,nrow=length(mids),ncol=5)
  if(crop==1){
    idx <- max(c(1:length(brks))[brks<min(x)]):min(c(1:length(brks))[brks>max(x)])
  } else idx <- 1:length(brks)
  tbrks <- brks[idx]
  tmids <- (tbrks[-1]+tbrks[-length(tbrks)])/2
  tlmat <- matrix(tmids,nrow=length(tmids),ncol=5)
  X<-col2rgb(colorRampPalette(c("blue","white","red"))(length(mids)))/255
  cols <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=alpha)
  
  X<-col2rgb(tim.colors(length(mids)))/255
  cols2<-rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=alpha)
  
  X<-col2rgb(colorRampPalette(c("black","white"))(length(mids)))/255
  cols3 <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=alpha)
  
  X<-col2rgb(colorRampPalette(c("white","dark green"))(length(mids)))/255
  cols4 <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=alpha)
  
  X<-col2rgb(colorRampPalette(c("dark orange","white","dark green"))(length(mids)))/255
  cols5 <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=alpha)
  
  # X<-col2rgb(brewer.pal(length(mids),"RdYlGn"))/255
  # cols6 <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=alpha)
  
  X<-col2rgb(colorRampPalette(c("#A50026","#FFFFBF","#006837"))(length(mids)))/255
  cols6 <- rgb(red=X[1,],green=X[2,],blue=X[3,],alpha=alpha)
  
  col.vec <- cols[cut(x,brks)]
  tim.col.vec <- cols2[cut(x,brks)]
  green.col.vec <- cols4[cut(x,brks)]
  rwg.col.vec <- cols5[cut(x,brks)]
  if(crop!=1) brks[c(1,length(brks))]<-mrng+c(-m/100,m/100)
  return(list(brks=brks,mids=mids,lmat=lmat,tbrks=tbrks,tmids=tmids,tlmat=tlmat,cols=cols,tim.cols=cols2,col.vec=col.vec,tim.col.vec=tim.col.vec,bw.cols=cols3,
              greens=cols4,rwg=cols6,rwg.col.vec=rwg.col.vec,green.col.vec=green.col.vec))
}














