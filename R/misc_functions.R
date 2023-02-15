geog <- "+proj=longlat +datum=WGS84"
my_aea <- "+proj=aea +lat_1=33 +lat_2=60 +lat_0=33 +lon_0=-120 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
date.seq <- function(x,y) format(seq(as.Date(x,"%Y%m%d"),as.Date(y,"%Y%m%d"),by=1),"%Y%m%d")


#' Wrapper for file.info that produces a more compact output ~~~~~~~
#'
#' @param fnames a vector of filenames
#' @param sort whether or not to sort the table by modification date
#' @return a data.frame of file info including size in mb and modification date
#' @export
fi <- function(fnames,sort=T){
  x <- file.info(fnames)
  x$size <- round(x$size/2^20,2)
  rownames(x) <- basename(fnames)
  colnames(x)[1] <- "size (mb)"
  if(sort) x <- x[order(x$mtime,decreasing=T),1:4]
  x
}


#' Wrapper for terra::extract which automatically reprojects points to match the raster ~~~~~~~
#'
#' @param img a rast raster
#' @param locs a data.frame with columns named 'lat' and 'lon'
#' @param ... additional arguemnts to terra::extract e.g. method='bilinear'
#' @return numeric vector of values at points
#' @export
textract <- function(img,locs,...){
  img.prj <- terra::crs(img,T)
  if(!inherits(locs,"data.frame")) locs <- as.data.frame(locs)
  if("latitude" %in% names(locs)) names(locs)[names(locs)=="latitude"] <- "lat"
  if("longitude" %in% names(locs)) names(locs)[names(locs)=="longitude"] <- "lon"
  lcs <- if(terra::is.lonlat(img)) terra::vect(locs,crs=geog) else terra::project(terra::vect(locs,crs=geog),img.prj)
  out <- terra::extract(img,lcs,ID=F,...)
  if("ID" %in% colnames(out)) out <- out[,-match("ID",colnames(out))]
  if(dim(img)[3]==1) as.vector(as.matrix(out)[,1]) else out
}

#' lists objects in the current environment, and their size ~~~~~~~
#'
#' @param ... additional arguemnts to terra::extract e.g. method='bilinear'
#' @return numeric vector of values at points
#' @export
objs<-function(n=10,verbose=T,envir=.GlobalEnv){
  z <- gc()[,c(2,4,6)];colnames(z) <- paste(c("used","trigger","max"),"(Mb)")
  if(verbose) print(z)
  x <- ls(envir=envir)
  y <- sapply(x, function(z) as.numeric(sub(" Mb","",format(object.size(get(z,envir=envir)),units="MB"))))
  # x <- ls()
  # y <- sapply(x, function(z) as.numeric(sub(" Mb","",format(object.size(get(z)),units="MB"))))
  if(verbose) print(data.frame(size=y[order(y,decreasing=T)][1:pmin(n,length(y))]))
  invisible(y)
}


#' wrapper for list.files() with additional features ~~~~~~~
#' 'pattern' can be a vector
#' 'any'
#' 'strict'
#' 'none'
#'
#' @param ... additional arguemnts to terra::extract e.g. method='bilinear'
#' @param
#' @return numeric vector of values at points
#' @export
listfiles <- function(...,include.dirs=F,full=T){
  # this is an enhanced version of list.files() that adds 4 additional arguments:
  # 'pattern' now allows a vector argument.  The output files will have each of these.
  # the new 'none' argument should also be a character vector.  It removes files with any one of it's patterns.
  # the new 'any' argument should be a character vector.  It further restricts the output to files that contain any one of it's patterns.
  # the new 'strict' argument modifies all to enforce that there should only be one of each.  If there are more or less than 1 it returns an NA.
  args <- list(...)
  # formals <- formals()
  # args <- c(args,subset(as.list(formals),names(formals)!="..."))
  argnames <-names(args)
  if(!"full"%in%argnames) args <- c(args,full=T)
  if(!"include.dirs"%in%argnames) args <- c(args,include.dirs=F)
  if(any(c("pattern","any","none","strict")%in%argnames)){
    tmpargs <- args
    if("any" %in% argnames){
      do.multi <- T
      multi <- args[["any"]]
      tmpargs <- subset(tmpargs,names(tmpargs)!="any")
    } else do.multi <- F
    if("none" %in% argnames){
      none <- args[["none"]]
      tmpargs <- subset(tmpargs,names(tmpargs)!="none")
    }
    if("strict" %in% argnames){
      strict <- as.logical(args[["strict"]])
      if(is.na(strict)) strict <- args[["strict"]]
      tmpargs <- subset(tmpargs,names(tmpargs)!="strict")
    } else strict<-F
    if(length(pat <- list(...)[["pattern"]])>1){
      tmpargs[["pattern"]]<-pat[1]
      fn <- do.call(list.files,tmpargs)
      for(i in 2:length(pat)) fn <- fn[grep(pat[i],basename(fn))]
    } else  fn <- do.call(list.files,tmpargs)
    if("none" %in% argnames) for(i in 1:length(none)){fn <- if(length(grep(none[i],basename(fn)))>0) fn[-grep(none[i],basename(fn))] else fn}
    if(do.multi){
      if(is.logical(strict)){
        if(strict & ("any" %in% argnames)) {fn <- unlist(lapply(multi,function(x) {y<-fn[grep(x,basename(fn))];ifelse(length(y)==1,y,NA)}));names(fn)<-multi}
        if(!strict & ("any" %in% argnames)) fn <- unlist(lapply(multi,function(x) fn[grep(x,basename(fn))]))
      } else if(strict=="latest") {
        fn <- unlist(lapply(multi,function(x) {y<-fn[grep(x,basename(fn))];y[order(file.info(y)$mtime,decreasing=T)][1]}))
        names(fn)<-multi
      }
    }
    return(fn)
  } else do.call(list.files,args)
}

#' Writes out a system call to gdal_translate or gdalwarp ~~~~~~~
#'
#' @param img a rast raster
#' @param locs a data.frame with columns named 'lat' and 'lon'
#' @param ... additional arguemnts to terra::extract e.g. method='bilinear'
#' @return numeric vector of values at points
#' @export
get_proj_string <- function(img,tgt,band=1,method="bilinear"){
  if(inherits(img,"RasterLayer")) img <- terra::rast(img)
  input.fn <- if(inherits(img,"character")) img else terra::sources(img)
  if(inherits(img,"character")) img <- terra::rast(img)
  if(inherits(tgt,"character")) tgt <- terra::rast(tgt)
  dext <- as.vector(terra::ext(tgt))
  dproj <- terra::crs(tgt,T)
  ddims <- c(terra::nrow(tgt),terra::ncol(tgt))

  sext <- as.vector(terra::ext(img))
  sproj <- terra::crs(img,T)
  sext <- as.vector(terra::ext(img))
  sdims <-  c(terra::nrow(img),terra::ncol(img))

  if(all(sext==dext) & all(ddims==sdims)) proj.string <- ""

  # if(all(c(sext[1]<=dext[1],sext[4]>=dext[4],sext[2]>=dext[2],sext[3]<=dext[3])) & isLonLat(img)){
  if(sproj==dproj){
    # If new image has an extent greater than template and is unprojected, gdal_translate works and is much faster ~~~
    ts <- paste(ddims[2:1],collapse=" ")
    proj.string <- paste0("gdal_translate -b ",band," -projwin ",paste(dext[c(1,4,2,3)],collapse=" ")," -outsize ",ts," -a_srs '",dproj,"' -r ",method," -ot Float32 -co 'COMPRESS=LZW' -co 'BIGTIFF=YES' -q -a_nodata -9999 ")
  } else {
    if(band>1) stop("can't warp single band of a multiband image")
    te <- paste(dext[c(1,3,2,4)],collapse=" ")
    tr <- paste(terra::res(tgt)[2:1],collapse=" ")
    ts <- paste(c(terra::ncol(tgt),terra::nrow(tgt)),collapse=" ")
    proj.string <- paste0("gdalwarp -wm 1000 -multi -wo 'NUM_THREADS=1' -t_srs '",dproj,"' -te ",te," -ts ",ts," -r ",method," -ot Float32 -co 'COMPRESS=LZW' -co 'BIGTIFF=YES' -q -dstnodata -9999 ")
  }
  paste0(proj.string,input.fn)
}

#' reprojects a raster using gdal commands ~~~~~~~
#'
#' @param img a spatRaster object or the name of a raster file
#' @param tgt a spatRaster object, the name of a raster file, or a numeric vector defining an extent.
#' @param fn.out a filename for the output raster
#' @param execute logical, whether to execute the command or return it as a string
#' @return either a spatRaster object invisiblly, or a text string with a gdal system command
#' @export
gwarp <- function(img,tgt,fn.out,execute=T){
  if(inherits(img,"character")){
    input.fn <- img
  } else {
    input.fn <- terra::sources(img)
    if(input.fn=="") {
      input.fn <- tempfile(fileext=".tif")
      writeRaster(img,input.fn,gdal=c("COMPRESS=LZW"))
    }
  }
  if(inherits(tgt,"SpatRaster")){
      cmd <- paste(get_proj_string(input.fn,tgt),fn.out)
  }
  if(inherits(tgt,"numeric")){
      img <- terra::rast(input.fn)
      res <- mean(terra::res(img))
      dims <- ceiling(round(c((tgt[2]-tgt[1])/res,(tgt[4]-tgt[3])/res),3))
      proj <- paste0("'",terra::crs(img,T),"'")
      cmd <- paste("gdal_translate -projwin",tgt[1],tgt[4],tgt[2],tgt[3],"-outsize",paste(dims,collapse=" "),"-a_srs",proj,"-co 'COMPRESS=LZW' -co 'BIGTIFF=YES'",input.fn,fn.out)
  }
  if(execute){
    system(cmd)
    invisible(rast(fn.out))
  } else return(cmd)
}

#' crops a set of spatRaster images ~~~~~~~
#'
#' @param img a spatRaster object or the name of a raster file
#' @param ext a spatRaster extent, or a numeric vector describing an extent.
#' @param dir an existing directory.  a subdirectory based on the center of the extent will be created
#' @return a vector of filenames for the subset rasters.
#' @export

make_subset <- function(my.ext,fn,dir="/mnt/DataDrive1/data/wg/raster/30m_subsets/",ext.name=NULL){
  '%dopar%' <- foreach::'%dopar%'
  '%do%' <- foreach::'%do%'
  img.names <- names(fn)
  if(!inherits(my.ext,"numeric")) my.ext <- as.numeric(my.ext)
  is.geog <- abs(my.ext[1])< 360 & abs(my.ext[2])< 360 & abs(my.ext[3])<90 & abs(my.ext[4])<90
  ext.cent <- abs(c((my.ext[1]+my.ext[2])/2,(my.ext[3]+my.ext[4])/2))
  ext.mult <- if(is.geog) 1000 else 10
  if(is.null(ext.name)) ext.name <- paste0(ext.mult*round(ext.cent[1]/ext.mult),"_",ext.mult*round(ext.cent[2]/ext.mult))
  ext.dir <- paste0(dir,ext.name,"/")
  if(!file.exists(ext.dir)) dir.create(ext.dir)
  outnames <- paste0(ext.dir,sub(".tif",paste0("_",ext.name,".tif"),basename(fn)))
  names(outnames) <- img.names
  if(any(ss.needed <- !file.exists(outnames))){
      print(ss.needed)
      cmds <- foreach::foreach(fn.in=fn,fn.out=outnames,.combine='c') %do% gwarp(fn.in,my.ext,fn.out,execute=F)
      if(test_cluster()) foreach::foreach(cmd=cmds[ss.needed]) %dopar% system(cmd) else foreach::foreach(cmd=cmds[ss.needed]) %do% system(cmd)
  }
  if((length(i <- grep("dem",img.names))==1) & (!"hill" %in% img.names)){
      hill.name <- sub("dem","hill",outnames[i])
      if(!file.exists(hill.name)){
          hill <- hillshade(terra::rast(outnames[i]))
          terra::writeRaster(hill,hill.name,gdal="COMPRESS=LZW")
      }
      outnames <- c(outnames,hill.name)
      names(outnames)[length(outnames)] <- "hill"
  }
  return(outnames)
}


#' Calculates R-squared
#'
#' This function tests to see if a dummy function can be executed over a doParallel cluster using
#' 'foreach() %dopar% ...' syntax.  Returns T or F based on whether the foreach executed without
#' error.  If verbose=T, the function prints the number of cores detected.
#'
#' @param y response varaible
#' @param yh prediction
#' @return numeric R-squared value
#' @export
rsq <- function(y,yh) ifelse(sum(ss<-complete.cases(cbind(y,yh)))<10,NA,pmax(1-sum((y-yh)^2,na.rm=TRUE)/sum((y[ss]-mean(y[ss],na.rm=TRUE))^2,na.rm=TRUE),0))


#' Converts a dataframe to a spatial points dataframe.
#'
#' @param xy data.frame or matrix with columns named 'lat' and 'lon' or 'latitude' and 'longitude'
#' @return same data.frame as a SpatialPointsDataFrame
#' @export
as.sp <- function(xy){
  # converts a dataframe to a spatial points dataframe #
  # dataframe must have columns named 'lat' and 'lon' or 'latitude' and 'longitude' #
  # wgs84 CRS assumed
  # written by Alan Swanson, 8-20-19 #
  xy <- as.data.frame(xy)
  if(!any(c("lat","latitude") %in% names(xy))) stop("data needs columns named 'lat' or 'latitude'")
  sp::coordinates(xy) <- if("lat" %in% names(xy)) c("lon","lat") else c("longitude","latitude")
  sp::proj4string(xy) <- sp::CRS("+proj=longlat +datum=WGS84")
  xy
}

#' Searches code directory for certain text strings and/or filenames.
#'
#' @param text text string to search for
#' @param name text to search filenames for
#' @param dir directory to look in
#' @param n number of instances to show
#' @param sort how to sort instances.
#' @return if a 'text' argument is given, prints a list of instances sorted as specified.  If a 'name' arg is given, prints filenames sorted based on modification times
#' @export
search.code2 <- function(text=NULL,name=NULL,dir="~/code/",n=30,sort=c("newest","random","oldest","line","latest"),include="all",mtl=220){
  sort <- match.arg(sort)
  if(is.null(text) & is.null(name)) stop("need to supply arguments for 'text' or 'name'")
  if(!is.null(text)) {
    x <- system(paste0("grep -rn --exclude-dir=.Rproj.user/ '",path.expand(dir),"' -e '",text,"'"),intern=T,ignore.stderr=T)
    fn <- sapply(strsplit(x,"\\:"),function(x) x[1])
    if(!is.null(name)){
      for(i in 1:length(name)){
        x <- x[grep(name[i],basename(fn))]
        fn <- fn[grep(name[i],basename(fn))]
      }
    }
    n.per.file <- table(fn)

    # remove duplicated references within a given file
    deets <- sapply(strsplit(x,"\\:"),function(y) paste(y[c(3:length(y))],collapse=":"))
    duped <- duplicated(paste(basename(fn),deets))
    x <- x[!duped]
    fn <- fn[!duped]
    deets <- deets[!duped]

    # deets <- sapply(strsplit(x,"\\:"),function(x) paste(x[2:length(x)],collapse=":"))
    line.numbers <- sapply(strsplit(x,"\\:"),function(x) paste(x[2],collapse=":"))
    # deets2 <- sapply(strsplit(x,"\\:"),function(x) paste(x[3:length(x)],collapse=":"))
    short.fn <- sapply(strsplit(basename(fn),"\\."), function(x) paste(x[1:(length(x)-1)],collapse="."))
    fn.dates <-sapply(strsplit(short.fn,"_"),function(x) {y<-x[!is.na(suppressWarnings(as.numeric(x))) & nchar(x)==6];if(length(y)==1) y else NA})
    mod.dates <- format(file.info(fn)$mtime,"%Y-%m-%d")

    # remove duplicate lines within a file

    #~     head(dts2 <- sapply(strsplit(y,"\\("),function(x) strsplit(x[2],"\\)")[[1]][1]))
    #~     Dts2 <- format(as.Date(dts2),"%Y%m%d")

    # mtm <- format(file.info(fn)$mtime,"%Y-%m-%d")
    n.per.file <- as.vector(n.per.file)[match(fn,names(n.per.file))]
    y <- paste0(basename(fn)," (n=",n.per.file,")(",mod.dates,")",line.numbers,": ",deets)

    # first instance within a file ~~~
    if(include!="all"){
      y <- y[!duplicated(fn)]
      mod.dates <- mod.dates[!duplicated(fn)]
      line.numbers <- line.numbers[!duplicated(fn)]
      deets <- deets[!duplicated(fn)]
    }
    y.crop <- paste0(substr(y,1,mtl),"\n")
    y.crop2 <- y.crop[grep(text,y.crop)]

    if(sort=="random") y.sort <- if(length(y.crop)>n) sample(y.crop,n) else y.crop
    if(sort=="oldest") y.sort <- y.crop[order(mod.dates)][1:pmin(n,length(y.crop))]
    if(sort=="newest" | sort=="latest") y.sort <- y.crop[order(mod.dates,decreasing=T)][1:pmin(n,length(y.crop))]
    if(sort=="line") {
      y.sort <- y.crop[order(line.numbers)][1:pmin(n,length(y.crop))]
      deets <- deets[order(line.numbers)][1:pmin(n,length(y.crop))]
      y.sort <- y.sort[!duplicated(deets)]
    }
    cat("found",length(x),"instances, showing",length(y.sort),":\n")
    cat(y.sort,"\n")
    invisible(y)
  } else {
    fn <- listfiles(path=dir,pattern=name,any=c(".r$",".R$"),recursive=T)
    x <- data.frame(name=fn[order(file.info(fn)$mtime,decreasing=T)],mtime=file.info(fn[order(file.info(fn)$mtime,decreasing=T)])$mtime)
    print(x[1:min(c(length(fn),n)),])
    invisible(x)
  }
}

#' Extracts date strings in YYYYMMDD format from filenames. These must be separated using underscores.
#'
#' @param filenames a character vector, presumably filenames.
#' @param dateformat whether to attempt to convert these to R date format
#' @param n manually specifies which underscore-separated element to extract
#' @return a text string or vector of dates, depending on 'dateformat'
#' @export
fdates <- function(filenames,dateformat=F,n=NULL){
  # this function extracts dates from filenames that contain a date in YYYYMMDD format.
  # the date portion of the filename must be separated by underscores from the rest of the
  # filename.  the function searches each filename for numeric elements with 8 digits.
  # the filenane extension is automatically removed.
  #
  # the 'filenames' argument is any vector of filenames
  # 'dateformat' causes the returned dates to be in R date format
  # 'n' overrides the automatic search and instead returns the nth element separated by underscores.
  if(length(filenames)==0) {cat("filenames supplied have length 0\n");return(NULL)}

  # strip extension ~~~~~~~~~~~~~~~~~~~~
  filenames <- sapply(strsplit(basename(filenames),"\\."), function(x) paste(x[1:(length(x)-1)],collapse="."))
  if(is.null(n)){
    dts<-sapply(strsplit(filenames,"_"),function(x) {y<-x[!is.na(suppressWarnings(as.numeric(x))) & nchar(x)==8];if(length(y)==1) y else NA})
  } else  dts<-sapply(strsplit(basename(filenames),"_"),function(x) x[n])
  if(dateformat) dts <- as.Date(dts,"%Y%m%d")
  dts
}


#' Makes a hillshade raster from a dem raster
#'
#' @param dem a digital elevation model in raster format.
#' @return a hillshade in raster format.
#' @export
hillshade <- function(dem)  terra::shade(terra::terrain(dem,"slope",unit="radians"),terra::terrain(dem,"aspect",unit="radians"))

# hillshade <- function(dem){
#   slope <- raster::terrain(dem,opt="slope")
#   aspect <- raster::terrain(dem,opt="aspect")
#   return(raster::hillShade(slope,aspect))
# }

#' Reprojects a DEM using an external call to the GDAL function 'gdalwarp'
#'
#' @param img a raster object to be reprojected, or the filename for one.
#' @param target a raster object with the desired projection and extent.
#' @param outname a filename for the output.  Can be left blank.
#' @param method interpolation method
#' @param wm memory to allocate for the gdal
#' @param nt number of threads to use in gdal
#' @param ot datatype of output following gdal conventions
#' @param NAflag nodata value given to output
#' @param overwrite whether or not to overwrite an existing outname
#' @return invisibly returns the reprojected rasters, and/or writes it out if 'outname' is given.
#' @export
gdalProjRaster <- function(img,target,method="bilinear",filename=NULL,wm=1000,nt=1,ot="Float32",NAflag=-9999,overwrite=F){
  x<-tempfile()
  if(is.character(img)) fname<-img else raster::writeRaster(img,fname<-paste(tempfile(),".tif",sep=""))
  timg <- if(is.character(target)) raster::raster(target) else target
  if(is.null(filename)) filename <-paste(tempfile(),".tif",sep="")
  dims <- c(terra::nrow(img),terra::ncol(img))
  extt <- terra::ext(timg)
  te <- paste(as.vector(terra::ext(timg))[c(1,3,2,4)],collapse=" ")
  # te <- paste(c(attr(extt,"xmin"),attr(extt,"ymin"),attr(extt,"xmax"),attr(extt,"ymax")),collapse=" ")
  tr <- paste(terra::res(timg)[2:1],collapse=" ")
  ts <- paste(c(terra::ncol(timg),terra::nrow(timg)),collapse=" ")
  t_srs <- terra::crs(timg,T)
  ov <- if(overwrite) " -overwrite " else ""
  system(paste0("gdalwarp -wm ",wm," -multi -wo 'NUM_THREADS=",nt,"' -t_srs '",t_srs,"' -te ",te," -ts ",ts," -r ",method," -ot ",ot," -co 'COMPRESS=LZW' -co 'BIGTIFF=YES' -q -dstnodata ",NAflag,ov," ",fname," ",filename))
  invisible(terra::rast(filename))
}

#' Projects a geographic extent
#'
#' @param x an extent or a raster object from which to obtain an extent
#' @param to a new projection, or a raster object from which to obtain a projection
#' @param from the old projection.  defaults to the projection of x, if x is a raster object, or geographic
#' @return a list with the corner coordinates in the new projection, an extent object in the new projection that contains the original extent, and a numeric vector of the same.
#' @export
xExtent <- function(x,to,from=NULL){
  if(inherits(x,"RasterLayer") | inherits(x,"RasterBrick")) {from <- raster::projection(x);z<- raster::extent(x)} else z<-x
  if(inherits(z,"numeric")) z <-  raster::extent(z)
  if(inherits(to,"RasterLayer") | inherits(to,"RasterBrick")) to <-  raster::projection(to)
  y <- data.frame(x=c(attr(z,"xmin"),attr(z,"xmin"),attr(z,"xmax"),attr(z,"xmax")),
                  y=c(attr(z,"ymin"),attr(z,"ymax"),attr(z,"ymax"),attr(z,"ymin")))
  # if(buffer>0) y <- y+matrix(c(-1,-1,1,1,-1,1,1,-1),ncol=2)*buffer
  if(is.null(from)) from <- "+proj=longlat +datum=WGS84"
  sp::coordinates(y)<-c("x","y")
  sp::proj4string(y)<-sp::CRS(from)
  zz <- data.frame(sp::spTransform(y,sp::CRS(to)))
  nxt <- c(min(zz[,1]),max(zz[,1]),min(zz[,2]),max(zz[,2]))
  return(list(zz[c(1:4,1),],raster::extent(nxt),nxt))
}

#' Resets all graphical parameters to default values and allows changes
#'
#' @param ... new parameters to pass to par()
#' @return an invisible output from par() call
#' @export
dpar <- function(...) {
  par(list(xlog=FALSE,ylog=FALSE,adj=0.5,ann=TRUE,ask=FALSE,bg='white',bty='o',cex=1,cex.axis=1,cex.lab=1,cex.main=1.2,cex.sub=1,
           col='black',col.axis='black',col.lab='black',col.main='black',col.sub='black',crt=0,err=0,family='sans',fg='black',fig=c(0,1,0,1),fin=c(7.014,14.056),
           font=1,font.axis=1,font.lab=1,font.main=2,font.sub=1,lab=c(5,5,7),las=0,lend='round',lheight=1,ljoin='round',lmitre=10,lty='solid',lwd=1,
           mai=c(1.36,1.093,1.093,0.56),mar=c(5.1,4.1,4.1,2.1),mex=1,mfcol=c(1,1),mfg=c(1,1,1,1),mfrow=c(1,1),mgp=c(3,1,0),mkh=0.001,new=FALSE,oma=c(0,0,0,0),
           omd=c(0,1,0,1),omi=c(0,0,0,0),pch=1,pin=c(5.361,11.602),plt=c(0.156,0.92,0.097,0.922),ps=16,pty='m',smo=1,srt=0,tck=NA,tcl=-0.5,usr=c(0,1,0,1),
           xaxp=c(0,1,5),xaxs='r',xaxt='s',xpd=FALSE,yaxp=c(0,1,5),yaxs='r',yaxt='s',ylbias=0.2))
  invisible(par(...))
}
