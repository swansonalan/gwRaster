
#' Makes GWR coefficient estimates on a grid
#'
#' This function fits a gwr model to a data.frame which includes geographic information
#'
#'
#' @param object a gwr.coef.fit object, typically created using gwr_coef_map()
#' @param img.fn filename of a raster to use as grid, or a raster object
#' @param verbose logical; whether or not to report progress
#' @param overwrite logical; whether or not to overwrite existing raster fit
#' @return a gwr object with the original inputs, original data, fit, statistics.
#' @export
gwr_fit_map <- function(x,img.fn,verbose=T,overwrite=F,overwrite.resamps=F){
  t0 <- unclass(Sys.time())
  if(!inherits(x,"gwr.coef.fit")) stop("x needs to be a gwr.coef.fit object")
  if(!test_cluster(F)) stop("need cluster")
  if(foreach::getDoParWorkers()==1)  stop("need cluster")

  rn <- strsplit(x@form,"~")[[1]][1]
  vn <- all.vars(as.formula(x@form))[-1]
  ss <- complete.cases(x@rdata[,c(rn,vn)])
  np <-  ncol( model.matrix(as.formula(x@form),x@rdata[1:3,]))
  if(!all(vn %in% names(img.fn))) stop(paste("missing image files for",paste(vn[!vn %in% names(img.fn)],collapse=", "),"for",x@form))
  img.fn2 <- img.fn[match(vn,names(img.fn))]
  if(inherits(try(rast(img.fn2)),"try-error")) stop("one or more raster inputs appear to be of the wrong extent/projection, or are missing, or are corrupt")

  coef.name <- x@coef.name
  out.dir <- paste0(dirname(img.fn2[1]),"/gwr_fits/")
  if(!file.exists(out.dir)) dir.create(out.dir)
  out.name <- paste0(out.dir,x@name,"_preds.tif")
  tmpnames <- paste0(out.dir,c("int",tolower(vn)),"_",x@name,"_coef_30m.tif")
  if(!file.exists(coef.name)) stop(paste(coef.name,"doesn't exist"))
  if(file.exists(out.name) & overwrite) file.remove(out.name)
  if(!file.exists(out.name)){
      # resample coefficients ~~~~~~~~~~~~~~~~~~~~~~~~~
      if(overwrite.resamps & sum(file.exists(tmpnames))>0)  file.remove(tmpnames[file.exists(tmpnames)])
      ss <- file.exists(tmpnames)
      if(sum(!ss)==0){
        if(verbose) cat(paste(basename(tmpnames),collapse=", "),"already made, re-using\n")
      } else {
        if(verbose & sum(ss)>0) cat(paste(basename(tmpnames[ss]),collapse=", "),"already exists, re-using\n")
        # coefs <- rast(coef.name)
        # tmplt <- rast(img.fn2[1])
        if(verbose) cat("resampling GWR coefficients to 30m grid\n")
        tr_args <- sapply((1:np)[!ss],function(i) get_proj_string(coef.name,img.fn2[1],band=i))
        if(test_cluster()){
          foreach::foreach(arg=paste(tr_args,tmpnames[!ss])) %dopar% system(arg)
        } else {
          foreach::foreach(arg=paste(tr_args,tmpnames[!ss])) %do% system(arg)
        }
        t1 <- unclass(Sys.time())

        # if(test_cluster()){
        #   # r <- rast(coef.name)
        #   # coefs <- wrap(coefs)
        #   # if(img.prj!=coef.prj) foreach::foreach(i=(1:np)[!ss],.packages='terra') %dopar% terra::project(terra::rast(coefs,lyrs=i),tmplt,"bilinear",filename=tmpnames[i],overwrite=T)
        #   # if(img.prj==coef.prj) foreach::foreach(i=(1:np)[!ss],.packages='terra') %dopar% terra::resample(terra::rast(coefs,lyrs=i),tmplt,"bilinear",filename=tmpnames[i],overwrite=T)
        #   if(img.prj!=coef.prj) foreach::foreach(i=(1:np)[!ss],.packages='terra') %dopar% terra::project(terra::rast(coef.name,lyrs=i),tmplt,"bilinear",filename=tmpnames[i],overwrite=T)
        #   if(img.prj==coef.prj) foreach::foreach(i=(1:np)[!ss],.packages='terra') %dopar% terra::resample(terra::rast(coef.name,lyrs=i),tmplt,"bilinear",filename=tmpnames[i],overwrite=T)
        # } else {
          # if(img.prj!=coef.prj) for(i in (1:np)[!ss]) terra::project(coefs[[i]],tmplt,"bilinear",filename=tmpnames[i],overwrite=T)
          # if(img.prj==coef.prj) for(i in (1:np)[!ss]) terra::resample(coefs[[i]],tmplt,"bilinear",filename=tmpnames[i],overwrite=T)
        # }
      }
      # gdal_calc to make final 30m raster ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      fn <- c(tmpnames,img.fn2)
      lets <- LETTERS[-8][1:length(fn)]
      bn <- paste(paste0("-",lets," ",fn),collapse=" ")
      form2 <- paste0("'maximum(A+",paste(paste0(lets[2:(np)],"*",lets[(np+1):(2*np-1)]),collapse="+"),",0)'")
      sys.command <- paste0("gdal_calc.py ",bn," --outfile=",out.name," --calc=",form2," --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
      if(verbose) cat("calculating final raster using",form2,"\n")
      system(sys.command) # 18min
      file.remove(tmpnames)
      cat("done\n\n")
  } else {
      if(verbose) cat(out.name,"already made, re-using\n")
      t1 <- unclass(Sys.time())
    }
  t2 <- unclass(Sys.time())
  x@rdata$rast.fit <- textract(terra::rast(out.name), x@rdata)
  rast.stats <- c(rsq=rsq(x@rdata[,rn],x@rdata$rast.fit),mae=mean(abs(x@rdata$rast.fit-x@rdata[,rn]),na.rm=T),mae_rast_vs_pt=mean(abs(x@rdata$full.fit-x@rdata$rast.fit),na.rm=T))
  tt <- round(c(t.resamp=t1-t0,t.calc=t2-t1,t.stats=unclass(Sys.time())-t2,t.tot=unclass(Sys.time())-t0),1)
  if(verbose) {cat("done\n");print(c(rast.stats,tt))}

  out <- new("gwr.raster.fit",x,img.fn=img.fn,fit.name=out.name,rast.stats=rast.stats)
  return(out)
}






