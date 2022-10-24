
#' Makes GWR coefficient estimates on a grid
#'
#' This function is for pan-sharpening a more accurate coarse grid using finer, but possibly biased, estimate of the same thing.
#' Steps:
#' 1. resample the coarse image to the fine resolution
#' 2. resample the fine image to the coarse resolution
#' 3. resample the coarse version of the fine back to fine res.
#' 4. compute ratio of the original fine imag to it's smoothed version (3)
#' 5. multiply (1) by (4)
#'
#'
#' @param img.coarse a filename for a coarse-res prediction
#' @param img.fine filename for a fine-res prediction
#' @param verbose logical; whether or not to report progress
#' @param overwrite logical; whether or not to overwrite existing raster fit
#' @return a gwr object with the original inputs, original data, fit, statistics.
#' @export
pan_sharpen <- function(gwr_coarse,gwr_fine,smooth.fact=8,verbose=T,overwrite=F,overwrite.smooth=F,suffix="pansharpen",gdal=T){
  t0 <- unclass(Sys.time())
  if(!inherits(gwr_coarse,"gwr.coef.fit")) stop("gwr_coarse needs to be a gwr.coef.fit object")
  if(!inherits(gwr_fine,"gwr.coef.fit")) stop("gwr_fine needs to be a gwr.coef.fit object")
  if(!test_cluster(F)) stop("need cluster")
  if(foreach::getDoParWorkers()==1)  stop("need cluster")
  rn <- strsplit(gwr_coarse@form,"~")[[1]][1]
  vn <- all.vars(as.formula(gwr_coarse@form))[-1]

  gwr_pansharpen <- gwr_coarse
  gwr_pansharpen@name <- paste0(gwr_coarse@name,"_pansharpened",smooth.fact,"x_with_",gwr_fine@name)


  gwr_pansharpen@fit.name <- paste0(dirname(gwr_coarse@fit.name),"/",gwr_pansharpen@name,"_preds.tif")
  gwr_pansharpen@fine.smooth.name <- sub(".tif",paste0("_smooth",smooth.fact,"x.tif"),gwr_fine@fit.name)

  if(file.exists(gwr_pansharpen@fit.name) & overwrite) file.remove(gwr_pansharpen@fit.name)
  if(file.exists(gwr_pansharpen@fine.smooth.name) & overwrite.smooth) file.remove(gwr_pansharpen@fine.smooth.name)
  if(!file.exists(gwr_pansharpen@fit.name)){
    fine.prj <- terra::crs(terra::rast(gwr_fine@fit.name),T)
    coarse.prj <- terra::crs(terra::rast(gwr_coarse@fit.name),T)
    # img.fine.smoothed <- terra::disagg(terra::aggregate(terra::rast(img.fine),smooth.fact),smooth.fact,"bilinear")
    # img.ps <- rast(img.coarse)*rast(img.fine)/img.fine.smoothed
    if(gdal){
      if(!file.exists(gwr_pansharpen@fine.smooth.name)){
        if(verbose) cat("making smoothed version of finescale\n")
        tmpname <- tempfile(fileext=".tif")
        img <- terra::rast(gwr_fine@fit.name)
        res <- mean(res(img))
        Ext <- as.vector(ext(img))
        dims <- c(terra::nrow(img),terra::ncol(img))
        # dims <- ceiling(round(c((Ext[2]-Ext[1])/res,(Ext[4]-Ext[3])/res),3))
        coarse_res <- res*smooth.fact
        coarse_dims <-  ceiling(c((Ext[2]-Ext[1])/res,(Ext[4]-Ext[3])/coarse_res))
        proj <- paste0("'",crs(img,T),"'")

        if(file.exists(tmpname)) file.remove(tmpname)
        cat("resampling fine gwr to coarse grid\n")
        sys.command <- paste("gdal_translate -q -a_nodata -9999 -r average -projwin",Ext[1],Ext[4],Ext[2],Ext[3],"-outsize",paste(coarse_dims,collapse=" "),"-a_srs",proj,"-co 'COMPRESS=LZW' -co 'BIGTIFF=YES'",gwr_fine@fit.name,tmpname)
        system(sys.command) # 66s

        cat("resampling back to fine grid\n")
        sys.command <- paste("gdal_translate -q -a_nodata -9999 -r bilinear -projwin",Ext[1],Ext[4],Ext[2],Ext[3],"-outsize",paste(dims[2:1],collapse=" "),"-a_srs",proj,"-co 'COMPRESS=LZW' -co 'BIGTIFF=YES'",tmpname,gwr_pansharpen@fine.smooth.name)
        system(sys.command) # 208s
        file.remove(tmpname)
      }
      t1 <- unclass(Sys.time())
      fn <- c(gwr_coarse@fit.name,gwr_fine@fit.name,gwr_pansharpen@fine.smooth.name)
      lets <- LETTERS[-8][1:length(fn)]
      bn <- paste(paste0("-",lets," ",fn),collapse=" ")
      # form <- paste0("'maximum(A*B/C,0)'")
      form <- "'numpy.where(C==0, 0, maximum(A*B/C,0))'"
      sys.command <- paste0("gdal_calc.py ",bn," --outfile=",gwr_pansharpen@fit.name," --calc=",form," --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
      cat("calculating final raster using",form,"\n")
      system(sys.command) # 6min
      cat("done\n\n")
    } else {
      if(!file.exists(gwr_pansharpen@fine.smooth.name)){
          if(verbose) cat("making smoothed version of finescale\n")
          if(fine.prj!=coarse.prj) img.fine.smoothed <- terra::project(terra::aggregate(terra::rast(gwr_fine@fit.name),smooth.fact),rast(gwr_coarse@fit.name),"bilinear",
                                                                       filename=gwr_pansharpen@fine.smooth.name,gdal=c("COMPRESS=LZW","BIGTIFF=YES","TFW=NO"))
          if(fine.prj==coarse.prj) img.fine.smoothed <- terra::resample(terra::aggregate(terra::rast(gwr_fine@fit.name),smooth.fact),rast(gwr_coarse@fit.name),"bilinear",
                                                                        filename=gwr_pansharpen@fine.smooth.name,gdal=c("COMPRESS=LZW","BIGTIFF=YES","TFW=NO"))
          t1 <- unclass(Sys.time())
        }
      if(verbose) cat("performing main calc\n")
      writeRaster(rast(gwr_coarse@fit.name)*rast(gwr_fine@fit.name)/rast(gwr_pansharpen@fine.smooth.name),filename=gwr_pansharpen@fit.name,gdal=c("COMPRESS=LZW","BIGTIFF=YES","TFW=NO"))
    }
  } else t1 <- unclass(Sys.time())
  if(verbose) cat("calculating stats\n")
  t2 <- unclass(Sys.time())
  gwr_pansharpen@rdata$fine.fit <- gwr_fine@rdata$rast.fit
  gwr_pansharpen@rdata$fine.smoothed.fit <- textract(rast(gwr_pansharpen@fine.smooth.name),gwr_pansharpen@rdata)
  gwr_pansharpen@rdata$coarse.fit <- gwr_coarse@rdata$rast.fit
  gwr_pansharpen@rdata$ps.fit <- textract(rast(gwr_pansharpen@fit.name),gwr_pansharpen@rdata)

  gwr_pansharpen@rast.stats <- c(rsq=rsq(gwr_pansharpen@rdata[,rn],gwr_pansharpen@rdata$ps.fit),mae=mean(abs(gwr_pansharpen@rdata$ps.fit-gwr_pansharpen@rdata[,rn]),na.rm=T),
                                 mae_rast_vs_loo=mean(abs(gwr_pansharpen@rdata$loo.fit-gwr_pansharpen@rdata$ps.fit),na.rm=T))
  tt <- round(c(t.resamp=t1-t0,t.calc=t2-t1,t.stats=unclass(Sys.time())-t2,t.tot=unclass(Sys.time())-t0),1)
  # if(verbose) cat(paste0("done, time used: ",tt[4]," seconds (",paste(tt[1:3],collapse="s;"),")\n\n"))
  if(verbose) {cat("done\n");print(c(gwr_pansharpen@rast.stats,tt))}
  return(gwr_pansharpen)
}












