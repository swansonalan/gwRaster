
#' Makes GWR coefficient estimates on a grid
#'
#' This function is for pan-sharpening a more accurate coarse grid using finer, but possibly biased, estimate of the same thing.
#' Based on ratio compnenent substitucion (RCS) pan-sharpening following
#' Mhangara, P., Mapurisa, W. and Mudau, N., 2020. Comparison of image fusion techniques using satellite pour l’Observation
#'   de la Terre (SPOT) 6 satellite imagery. Applied Sciences, 10(5), p.1881.
#' see 'https://www.r-bloggers.com/2015/02/pan-sharpening-using-r/'
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
pan_sharpen <- function(gwr_coarse,gwr_fine,smooth.fact=8,verbose=T,overwrite=F,overwrite.smooth=F,gdal=T,ps.strength=1,resamp.method="bilinear",do.stats=F){
  t0 <- unclass(Sys.time())
  if(!inherits(gwr_coarse,"gwr.coef.fit")) stop("gwr_coarse needs to be a gwr.coef.fit object")
  if(!inherits(gwr_fine,"gwr.coef.fit")) stop("gwr_fine needs to be a gwr.coef.fit object")
  # if(!test_cluster(F)) stop("need cluster")
  # if(foreach::getDoParWorkers()==1)  stop("need cluster")
  rn <- strsplit(gwr_coarse@form,"~")[[1]][1]
  vn <- all.vars(as.formula(gwr_coarse@form))[-1]
  '%dopar%' <- foreach::'%dopar%'
  '%do%' <- foreach::'%do%'


  gwr_pansharpen <- gwr_coarse
  if(file.exists(gwr_fine@ps.multiplier.name)){
      x <- strsplit(gwr_fine@name,"_")[[1]]
      fine.name <- paste(x[1:(grep("pansharpen",x)-1)],collapse="_")
      gwr_pansharpen@name <- paste0(gwr_coarse@name,sub("pansharpener","_pansharpened",x[grep("pansharpener",x)]),"_",x[length(x)],"_",fine.name,"_strength",round(ps.strength*100))
      gwr_pansharpen@fit.name <- paste0(dirname(gwr_coarse@fit.name),"/",gwr_pansharpen@name,"_preds.tif")
      gwr_pansharpen@fine.smooth.name <- gwr_fine@fine.smooth.name
      gwr_pansharpen@ps.multiplier.name <- gwr_fine@ps.multiplier.name
  } else {
      gwr_pansharpen@name <- paste0(gwr_coarse@name,"_pansharpened",smooth.fact,"x_",resamp.method,"_with_",gwr_fine@name,"_strength",round(ps.strength*100))
      gwr_pansharpen@fit.name <- paste0(dirname(gwr_coarse@fit.name),"/",gwr_pansharpen@name,"_preds.tif")
      gwr_pansharpen@fine.smooth.name <- paste0(dirname(gwr_fine@fit.name),"/",gwr_fine@name,"_smoothed",smooth.fact,"x_",resamp.method,".tif")
      gwr_pansharpen@ps.multiplier.name <- paste0(dirname(gwr_fine@fit.name),"/",gwr_fine@name,"_pansharpener",smooth.fact,"x_",resamp.method,".tif")
  }



  # gwr_pansharpen@fit.name <- paste0(dirname(gwr_coarse@fit.name),"/",gwr_pansharpen@name,"_preds.tif")
  # if(ps.strength!=1) gwr_pansharpen@fit.name <- sub(".tif",paste0("_",round(ps.strength*100),"strength.tif"),gwr_pansharpen@fit.name)
  # gwr_pansharpen@fine.smooth.name <- sub(".tif",paste0("_smooth",smooth.fact,"x.tif"),gwr_fine@fit.name)
  # gwr_pansharpen@ps.multiplier.name <- sub(".tif",paste0("_pansharpener",smooth.fact,"x.tif"),gwr_fine@fit.name)

  if(file.exists(gwr_pansharpen@fit.name) & overwrite) file.remove(gwr_pansharpen@fit.name)
  if(file.exists(gwr_pansharpen@fine.smooth.name) & overwrite.smooth) file.remove(gwr_pansharpen@fine.smooth.name)
  if(file.exists(gwr_pansharpen@ps.multiplier.name) & overwrite.smooth) file.remove(gwr_pansharpen@ps.multiplier.name)
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
        sys.command <- paste("gdal_translate -q -a_nodata -9999 -r ",resamp.method," -projwin",Ext[1],Ext[4],Ext[2],Ext[3],"-outsize",paste(dims[2:1],collapse=" "),"-a_srs",proj,"-co 'COMPRESS=LZW' -co 'BIGTIFF=YES'",tmpname,gwr_pansharpen@fine.smooth.name)
        system(sys.command) # 208s
        file.remove(tmpname)
      }
      t1 <- unclass(Sys.time())
      bn1 <- paste0("-A ",gwr_fine@fit.name," -B ",gwr_pansharpen@fine.smooth.name)
      bn2 <- paste0("-A ",gwr_coarse@fit.name," -B ",gwr_pansharpen@ps.multiplier.name)
      form <- if(ps.strength==1) 'A*B' else paste0("'A*((B-1)*",ps.strength,"+1)'")
      sys.cmd1 <- paste0("gdal_calc.py ",bn1," --outfile=",gwr_pansharpen@ps.multiplier.name," --calc='numpy.where(B==0, 1, A/B)' --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
      sys.cmd2 <- paste0("gdal_calc.py ",bn2," --outfile=",gwr_pansharpen@fit.name," --calc=",form," --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
      cat("calculating pansharpener\n")
      if(!file.exists(gwr_pansharpen@ps.multiplier.name)) system(sys.cmd1)
      cat("calculating final raster\n")
      system(sys.cmd2)

      # fn <- c(gwr_coarse@fit.name,gwr_fine@fit.name,gwr_pansharpen@fine.smooth.name)
      # lets <- LETTERS[-8][1:length(fn)]
      # bn <- paste(paste0("-",lets," ",fn),collapse=" ")
      # form <- "'numpy.where(C==0, 0, maximum(A*B/C,0))'"
      # sys.command <- paste0("gdal_calc.py ",bn," --outfile=",gwr_pansharpen@fit.name," --calc=",form," --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
      # cat("calculating final raster using",form,"\n")
      # system(sys.command) # 6min
      cat("done\n\n")
    } else {
      if(!file.exists(gwr_pansharpen@fine.smooth.name)){
          if(verbose) cat("making smoothed version of finescale\n")
          if(fine.prj!=coarse.prj) img.fine.smoothed <- terra::project(terra::aggregate(terra::rast(gwr_fine@fit.name),smooth.fact),terra::rast(gwr_coarse@fit.name),"bilinear",
                                                                       filename=gwr_pansharpen@fine.smooth.name,gdal=c("COMPRESS=LZW","BIGTIFF=YES","TFW=NO"))
          if(fine.prj==coarse.prj) img.fine.smoothed <- terra::resample(terra::aggregate(terra::rast(gwr_fine@fit.name),smooth.fact),terra::rast(gwr_coarse@fit.name),"bilinear",
                                                                        filename=gwr_pansharpen@fine.smooth.name,gdal=c("COMPRESS=LZW","BIGTIFF=YES","TFW=NO"))
          t1 <- unclass(Sys.time())
        }
      if(verbose) cat("performing main calc\n")
      terra::writeRaster(terra::rast(gwr_coarse@fit.name)*terra::rast(gwr_fine@fit.name)/terra::rast(gwr_pansharpen@fine.smooth.name),filename=gwr_pansharpen@fit.name,gdal=c("COMPRESS=LZW","BIGTIFF=YES","TFW=NO"))
      }
  } else t1 <- unclass(Sys.time())
  if(verbose) cat("calculating stats\n")
  t2 <- unclass(Sys.time())
  # gwr_pansharpen@rdata$fine.fit <- gwr_fine@rdata$rast.fit
  gwr_pansharpen@rdata$coarse.fit <- gwr_coarse@rdata$rast.fit
  gwr_pansharpen@rdata$ps.fit <- if(do.stats) textract(terra::rast(gwr_pansharpen@fit.name),gwr_pansharpen@rdata) else NA
  gwr_pansharpen@rdata$fine.smoothed.fit <- if(do.stats) textract(rast(gwr_pansharpen@fine.smooth.name),gwr_pansharpen@rdata) else NA
  gwr_pansharpen@rast.stats <- if(do.stats) {c(rsq=rsq(gwr_pansharpen@rdata[,rn],gwr_pansharpen@rdata$ps.fit),mae=mean(abs(gwr_pansharpen@rdata$ps.fit-gwr_pansharpen@rdata[,rn]),na.rm=T),
                                     mae_rast_vs_loo=mean(abs(gwr_pansharpen@rdata$loo.fit-gwr_pansharpen@rdata$ps.fit),na.rm=T))} else 0
  tt <- round(c(t.resamp=t1-t0,t.calc=t2-t1,t.stats=unclass(Sys.time())-t2,t.tot=unclass(Sys.time())-t0),1)
  if(verbose & do.stats) {cat("done\n");print(c(gwr_pansharpen@rast.stats,tt))}
  if(verbose & !do.stats) {cat("done\n");print(tt)}
  return(gwr_pansharpen)
}












