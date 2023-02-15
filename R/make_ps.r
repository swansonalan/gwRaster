#' Makes a pansharpener raster, which is the ratio of a finescale fit divided by a smoothed version of the same fit.
#' This is accomplished by aggregating a
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
make_ps <- function(gwr_fine,smooth.fact=8,verbose=T,overwrite=F,overwrite.smooth=F,resamp.method="bilinear",do.stats=F){
  t0 <- unclass(Sys.time())
  # if(!inherits(gwr_coarse,"gwr.coef.fit")) stop("gwr_coarse needs to be a gwr.coef.fit object")
  if(!inherits(gwr_fine,"gwr.coef.fit")) stop("gwr_fine needs to be a gwr.coef.fit object")
  # if(!test_cluster(F)) stop("need cluster")
  # if(foreach::getDoParWorkers()==1)  stop("need cluster")
  # rn <- strsplit(gwr_coarse@form,"~")[[1]][1]
  # vn <- all.vars(as.formula(gwr_coarse@form))[-1]
  '%dopar%' <- foreach::'%dopar%'
  '%do%' <- foreach::'%do%'

  gwr_pansharpen <- gwr_fine
  gwr_pansharpen@name <- paste0(gwr_fine@name,"_pansharpener",smooth.fact,"x_",resamp.method)


  gwr_pansharpen@fit.name <- paste0(dirname(gwr_fine@fit.name),"/",gwr_fine@name,"_pansharpener",smooth.fact,"x_",resamp.method,".tif")
  gwr_pansharpen@fine.smooth.name <- paste0(dirname(gwr_fine@fit.name),"/",gwr_fine@name,"_smoothed",smooth.fact,"x_",resamp.method,".tif")
  gwr_pansharpen@ps.multiplier.name <-  gwr_pansharpen@fit.name

  gaussian.smooth.name <- paste0(dirname(gwr_fine@fit.name),"/",gwr_fine@name,"_smoothed",smooth.fact,"x_gaussian.tif")
  gwr_pansharpen@gaussian.smooth.cmd <- paste0("otbcli_Smoothing -in ",gwr_fine@fit.name," -out ",gaussian.smooth.name," float -type gaussian -type.gaussian.maxerror 0.5 -type.gaussian.maxwidth ",smooth.fact*2," -type.gaussian.stdev ",smooth.fact)

  if(file.exists(gwr_pansharpen@fit.name) & overwrite) file.remove(gwr_pansharpen@fit.name)
  if(file.exists(gwr_pansharpen@fine.smooth.name) & overwrite.smooth) file.remove(gwr_pansharpen@fine.smooth.name)
  if(file.exists(gwr_pansharpen@ps.multiplier.name) & overwrite.smooth) file.remove(gwr_pansharpen@ps.multiplier.name)
  if(!file.exists(gwr_pansharpen@fit.name)){
      if(!file.exists(gwr_pansharpen@fine.smooth.name)){
        if(resamp.method!="gaussian"){
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

        if(resamp.method=="gaussian"){
          if(inherits(try(system(gwr_pansharpen@gaussian.smooth.cmd,intern=T)),"try-error")){
            cat("Error: orfeo toolbox doesn't seem to be installed properly.\n")
            cat("try pasting the following command into an orfeo environment then running again:\n")
            print(gwr_pansharpen@gaussian.smooth.cmd)
            return(gwr_pansharpen)
          }
        }
    }
    t1 <- unclass(Sys.time())
    bn1 <- paste0("-A ",gwr_fine@fit.name," -B ",gwr_pansharpen@fine.smooth.name)
    sys.cmd1 <- paste0("gdal_calc.py ",bn1," --outfile=",gwr_pansharpen@ps.multiplier.name," --calc='numpy.where(B==0, 1, A/B)' --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
    cat("calculating pansharpener\n")
    system(sys.cmd1)
  } else t1 <- unclass(Sys.time())
  if(verbose) cat("calculating stats\n")
  t2 <- unclass(Sys.time())
  gwr_pansharpen@rdata$fine.fit <- gwr_fine@rdata$rast.fit
  if(do.stats){
      gwr_pansharpen@rdata$fine.smoothed.fit <- if(do.stats) textract(rast(gwr_pansharpen@fine.smooth.name),gwr_pansharpen@rdata) else 0
      gwr_pansharpen@rdata$ps <- if(do.stats) textract(rast(gwr_pansharpen@fit.name),gwr_pansharpen@rdata) else 0
  }
  gwr_pansharpen@rast.stats <- 0
  tt <- round(c(t.resamp=t1-t0,t.calc=t2-t1,t.tot=unclass(Sys.time())-t0),1)
  # if(verbose) cat(paste0("done, time used: ",tt[4]," seconds (",paste(tt[1:3],collapse="s;"),")\n\n"))
  if(verbose) {cat("done\n");print(tt)}
  return(gwr_pansharpen)
}
