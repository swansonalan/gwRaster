
# functions ~~~

# gwr.obj object:  slots=form,name,rdata,out.dir

# gwr.loocv.fit <- fct(form,rdata,idx,name) # return(list(form,rdata,idx,loostats,fit,name))
# adds slots +loostats +loofit +coefs
# gwr.coef.map <- fct(form,rdata,rgrid,name)  # return(list(form,rdata,idx,loostats,fit,name))
# adds slots +tmplt.name +coef.name
# gwr.map <- fct(form,rdata,rgrid,img.names,name,sub.name) # return(list(form,rdata,idx,loostats,name)))
# adds slots for sub.name, which has slots for +fit.name +fitstats +hill.name +rast.fit
# sub.hill <- fct(ext,res,name,out.dir) return(hill.name)
# diag.plot <- fct(gwr.obj,hill)
# coef.plot <- fct(gwr.obj,hill,band)
# fit.plot <- fct(gwr.obj,hill)




library(geodist);library(foreach);library(doParallel);library(abind);library(raster)

stopCluster(cl)
cl <- makeCluster(49);registerDoParallel(cl)

server <- strsplit(Sys.info()[[4]],"\\.")[[1]][1]
rdir <- switch(server,bucephalus="/mnt/DataDrive1/data/wg/raster/",vulcan="/mnt/ScratchDrive/data/weather_generators/interpolation/data/")
xfer.dir <- "/mnt/orthanc1/water_balance/"
data.dir <- xfer.dir
out.dir <- paste0(rdir,"swe_gwr_fits/")

fn.8s <- paste0(rdir,c("swe_annual_mean_wyr0901_9221_annual_wus_8s.tif","et0_annual_sum_9221_wus_8s.tif","def_annual_sum_9221_wus_8s.tif"))
fn.30m <- paste0(rdir,c("dem_west_30m_crop.tif","gsr_west_30m_albers_reproj_rescale_crop.tif","swe_annual_mean_wyr0901_9221_annual_wus_30m.tif",
                        "et0_annual_sum_9221_wus_30m.tif","def_annual_sum_9221_wus_30m.tif",
                        "dem_x_gsr_wus_30m.tif","dem_squared_wus_30m.tif","swe_9221_squared_wus_30m.tif"))

names(fn.30m) <- c("dem30","gsr30","swe30","et30","def30","dem_x_gsr","dem_squared","swe_squared")

img <- raster(fn.30m["swe30"])^2
writeRaster(img,paste0(rdir,"swe_9221_squared_wus_30m.tif"),options=c("COMPRESS=LZW","BIGTIFF=YES","TFW=NO"))


mask.name<- paste0(rdir,"wus_dem_5km_mask_crop.tif")
mask.name10 <- paste0(rdir,"wus_dem_10km_mask_crop.tif")
mask.name20 <- paste0(rdir,"wus_dem_20km_mask_crop.tif")
dem.name.1k <- paste0(rdir,"dem_wus_1km_masked_crop.tif")

# all.fn <- c(mask.name,mask.name10,mask.name20,dem.name.1k)
# old.fn <- sub(".tif","old.tif",all.fn)
# file.rename(all.fn,old.fn)
# foreach(fn.in=old.fn,fn.out=all.fn) %do% {
#   z <- raster(fn.in)
#   projection(z) <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
#   writeRaster(z,fn.out,overwrite=T,options=c("COMPRESS=LZW","TIFW=NO"))
#   }
#

# define small raster subsets over which to apply gwr fit ~~~~~~~~~~~~~~~~~~~
subsets <- list(broot_sm=c(-114.43,-114.29,46.36,46.46),
                broot_lg=c(-114.54,-114.16,46.15,46.45),
                missoula=c(-114.42,-113.67,46.44,47.03),
                missions=c(-113.97,-113.81,47.25,47.37),
                glacier=c(-113.78,-113.68,48.71,48.78),
                shuksan=c(-121.64,-121.53,48.8,48.855),
                olympus=c(-123.77,-123.63,47.766,47.838),
                glacier_lg=c(-114.168,-113.52,48.55,49),
                n_cascades=c(-121.535,-121.07,48.702,49))
ss.dir <- "/mnt/DataDrive1/data/wg/raster/30m_subsets/"
out.names <- lapply(names(subsets),function(x) {y <- paste0(ss.dir,x,"/",sub(".tif",paste0("_",x,".tif"),basename(fn.30m)));
names(y) <- c("dem30","gsr30","swe30","et30","def30","dem_x_gsr","dem_squared","swe_squared");
y})
names(out.names) <- names(subsets)
dem.names <- sapply(out.names,function(x) x[1])
hill.names <- sub("dem","hill",dem.names)
img.FN <- foreach(fn=out.names,hill=hill.names) %do% c(fn,hill=hill)
save(img.FN,fn.30m,fn.8s,file=paste0(xfer.dir,"wus_gwr_raster_filenames_bucephalus.Rdata"))



hill <- hillshade(terra::rast(dem.name.1k))

print(load(paste0(data.dir,"cmip6_wb_gwr_data_100722.Rdata")))
rdata$swe_squared <- terra::extract(rast(fn.30m["swe_squared"]),project(vect(rdata,crs=geog),crs(rast(fn.30m[1]))))[,2]
gwr0 <- gwr_loocv("swe30~dem30+dem_squared+gsr30",rdata,bw=100,adaptive=T,out.dir=out.dir,n.samp=7000)
plot(gwr0)
coef0 <- gwr_coef_map(gwr0,tmplt.name=mask.name10,overwrite=F,verbose=T)
plot(coef0,hill=hill,crop=.95)
bs.fit <- gwr_fit_map(coef0,img.FN[[1]])
bl.fit <- gwr_fit_map(coef0,img.FN[[2]])
mission.fit <- gwr_fit_map(coef0,img.FN[[4]])
glacier.fit <- gwr_fit_map(coef0,img.FN[[5]])
shuksan.fit <- gwr_fit_map(coef0,img.FN[[6]])
oly.fit <- gwr_fit_map(coef0,img.FN[[7]])
gl.fit <- gwr_fit_map(coef0,img.FN[[8]])
nc.fit <- gwr_fit_map(coef0,img.FN[[9]])

plot(shuksan.fit,title="Mt Shuksan")
plot(oly.fit,title="Mt Olympus")

gwr1 <- gwr_loocv("swe2050~dem30+dem_squared+gsr30",rdata,bw=100,adaptive=T,out.dir=out.dir,n.samp=7000)
plot(gwr1)
coef1 <- gwr_coef_map(gwr1,tmplt.name=mask.name10,overwrite=F,verbose=T)
plot(coef1,hill=hill,crop=.95)
bs.fit <- gwr_fit_map(coef1,img.FN[[1]])
bl.fit <- gwr_fit_map(coef1,img.FN[[2]])
mission.fit <- gwr_fit_map(coef1,img.FN[[4]])
glacier.fit <- gwr_fit_map(coef1,img.FN[[5]])
shuksan.fit <- gwr_fit_map(coef1,img.FN[[6]])
oly.fit <- gwr_fit_map(coef1,img.FN[[7]])
gl.fit1 <- gwr_fit_map(coef1,img.FN[[8]])
nc.fit1 <- gwr_fit_map(coef1,img.FN[[9]])

gwr2 <- gwr_loocv("swe2050~swe30+dem30+dem_squared+gsr30",rdata,bw=100,adaptive=T,out.dir=out.dir,n.samp=7000)
plot(gwr2)
coef2 <- gwr_coef_map(gwr2,tmplt.name=mask.name10,overwrite=F,verbose=T)
plot(coef2,hill=hill,crop=.95)
bs.fit <- gwr_fit_map(coef2,img.FN[[1]])
bl.fit <- gwr_fit_map(coef2,img.FN[[2]])
mission.fit <- gwr_fit_map(coef2,img.FN[[4]])
glacier.fit <- gwr_fit_map(coef2,img.FN[[5]])
shuksan.fit <- gwr_fit_map(coef2,img.FN[[6]])
oly.fit <- gwr_fit_map(coef2,img.FN[[7]])
gl.fit2 <- gwr_fit_map(coef2,img.FN[[8]])
nc.fit2 <- gwr_fit_map(coef2,img.FN[[9]])

dpar(mar=c(1,1,3,1))
plot(gl.fit)

nc.fit <- gwr_fit_map(coef1,img.FN[[9]])

gwr3 <- gwr_loocv("swe2050~swe30+dem30+dem_squared+gsr30",rdata,bw=100,adaptive=F,out.dir=out.dir,n.samp=7000)
plot(gwr3)
coef3 <- gwr_coef_map(gwr3,tmplt.name=mask.name10,overwrite=F,verbose=T)
plot(coef3,hill=hill,crop=.95)
bs.fit <- gwr_fit_map(coef3,img.FN[[1]])
bl.fit <- gwr_fit_map(coef3,img.FN[[2]])
mission.fit <- gwr_fit_map(coef3,img.FN[[4]])
glacier.fit <- gwr_fit_map(coef3,img.FN[[5]])
shuksan.fit <- gwr_fit_map(coef3,img.FN[[6]])
oly.fit <- gwr_fit_map(coef3,img.FN[[7]])
gl.fit3 <- gwr_fit_map(coef3,img.FN[[8]])
nc.fit3 <- gwr_fit_map(coef3,img.FN[[9]])


gwr4 <- gwr_loocv("swe2050~swe30+swe_squared+dem30+dem_squared+gsr30",rdata,bw=100,adaptive=F,out.dir=out.dir,n.samp=7000)
plot(gwr4)
coef4 <- gwr_coef_map(gwr4,tmplt.name=mask.name10,overwrite=F,verbose=T)
plot(coef4,hill=hill,crop=.95)
bs.fit4 <- gwr_fit_map(coef4,img.FN[[1]])
bl.fit4 <- gwr_fit_map(coef4,img.FN[[2]])
gl.fit4 <- gwr_fit_map(coef4,img.FN[[8]])

gwr5b <- gwr_loocv("swe2050~swe30+swe_squared",rdata,bw=75,adaptive=F,out.dir=out.dir,n.samp=7000,buffer=10)
plot(gwr5b)
coef5b <- gwr_coef_map(gwr5b,tmplt.name=mask.name10,overwrite=T,verbose=T)
plot(coef5b,hill=hill,crop=.95)
bs.fit4 <- gwr_fit_map(coef4,img.FN[[1]])
bl.fit4 <- gwr_fit_map(coef4,img.FN[[2]])
gl.fit4 <- gwr_fit_map(coef4,img.FN[[8]])

gwr6 <- gwr_loocv("swe2050~swe30",rdata,bw=50,adaptive=F,out.dir=out.dir,n.samp=70000,buffer=20)
plot(gwr6)
coef6 <- gwr_coef_map(gwr6,tmplt.name=mask.name10,overwrite=T,verbose=T)
plot(coef6,hill=hill,crop=.95)

gl.fit6 <- gwr_fit_map(coef6,img.FN[[8]])

img.fine <- gl.fit1@fit.name
img.coarse <- gl.fit6@fit.name

# make raster subsets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prj <- projection(raster(fn.30m[1]))
subsets.prj <- lapply(subsets,function(x) xExtent(x,prj)[[3]])
out.dirs <- unique(dirname(unlist(out.names)))
for(dn in out.dirs) dir.create(dn)
in.names <- rep(fn.30m,length(subsets))
exts <- rep(subsets.prj,each=length(fn.30m))
table(ss <- !file.exists(unlist(out.names)))
foreach(fn.in=in.names[ss],fn.out=unlist(out.names)[ss],ext=exts[ss],.packages='raster') %dopar% writeRaster(crop(raster(fn.in),ext),fn.out,options=c("COMPRESS=LZW","TFW=NO"))
dem.names <- sapply(out.names,function(x) x[1])
hill.names <- sub("dem","hill",dem.names)
table(ss <- !file.exists(hill.names))
foreach(dn=dem.names[ss],hn=hill.names[ss],.packages='terra') %dopar% writeRaster(hillshade(terra::rast(dn)),hn,gdal=c("COMPRESS=LZW","TFW=NO"))
dem.names <- sapply(out.names,function(x) x[1])
hill.names <- sub("dem","hill",dem.names)
img.FN <- foreach(fn=out.names,hill=hill.names) %do% c(fn,hill=hill)






Pos <- setClass("Pos", slots = c(latitude = "numeric", longitude = "numeric", altitude = "numeric"))

Pos(latitude = c(1,2,3), longitude = c(4,5,6), altitude = c(7,8,9))


#' An S4 class to Logistic Regression.
#'
#' @export
#'
#' @slot coefficients Coefficients
#' @slot var Variance Covariance Matrix
#' @slot deviance Deviance
#' @slot predictors Predictors of the model
#' @slot iterations No of iterations for convergence

setClass(
  Class = "lreg5"
  , slots =  list(
    coefficients="numeric"
    , var="matrix"
    , deviance="numeric"
    , predictors="character"
    , iterations="numeric"
  )
)


lreg5 <-
  function(X, y, predictors=colnames(X), max.iter=10,
           tol=1E-6, constant=TRUE, ...) {
    if (!is.numeric(X) || !is.matrix(X))
      stop("X must be a numeric matrix")
    if (!is.numeric(y) || !all(y == 0 | y == 1))
      stop("y must contain only 0s and 1s")
    if (nrow(X) != length(y))
      stop("X and y contain different numbers of observations")
    if (constant) {
      X <- cbind(1, X)
      colnames(X)[1] <- "Constant"
    }
    b <- b.last <- rep(0, ncol(X))
    it <- 1
    while (it <= max.iter){
      p <- as.vector(1/(1 + exp(-X %*% b)))
      var.b <- solve(crossprod(X, p * (1 - p) * X))
      b <- b + var.b %*% crossprod(X, y - p)
      if (max(abs(b - b.last)/(abs(b.last) + 0.01*tol)) < tol) break
      b.last <- b
      it <- it + 1
    }
    if (it > max.iter) warning("maximum iterations exceeded")
    dev <- -2*sum(y*log(p) + (1 - y)*log(1 - p))
    result <- new("lreg5", coefficients=as.vector(b), var=var.b,
                  deviance=dev, predictors=predictors, iterations=it)
    result
  }

setMethod("show", signature(object="lreg5"),
          definition=function(object) {
            coef <- object@coefficients
            names(coef) <- object@predictors
            print(coef)
          }
)


setMethod("summary", signature(object="lreg5"),
          definition=function(object, ...) {
            b <- object@coefficients
            se <- sqrt(diag(object@var))
            z <- b/se
            table <- cbind(b, se, z, 2*(1-pnorm(abs(z))))
            colnames(table) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
            rownames(table) <- object@predictors
            printCoefmat(table)
            cat("\nDeviance =", object@deviance,"\n")
          }
)












