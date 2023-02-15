
######################################################################################
# test pan-sharpening using aggregated 8s data ~~~~~~~~~~~~~~~
######################################################################################

setwd("/home/aswanson/code/gwRaster")
# devtools::load_all()
library(geodist);library(foreach);library(doParallel);library(abind);library(terra);library(gwRaster)
server <- strsplit(Sys.info()[[4]],"\\.")[[1]][1]
rdir <- switch(server,bucephalus="/mnt/DataDrive1/data/wg/raster/",vulcan="/mnt/ScratchDrive/data/weather_generators/interpolation/data/")
fig.dir <- "/mnt/orthanc1/water_balance/wb_gwr/"
xfer.dir <- "/mnt/orthanc1/water_balance/"
data.dir <- xfer.dir
out.dir <- paste0(rdir,"swe_gwr_fits/")
tab.dir <- paste0(rdir,"gwr_tab/")
ss.dir <- "/mnt/DataDrive1/data/wg/raster/30m_subsets/"
mask.name10 <- paste0(rdir,"wus_dem_10km_mask_crop.tif")
dem.name.1k <- paste0(rdir,"dem_wus_1km_masked_crop.tif")
hill <- hillshade(rast(dem.name.1k))
# writeRaster(hill,file=paste0(rdir,"hill_wus_1km_masked_crop.tif"),gdal=c("COMPRESS=LZW","TFW=NO"))
print(load(paste0(xfer.dir,"wus_gwr_raster_filenames_bucephalus_102822.Rdata")));t(abind(subsets,along=2))
print(load(paste0(data.dir,"cmip6_wb_gwr_data_102822.Rdata")))

list.files(path=rdir,pattern="240m.tif$")

#
fn <- fn.30m[c("dem30","gsr30","et09221")]
tmplt <- "/mnt/DataDrive1/data/wg/raster//dem_west_240m_crop.tif"
# tmplt8 <- aggregate(rast(tmplt),8)
# writeRaster(tmplt8,paste0(rdir,"test/dem_west_2km_crop.tif"))
tmplt8 <- paste0(rdir,"test/dem_west_2km_crop.tif")

# dir.create(paste0(rdir,"test"))
fn240 <- paste0(rdir,"test/",sub("30m","240m",basename(fn)))
fn2km <- paste0(rdir,"test/",sub("30m","2km",basename(fn)))
fn2km_240 <- paste0(rdir,"test/",sub("30m","2km-to-240m",basename(fn)))
Out <- foreach(fn.in=fn,fn.out=fn240,fn.out2=fn2km,fn.out3=fn2km_240,.packages='raster',.export="gwarp",.errorhandling='pass') %dopar% {
  img <- gwarp(fn.in,tmplt,fn.out,execute=T)  # make 240m raster
  out <- textract(img,rdata,method="bilinear")
  img2 <-gwarp(fn.out,tmplt8,fn.out2)        # make 2km raster
  out2 <- textract(img2,rdata,method="bilinear")
  img3 <- gwarp(fn.out2,tmplt,fn.out3)    # resample back to 240m
  return(cbind(x240=out,x1920=out2))
}
out <- abind(Out,along=2)
colnames(out) <- paste0(rep(c("dem","gsr","et0"),each=2),c("240m","2km"))
out[out< -100] <- NA
rdata <- cbind(rdata,out)
names()
save(rdata,fn,fn240,fn2km,fn2km_240,file=paste0(rdir,"test/testdata.Rdata"))
# coarse model: et0~et02km
# fine model: et0~dem240+gsr240

print(load(paste0(rdir,"test/testdata.Rdata")))
rdata <- rdata[complete.cases(rdata),]
set.seed(1234)
idx <- sort(sample(which(complete.cases(rdata)),5000))

cl <- makeCluster(80)
registerDoParallel(cl)

# model 240m et0 using same resample to 2km ~~~~~~~
resp <- "et0240m"
gf250 <- gwr_loocv(paste0(resp,"~dem240m+gsr240m"),rdata,bw=250,idx=idx,adaptive=F,buffer=5,fast=T) # 0.866 ****************
plot(gf250)
cf250 <- gwr_coef_map(gf250,tmplt.name=mask.name10,overwrite=T,verbose=T,fast=T)
plot(cf250,hill=hill,crop=.95)
summary(values(rast(cf250@coef.name)))

gc50a <- gwr_loocv("et0240m~et02km",rdata,bw=50,idx=1:nrow(rdata),adaptive=F,buffer=5) # 0.962 ***********
plot(cc50a,hill=hill,crop=.95)
gc_coarse <- gc50a
gc_fine <- gf250
coef_coarse <- cc50a
coef_fine <- cf250
save(gc_coarse,gc_fine,coef_coarse,coef_fine,resp,file=paste0(tab.dir,"def2050_gwr_102722.Rdata"))

et02050 <- list(gwr_coarse=gc50a,gwr_fine=gf250,coef_coarse=cc50a,coef_fine=cf250,response=resp)
save(et02050,gc50a,file=paste0(rdir,"test/et0_9221_gwr_103122.Rdata"))

print(load(paste0(rdir,"test/et0_9221_gwr_103122.Rdata")))
print(load(paste0(rdir,"test/testdata.Rdata")))
et02050$coef_coarse@form
et02050$coef_fine@form
names(fn240) <- c("dem240m","gsr240m","et0240m")
names(fn2km_240) <-c("dem2km","gsr2km","et02km")
gwr_coarse <- gwr_fit_map(et02050$coef_coarse,fn2km_240,overwrite=T,overwrite.resamps=T)
gwr_fine <- gwr_fit_map(et02050$coef_fine,fn240,overwrite=T)
gwr_ps <- make_ps(gwr_fine,smooth.fact=8,resamp.method="cubicspline")
gwr_ps2 <- make_ps(gwr_fine,smooth.fact=8,resamp.method="bilinear")

fit.name <- "/mnt/DataDrive1/data/wg/raster/test/gwr_fits/et0240m_vs_dem240m_gsr240m_bw250_preds.tif"
for(sf in c(4,6,8)){
    smooth.name <- paste0("/mnt/DataDrive1/data/wg/raster/test/gwr_fits/et0240m_vs_dem240m_gsr240m_bw250_smoothed",sf,"_gaussian.tif")
    ps.name  <- paste0("/mnt/DataDrive1/data/wg/raster/test/gwr_fits/et0240m_vs_dem240m_gsr240m_bw250_pansharpener",sf,"_gaussian.tif")
    sys.cmd1 <- paste0("otbcli_Smoothing -in ",fit.name," -out ",smooth.name," float -type gaussian -type.gaussian.stdev ",sf)
    system(sys.cmd1)
    bn <- paste("-A",fit.name,"-B",smooth.name)
    sys.cmd2 <- paste0("gdal_calc.py ",bn," --outfile=",ps.name," --calc='numpy.where(B==0, 1, A/B)' --NoDataValue=-9999 --co COMPRESS=LZW --co BIGTIFF=YES --overwrite --quiet")
    system(sys.cmd2)
}

rd <- gc50a@rdata
ps.names <- c(paste0("/mnt/DataDrive1/data/wg/raster/test/gwr_fits/et0240m_vs_dem240m_gsr240m_bw250_pansharpener",c(2,4,6,8),"_gaussian.tif"),
              gwr_ps@ps.multiplier.name,gwr_ps2@ps.multiplier.name)
out <- textract(rast(ps.names),rd,method="bilinear")
colnames(out) <- paste0("ps_",fdates(colnames(out),n=7),sub("pansharpener","",fdates(colnames(out),n=6)))
colnames(out)[1:4] <- paste0(colnames(out)[1:4],"x")
rd <- cbind(rd[,1:9],out)
x <- seq(-1,1.5,by=.1)
rv <- foreach(i=10:15,.combine='cbind') %do% {
  foreach(k=x,.combine='c') %do% rsq(rd$et0240m,rd$loo.fit*((rd[,i]-1)*k+1))
}
colnames(rv)<-colnames(rd)[10:15]
rv <- cbind(k=x,rv)
apply(rv[,-1],2,function(x) c(rv[which.max(x),1],max(x)))

x <- seq(0,.75,by=.01)
rv2 <- foreach(i=10:15,.combine='cbind') %do% {
  foreach(k=x,.combine='c') %do% rsq(rd$et0240m,rd$loo.fit*((rd[,i]-1)*k+1))
}
colnames(rv2)<-colnames(rd)[10:15]
rv2 <- cbind(k=x,rv2)
apply(rv2[,-1],2,function(x) c(rv2[which.max(x),1],max(x)))
save(rd,rv,rv2,file=paste0(rdir,"test/resultsdata_110222.Rdata"))


# pan-sharpening ~~~
fn <- c(dem240m=fn240[1],gsr240m=fn240[2],et02km=fn2km_240[3])
print(load(paste0(tab.dir,"et02050_gwr_102822.Rdata")))
gwr_coarse <- gwr_fit_map(et02050$coef_coarse,fn,overwrite=T,overwrite.resamps=T)
gwr_fine <- gwr_fit_map(et02050$coef_fine,fn,overwrite=T)
gwr_ps <- pan_sharpen(gwr_coarse,gwr_fine,overwrite=T)
save(et02050,gwr_coarse,gwr_fine,gwr_ps,file=paste0(rdir,"test/et0_9221_gwr_103122.Rdata"))

# model 2050 240m et0 using 92-21 resample to 2km ~~~~~~~
resp <- "et2050"
gf250 <- gwr_loocv(paste0(resp,"~dem240m+gsr240m"),rdata,bw=250,idx=idx,adaptive=F,buffer=5,fast=T) # 0.866 ****************
plot(gf250)
cf250 <- gwr_coef_map(gf250,tmplt.name=mask.name10,overwrite=T,verbose=T,fast=T)
plot(cf250,hill=hill,crop=.95)
summary(values(rast(cf250@coef.name)))

gc50a <- gwr_loocv("et2050~et02km",rdata,bw=50,idx=1:nrow(rdata),adaptive=F,buffer=5) # 0.962 ***********
cc50a <- gwr_coef_map(gc50a,tmplt.name=mask.name10,overwrite=F,fast=T,verbose=T)
plot(cc50a,hill=hill,crop=.95)
et02050 <- list(gwr_coarse=gc50a,gwr_fine=gf250,coef_coarse=cc50a,coef_fine=cf250,response=resp)
save(et02050,file=paste0(rdir,"test/et0_2050_gwr_103122.Rdata"))

# pan-sharpening ~~~
fn <- c(dem240m=fn240[1],gsr240m=fn240[2],et02km=fn2km_240[3])
print(load(paste0(rdir,"test/et0_2050_gwr_103122.Rdata")))
gwr_coarse <- gwr_fit_map(et02050$coef_coarse,fn,overwrite=T,overwrite.resamps=T)
gwr_fine <- gwr_fit_map(et02050$coef_fine,fn,overwrite=T)
gwr_ps <- pan_sharpen(gwr_coarse,gwr_fine,overwrite=T)
gwr_ps16 <- pan_sharpen(gwr_coarse,gwr_fine,smooth.fact=16,overwrite=T)
save(et02050,gwr_coarse,gwr_fine,gwr_ps,gc50a,file=paste0(rdir,"test/et0_2050_gwr_103122.Rdata"))

rsq(gwr_ps@rdata$et2050,gwr_coarse@rdata$rast.fit)

rsq(gwr_ps@rdata$et2050,gwr_ps@rdata$fine.fit)
rsq(gwr_ps@rdata$et2050,gwr_ps@rdata$coarse.fit)
rsq(gwr_ps@rdata$et2050,gwr_ps@rdata$ps.fit)
rsq(gwr_ps@rdata$et2050,gwr_coarse@rdata$loo.fit)
rsq(gwr_ps@rdata$et2050,gwr_coarse@rdata$loo.fit*gwr_ps@rdata$fine.fit/gwr_ps@rdata$fine.smoothed.fit)
rdata <- gwr_ps

rd <- cbind(gc50a@rdata,gwr_ps@rdata[,c("fine.fit","fine.smoothed.fit","coarse.fit","ps.fit")])
rd$ps <- rd$fine.fit/rd$fine.smoothed.fit
rd$ps16 <- gwr_ps16@rdata$fine.fit/gwr_ps16@rdata$fine.smoothed.fit
rd$ps16[gwr_ps16@rdata$fine.smoothed.fit>3000] <- NA

dpar(mfrow=c(3,1))
hist(rd$loo.err,breaks=40)
hist(rd$ps,breaks=40)
plot(loo.err~ps,data=rd[sample(nrow(rd),5000),])
plot(loo.err~ps,data=rd[abs(rd$loo.err)>75,],ylim=c(-250,250))
points(loo.err~ps,data=rd[sample(nrow(rd),5000),])
brks <- c(seq(-225,250,by=25),755)
z<-hist(rd$loo.err,breaks=brks,plot=F)
ss <- foreach(i=1:(length(brks)-2),.combine='c') %do% {
  idx <- which(rd$loo.err>brks[i] & rd$loo.err<brks[i+1])
  if(length(idx)<1000) idx else sample(idx,1000)
}
brks <- c(seq(500,2500,by=25),755)
ss <- foreach(i=1:(length(brks)-2),.combine='c') %do% {
  idx <- which(rd$loo.err>brks[i] & rd$loo.err<brks[i+1])
  if(length(idx)<1000) idx else sample(idx,1000)
}

rd$loo.fit.ps <-  rd$loo.fit*rd$ps
rd$loo.fit.ps2 <- rd$loo.fit*((rd$ps-1)/2+1)
rd$loo.err.ps <- rd$loo.fit.ps-rd$et2050
rd$loo.err.ps2 <- rd$loo.fit.ps2-rd$et2050
rsq(rd$et2050,rd$loo.fit)
rsq(rd$et2050,rd$loo.fit.ps)
rsq(rd$et2050,rd$loo.fit.ps2)
rd$neg.loo.err <- -1*rd$loo.err
m1 <- lm(loo.err~ps,data=rd)
m2 <- lm(neg.loo.err~ps-1,data=rd)
rd$loo.fit.ps3 <- rd$loo.fit+fitted(m1)
rd$ratio <- rd$loo.fit/rd$ps
rd$log.ratio <- log(rd$loo.fit/rd$ps)
m3 <- lm(ratio~ps,data=rd)
m4 <- lm(ratio~ps-1,data=rd)
plot(ratio~ps,data=rd[ss,])
abline(coef(m3))
plot(log.ratio~ps,data=rd[ss,])

rd$ratio <- rd$loo.fit/rd$ps
m3 <- lm(ratio~ps,data=rd)
rd$loo.fit.ps <-  rd$loo.fit*rd$ps
rd$z <- rd$loo.fit*rd$ps
plot(et2050~loo.fit.ps,data=rd[ss,])
m4 <- lm(et2050~loo.fit.ps,data=rd)
m5 <- lm(et2050~loo.fit.ps-1,data=rd)
abline(coef(m4),col=2,lty=2,lwd=2)
abline(c(0,coef(m5)))
rd$log.et2050 <- log(rd$et2050)
rd$log.loo.fit <- log(rd$loo.fit)
rd$log.ps <- log(rd$ps)
m6 <- lm(log.et2050~log.loo.fit,data=rd)
rd$log.coarse.err <- fitted(m6)-rd$log.et2050
dpar(mfrow=c(2,1))
plot(log.et2050~log.loo.fit,data=rd[ss,]);abline(coef(m6),lty=2,col=2,lwd=2)
plot(log.coarse.err~log.ps,data=rd[ss,])
summary(m7 <- lm(log.coarse.err~log.ps,data=rd))
summary(m8 <- lm(log.coarse.err~log.ps-1,data=rd))
summary(m9 <- lm(log.et2050~log.loo.fit+log.ps,data=rd,na.action=na.exclude))
rsq(rd$log.et2050,fitted(m9))
rsq(rd$et2050,exp(fitted(m9)))
rd$loo.fit.ps4 <- exp(fitted(m9))
(k4 <- coef(lm(rd$et2050~rd$loo.fit.ps4-1,data=rd)))
rd$loo.fit.ps5 <- k4*exp(fitted(m9))
rd$loo.err.ps5 <- rd$loo.fit.ps5-rd$et2050
x <- foreach(k=seq(.1,1.5,by=.01),.combine='c') %do% rsq(rd$et2050,rd$loo.fit*(k*(rd$ps-1)+1))
z <- data.frame(k=seq(.1,1.5,by=.01),r=x)
plot(r~k,data=z,type="l",ylim=c(.975,1))
abline(h=rsq(rd$et2050,rd$loo.fit),col=2,lty=2)

k <- z$k[which.max(z$r)] # .42
rd$loo.fit.ps3 <- rd$loo.fit*(k*(rd$ps-1)+1)
rd$loo.err.ps3 <- rd$loo.fit.ps3-rd$et2050
rsq(rd$et2050,rd$loo.fit.ps5)

m10 <- lm(et2050~loo.fit+ps,data=rd,na.action=na.exclude)
rd$loo.fit.ps6 <- fitted(m10)
rd$loo.err.ps6 <- rd$loo.fit.ps6-rd$et2050
rsq(rd$et2050,rd$loo.fit.ps6) # .9926


png(paste0(fig.dir,"pansharpeing_err_effect_110122.png"),9,9,units="in",res=180)
dpar(mfrow=c(2,2),cex.main=.8)
plot(loo.err~ps,data=rd[ss,],xlab="pansharpening coef",ylab="LOO error (mm)",main="coarse model, R^2=0.991")
abline(v=1,h=0,lty=2,col=2,lwd=2)
plot(loo.err.ps~ps,data=rd[ss,],xlab="pansharpening coef",ylab="LOO error (mm)",main="pan-sharpened coarse model, R^2=.990")
abline(v=1,h=0,lty=2,col=2,lwd=2)
plot(loo.err.ps3~ps,data=rd[ss,],xlab="pansharpening coef",ylab="LOO error (mm)",main="42% pan-sharpened coarse model, R^2=.993")
abline(v=1,h=0,lty=2,col=2,lwd=2)
plot(loo.err.ps5~ps,data=rd[ss,],xlab="pansharpening coef",ylab="LOO error (mm)",main="log-est pan-sharpened coarse model, R^2=.993")
abline(v=1,h=0,lty=2,col=2,lwd=2)
dev.off()


#######################################################################################################
# glacier example ~~~~~~~~~~~~~~~~
#######################################################################################################

rdata <- rdata[,-grep("_mo",names(rdata))]
rdata$et2050 <- rdata$et2050/365.24
summary(rdata)




gwr_ps16 <- pan_sharpen(gwr_coarse,gwr_fine,overwrite=T,smooth.fact=16)
gwr_ps32 <- pan_sharpen(gwr_coarse,gwr_fine,overwrite=T,smooth.fact=32)
gwr_ps64 <- pan_sharpen(gwr_coarse,gwr_fine,overwrite=T,smooth.fact=64)

hill <- rast(img.FN[["glacier"]]["hill"])
dpar(mfrow=c(2,2),mar=c(1,1,3,1))
rplot(rast(gwr_fine@fit.name),hill=hill,main="fine gwr")
rplot(rast(gwr_ps@fine.smooth.name),hill=hill,main="smoothed fine gwr")
rplot(rast(gwr_ps@ps.multiplier.name),hill=hill,main="pansharpener")
rplot(rast(gwr_ps@fit.name),hill=hill,main="pansharped")

dpar(mfrow=c(2,2),mar=c(1,1,3,1))
rplot(rast(gwr_fine@fit.name),main="fine gwr")
rplot(rast(gwr_ps@fine.smooth.name),main="smoothed fine gwr")
rplot(rast(gwr_ps@ps.multiplier.name),main="pansharpener")
rplot(rast(gwr_ps@fit.name),main="pansharpened fit")

dpar(mfrow=c(2,2),mar=c(1,1,3,1))
rplot(rast(gwr_ps@ps.multiplier.name),main="8x pansharpener")
rplot(rast(gwr_ps16@ps.multiplier.name),main="16x pansharpener")
rplot(rast(gwr_ps32@ps.multiplier.name),main="64x pansharpener")
rplot(rast(gwr_ps64@ps.multiplier.name),main="64x pansharpener")


