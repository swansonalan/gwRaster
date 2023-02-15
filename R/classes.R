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

check_gwr <- function(object){
  errors <- character()
  vn <- all.vars(as.formula(object@form))
  if(!all(vn %in% names(object@rdata))){
    msg <- paste("data.frame 'rdata' is missing variables found in ",object@form)
    errors <- c(errors,msg)
  }
  if (length(errors) == 0) TRUE else errors
}


#' An S4 class to for fitting gwr
#'
#' @slot form a standard R formula in character format e.g. z~x+y
#' @slot rdata a data.frame containing the variables in 'form', along with columns named 'lat' and 'lon'.
#' @slot name a character vector of length 1 used in any filenames created.
#' @slot out.dir a file path for outputs. Should end with a '/'
#' @slot fit a vector of loocv fits matching rows of rdata
#' @slot coefs an array of fitted loocv coefficient values.  First dimension corresponds to the points actually evaluated.
#' @slot loostats a vector of loocv fit statistics.
#' @slot bw fixed bandwidth in km.  Can be left blank if 'abw' is set.  Error will be triggered if both 'bw' and 'abw' are set.
#' @slot abw adaptive bandwidth size.  Sets bandwidth as radius required to capture n neighbors.  Can be left blank if 'bw' is specified.
#' @export
setClass(
  Class = "gwr",
  representation(
    form="character",
    rdata="data.frame",
    name="character",
    # out.dir="character",
    # fit="numeric",
    # coefs="array",
    loostats="numeric",
    bw="numeric",
    adaptive="logical",
    fit_idx="numeric",
    bw_used="numeric",
    buffer="numeric"
    ),
  prototype(name="dummy",bw=NaN,adaptive=FALSE),
  validity=check_gwr
)


#' An S4 class to for a raster of GWR coefficients
#'
#' @slot tmplt.name filename of a raster grid over which to fit GWR coefficient values
#' @slot coef.name filename of the coefficient raster
#' @slot band.name character vector of names for each of the bands of the coefficient raster
#' @slot hillshade.name filename of a hillshade raster with same extent as the coefficient raster
#' @export
setClass(
  Class = "gwr.coef.fit",
  representation(
      tmplt.name="character",
      coef.name="character",
      band.names="character",
      hillshade.name="character",
      out.dir="character"
    ),
  contains="gwr"
)

#' An S4 class to for a raster of GWR coefficients
#'
#' @slot tmplt.name filename of a raster grid over which to fit GWR coefficient values
#' @slot coef.name filename of the coefficient raster
#' @slot band.name character vector of names for each of the bands of the coefficient raster
#' @slot hillshade.name filename of a hillshade raster with same extent as the coefficient raster
#' @export
setClass(
  Class = "gwr.raster.fit",
  representation(
    img.fn="character",
    fit.name="character",
    fine.smooth.name="character",
    gaussian.smooth.cmd="character",
    ps.multiplier.name="character",
    temp.coef.names="character",
    rast.stats="numeric"
  ),
  contains="gwr.coef.fit"
)
