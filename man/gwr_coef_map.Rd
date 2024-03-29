% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwr_coef_map.R
\name{gwr_coef_map}
\alias{gwr_coef_map}
\title{Makes GWR coefficient estimates on a grid}
\usage{
gwr_coef_map(
  x,
  tmplt.name,
  name = "",
  out.dir = NULL,
  buffer = 3,
  verbose = T,
  overwrite = F,
  fast = T
)
}
\arguments{
\item{tmplt.name}{filename of a raster to use as grid, or a raster object}

\item{verbose}{logical; whether or not to report progress}

\item{object}{a gwr object, typically created using gwr_loocv()}
}
\value{
a gwr object with the original inputs, original data, fit, statistics.
}
\description{
This function fits a gwr model to a data.frame which includes geographic information
}
