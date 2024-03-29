% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwr_fit_map.R
\name{gwr_fit_map}
\alias{gwr_fit_map}
\title{Makes GWR coefficient estimates on a grid}
\usage{
gwr_fit_map(
  x,
  img.fn,
  verbose = T,
  overwrite = F,
  overwrite.resamps = F,
  do.stats = F
)
}
\arguments{
\item{img.fn}{filename of a raster to use as grid, or a raster object}

\item{verbose}{logical; whether or not to report progress}

\item{overwrite}{logical; whether or not to overwrite existing raster fit}

\item{object}{a gwr.coef.fit object, typically created using gwr_coef_map()}
}
\value{
a gwr object with the original inputs, original data, fit, statistics.
}
\description{
This function fits a gwr model to a data.frame which includes geographic information
}
