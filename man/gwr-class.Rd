% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{gwr-class}
\alias{gwr-class}
\title{An S4 class to for fitting gwr}
\description{
An S4 class to for fitting gwr
}
\section{Slots}{

\describe{
\item{\code{form}}{a standard R formula in character format e.g. z~x+y}

\item{\code{rdata}}{a data.frame containing the variables in 'form', along with columns named 'lat' and 'lon'.}

\item{\code{name}}{a character vector of length 1 used in any filenames created.}

\item{\code{out.dir}}{a file path for outputs. Should end with a '/'}

\item{\code{fit}}{a vector of loocv fits matching rows of rdata}

\item{\code{coefs}}{an array of fitted loocv coefficient values.  First dimension corresponds to the points actually evaluated.}

\item{\code{loostats}}{a vector of loocv fit statistics.}

\item{\code{bw}}{fixed bandwidth in km.  Can be left blank if 'abw' is set.  Error will be triggered if both 'bw' and 'abw' are set.}

\item{\code{abw}}{adaptive bandwidth size.  Sets bandwidth as radius required to capture n neighbors.  Can be left blank if 'bw' is specified.}
}}

