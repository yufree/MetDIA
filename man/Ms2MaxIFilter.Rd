% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Ms2MaxIFilter.R
\name{Ms2MaxIFilter}
\alias{Ms2MaxIFilter}
\title{Remove fragment ions with correlation below 0.6
 replace intensity of fragment ion whose correlation is below 0.6 with 0}
\usage{
Ms2MaxIFilter(ms2maxi, corresult)
}
\arguments{
\item{ms2maxi}{a list contains fragment ion intensity at the summit}

\item{corresult}{a list contains correaltion between each fragment ion and precursor ion}
}
\description{
Remove fragment ions with correlation below 0.6
 replace intensity of fragment ion whose correlation is below 0.6 with 0
}

