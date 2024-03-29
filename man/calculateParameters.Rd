% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculators.R
\name{calculateParameters}
\alias{calculateParameters}
\title{Calculate stimulus parameters from a split CalciumExperiment (i.e. MultiAssayExperiment) object.}
\usage{
calculateParameters(mcexp, ...)
}
\arguments{
\item{mcexp}{MultiAssayExperiment object containing multiple CalciumExperiment objects (one per stimulus) to calculate from.}

\item{...}{}
}
\value{
A MultiAssayExperiment object with calculated base parameters for each CalciumExperiment object.
}
\description{
Function that takes a MultiAssayExperiment object containings CalciumExperiment objects and calculates parameters for each assay (stimulus).
}
\author{
Muad Abd El Hay
}
