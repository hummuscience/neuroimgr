% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculators.R
\name{getBaseParameters}
\alias{getBaseParameters}
\title{Calculate base parameters.}
\usage{
getBaseParameters(cexp, slot = "normalized", binslot = "l0spikes", ...)
}
\arguments{
\item{cexp}{CalciumExperiment object to calculate base parameters from. Should be normalized and binarized before running.}

\item{slot}{}

\item{binslot}{}

\item{...}{}
}
\value{
A CalciumExperiment object with calculated base parameters added to the colData DataFrame
}
\description{
Calculates base parameters of a CalciumExperiment object such as mean, max, min responses as well as number of predicted spikes and location of first spike.
}
\author{
Muad Abd El Hay
}
