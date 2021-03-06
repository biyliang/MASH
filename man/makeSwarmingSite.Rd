% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MICRO-LANDSCAPE.R
\name{makeSwarmingSite}
\alias{makeSwarmingSite}
\title{Generate a Sugar Feeding Site Object}
\usage{
makeSwarmingSite(siteIx, siteXY, searchWt, landHaz, nMales = 5)
}
\arguments{
\item{siteIx}{index of this site}

\item{siteXY}{length 2 numeric vector of xy-coordinates}

\item{searchWt}{search weight for MvOb kernel}

\item{landHaz}{resting or landing hazard}

\item{nMales}{number of males to pre-allocate for MatingQ}
}
\value{
a named list
}
\description{
This wil create a single sugar feeding site object for MICRO-LANDSCAPE
}
\examples{
makeSwarmingSite(siteIx, siteXY, searchWt, landHaz, nMales = 5)
}
