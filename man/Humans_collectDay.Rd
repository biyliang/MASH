% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PATHOGEN-PFSI-Auxiliary.R
\name{Humans_collectDay}
\alias{Humans_collectDay}
\title{Find State Status for All Humans on Specific Day as a data.table}
\usage{
Humans_collectDay(history, day)
}
\arguments{
\item{day}{what day}

\item{human}{a human history object}
}
\value{
data.table of states
}
\description{
Given all human histories (which can be extracted via \code{\link{convertHumanEvent_PfSI}} or \code{\link{importHumanEvent_PfSI}}),
retrieve their states on that day as a \code{\link{data.table}}.
}
\examples{
Humans_getInfoDay(history = humanHistories,day = 1)
}
