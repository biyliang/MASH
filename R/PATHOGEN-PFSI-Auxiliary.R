#################################################################
#
#   MASH-MAP relevant functions
#   Sean Wu & Su Yun Kang
#   May 19, 2017
#
#################################################################


#################################################################
# Simulate Biting Processes with desired statistical properties
#################################################################

#' Gamma Distributed Biting Rates
#'
#' Produce gamma-distributed rates for each person, with a chosen mean value of Poisson process parameterized by gamma rates, and selected shape (controls variance). This function calls \code{\link{rgamma}} to sample random numbers.
#'
#' @param N number of humans
#' @param mean mean value
#' @param shape controls variance in gamma distributed rates; larger values produce distributions with lower variance.
#' @return distribution of biting rates
gammaRates <- function(N, mean = 5,shape = 10){
  scale = mean / shape
  1/rgamma(n = N,shape = shape,scale = scale)
}

#' Simulate Biting Process from Population Mean Waiting Time Between Bites
#'
#' Simulate a biting process to result in a negative-binomially distributed cumulative number of bites per person over \code{days}.
#' The biting process for each person is a Poisson process with individual gamma-distributed rates.
#' Because the negative binomial distribution arises as a compound Poisson-Gamma distribution, it will also arise from independent Poisson processes each with gamma-distributed rates,
#' because the expected value of a Poisson process is also Poisson.
#' This calls \code{\link{gammaRates}}.
#'
#' @param nH number of humans
#' @param meanWaitingTime population mean waiting time between bites
#' @param days number of days
#' @param shape shape parameter of gamma-distributed variation in individual biting rate (higher values lead to more normally distributed cumulative biting counts)
#' @param plot visualize process?
#' @return a list of length \code{nH}; each element is the arrivial times of bites on that individual
simBites_waitingTime <- function(nH, meanWaitingTime, days = 365, shape = 5, plot = TRUE){

  rRates = gammaRates(N = nH,mean = meanWaitingTime,shape = shape)

  # each human's biting times are a Poisson process with their own Gamma-distributed rate parameter
  bitingTimes = vector(mode="list",length = nH)

  for(i in 1:nH){
    t = 0
    times = NULL
    print(paste0("i: ",i))
    while(t < days){
      t = t + rexp(n = 1,rate = rRates[i])
      times = c(times,t)
    }
    bitingTimes[[i]] = times
  }

  if(plot){
    par(mfrow=c(1,2))

    hist(x = 1/rRates,probability = TRUE,main = paste0("Mean Waiting Times Between Bites \n(1 / rate): ",signif(mean(1/rRates),digits = 3)),xlab= "Time (Days)")

    nBites = sapply(bitingTimes,length)
    hist(nBites,breaks = 20,probability = TRUE,main = paste0("Number of Bites in ",days," Days","\n Mean number of bites: ",mean(nBites)),xlab="Count")
    fitNbinom = MASS::fitdistr(x = nBites,densfun = "negative binomial")
    nBinomDens = dnbinom(x = 0:max(nBites),size = fitNbinom$estimate[["size"]],mu = fitNbinom$estimate[["mu"]])
    lines(x = 0:max(nBites),y = nBinomDens,col="red",lwd=1.75)

    par(mfrow=c(1,1))
  }

  return(bitingTimes)
}


#' Simulate Biting Process from Population Mean Number of Bites
#'
#' Simulate a biting process to result in a negative-binomially distributed cumulative number of bites per person over \code{days}.
#' The biting process for each person is a Poisson process with individual gamma-distributed rates.
#' Because the negative binomial distribution arises as a compound Poisson-Gamma distribution, it will also arise from independent Poisson processes each with gamma-distributed rates,
#' because the expected value of a Poisson process is also Poisson.
#' This calls \code{\link{gammaRates}}.
#'
#' @param nH number of humans
#' @param meanNumberBites population mean number of bites (over time period equal to \code{days})
#' @param days number of days
#' @param shape shape parameter of gamma-distributed variation in individual biting rate (higher values lead to more normally distributed cumulative biting counts)
#' @param plot visualize process?
#' @return a list of length \code{nH}; each element is the arrivial times of bites on that individual
simBites_MeanBites <- function(nH, meanNumberBites, days = 365, shape = 5, plot = TRUE){

  meanRate = meanNumberBites/days # mean biting rate
  rRates = gammaRates(N = nH,mean = 1/meanRate,shape = shape) # generate gamma-distributed rates

  # each human's biting times are a Poisson process with their own Gamma-distributed rate parameter
  bitingTimes = vector(mode="list",length = nH)

  for(i in 1:nH){
    t = 0
    times = NULL
    print(paste0("i: ",i))
    while(t < days){
      t = t + rexp(n = 1,rate = rRates[i])
      times = c(times,t)
    }
    bitingTimes[[i]] = times
  }

  if(plot){
    par(mfrow=c(1,2))

    hist(x = 1/rRates,probability = TRUE,main = paste0("Mean Waiting Times Between Bites \n(1 / rate): ",signif(mean(1/rRates),digits = 3)),xlab= "Time (Days)")

    nBites = sapply(bitingTimes,length)
    hist(nBites,breaks = 20,probability = TRUE,main = paste0("Number of Bites in ",days," Days","\n Mean number of bites: ",mean(nBites)),xlab="Count")
    fitNbinom = MASS::fitdistr(x = nBites,densfun = "negative binomial")
    nBinomDens = dnbinom(x = 0:max(nBites),size = fitNbinom$estimate[["size"]],mu = fitNbinom$estimate[["mu"]])
    lines(x = 0:max(nBites),y = nBinomDens,col="red",lwd=1.75)

    par(mfrow=c(1,1))
  }

  return(bitingTimes)
}


#################################################################
# Convienience Functions for analysis
#################################################################

#' Find State Status for Specific Human on Specific Day
#'
#' Given a single human history (which can be extracted via \code{\link{convertHumanEvent_PfSI}} or \code{\link{importHumanEvent_PfSI}}),
#' retrieve their state on that day.
#'
#' @param human a single human history
#' @param day what day
#' @return character vector of state
#' @examples
#' Humans_getInfoDay(human = humanHistories[[1]],day = 1)
Humans_getInfoDay <- function(human, day) {
  index <- max(which(human$eventT == max(human$eventT[human$eventT <= day])))
  human$events[index]
}

#' Find State Status for All Humans on Specific Day as a data.table
#'
#' Given all human histories (which can be extracted via \code{\link{convertHumanEvent_PfSI}} or \code{\link{importHumanEvent_PfSI}}),
#' retrieve their states on that day as a \code{\link{data.table}}.
#'
#' @param human a human history object
#' @param day what day
#' @return data.table of states
#' @examples
#' Humans_getInfoDay(history = humanHistories,day = 1)
Humans_collectDay <- function(history, day) {
  data.table::data.table(day = day, event = sapply(history, Humans_getInfoDay, day = day))
}

#' Return Human State Histories as a data.table
#'
#' Given all human histories (which can be extracted via \code{\link{convertHumanEvent_PfSI}} or \code{\link{importHumanEvent_PfSI}}),
#' retrieve the cohort state histories as a \code{\link{data.table}}.
#'
#' @param history a human history object
#' @param day_integers integer vector of days
#' @return data.table of states
#' @examples
#' Humans_getInfoDay(history = humanHistories,day = -1L:tMax)
Humans_collectOutput <- function(history, day_integers) {
  results <- data.table::rbindlist(lapply(day_integers, Humans_collectDay, history = history))
  data.table::dcast(results, day~event, fun.aggregate = length, value.var = "event")
}
