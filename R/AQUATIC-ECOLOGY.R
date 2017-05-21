#################################################################
#
#   MASH/MBITES
#   Aquatic Ecology
#   Main Structure Definitions
#   R version
#   Sean Wu
#   February 2, 2017
#
#################################################################


#############################################
#  Setup the ImagoQ and associated functions
#############################################

#' Empty Adult Slot for ImagoQ
#'
#' This function is a low-level utility to generate a null placeholder value for adults in ImagoQ. It is called by allocImagoQ and extendImagoQ.
#'
#' @param N number of emerging adults from this element of queue
#' @param tm time of emergence
#' @param ix index of aquatic habitat
#' @param dam id of mother
#' @param sire id of father
#' @return a list of named elements
#' @examples
#' newAdults(N=0,tm=0,ix=0,dam=0,sire=0)
newAdults = function(N=0,tm=0,ix=0,dam=0,sire=0){
  list(
  N=N,
  tm=tm,
  ix=ix,
  dam=dam,
  sire=sire
)}

#' Zero out a space in ImagoQ
#'
#' This function is a low-level utility to zero out a position in ImagoQ after emergence. It is called by addCohort.
#'
#' @param ixQ index of position in ImagoQ to zero out
#' @param ix index of aquatic habitat of emergence
#' @return modifies global LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]
#' @examples
#' zeroAdult(ixQ=1,ix=1)
zeroAdult <- function(ixQ, ix){
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$N     <<- 0
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$tm    <<- 0
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$ix    <<- 0
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$dam   <<- 0
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$sire  <<- 0
}


#' Insert Adults into ImagoQ
#'
#' This function is a low-level utility to add a cohort of adults to ImagoQ. It is called by addAdults2Q.
#'
#' @param ixQ index of position in ImagoQ to zero out
#' @param ix index of aquatic habitat of emergence
#' @param N number of emerging adults from this element of queue
#' @param tm time of emergence
#' @param ix index of aquatic habitat
#' @param dam id of mother
#' @param sire id of father
#' @return modifies global LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]
#' @examples
#' addAdult(ixQ,N,tm,ix,dam,sire)
addAdult <- function(ixQ, ix, N, tm, dam, sire){
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$N     <<- N
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$tm    <<- tm
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$ix    <<- ix
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$dam   <<- dam
  LANDSCAPE$aquaSites[[ix]]$ImagoQ[[ixQ]]$sire  <<- sire
}


#allocImagoQ: allocate ImagoQ
allocImagoQ <- function(N){
  ImagoQ = replicate(n=N,expr=newAdults(),simplify=FALSE)
  return(ImagoQ)
}

#extendImagoQ: extend ImagoQ
extendImagoQ <- function(ix,offset = NULL){

  if(is.null(offset)){
    offset = ceiling(length(LANDSCAPE$aquaSites[[ix]]$ImagoQ)*1.5)
  }

  LANDSCAPE$aquaSites[[ix]]$ImagoQ <<- c(LANDSCAPE$aquaSites[[ix]]$ImagoQ,replicate(n = offset,expr = newAdults(),simplify = FALSE))
}


#addAdults2QEmerge: add emerging adults to ImagoQ for "emerge" module
#lambda: vector of length equal to nTypes of emerging adults
#t: time of emergence
#ix: index of aquatic habitat
#dam: maternal id
#sire: paternal id
addAdults2Q <- function(lambda, tm, ix, dam, sire){

  emptyIx = which(sapply(LANDSCAPE$aquaSites[[ix]]$ImagoQ,function(x){x$N==0}))
  if(length(emptyIx)==0){
    extendImagoQ(ix = ix)
    emptyIx = which(sapply(LANDSCAPE$aquaSites[[ix]]$ImagoQ,function(x){x$N==0}))
  }
  ixQ = emptyIx[1]
  addAdult(ixQ = ixQ, ix = ix,N = lambda,tm = tm, dam = dam,sire = sire)

}


# addMosquito: add a single mosquito to the cohort and adjust size / null indices
# t: time of emergence
# ix: index of home aquatic habitat
# state: state of mosy upon emergence
# EIP:
# inPointSet:
# female:
addMosquito <- function(tm, ix, state = "M", EIP = 10, inPointSet = "l", female = TRUE){
  if(female){
    # extend the mosquito cohort object
    if(length(MPopF$nullIx) < 2){
      enlargeMosyPop(female = TRUE)
    }
    MPopIx = MPopF$nullIx[1]
    # check for overwrite
    if(!is.null(MPopF$mosy[[MPopIx]]$id)){
      stop("overwriting non-NULL mosquito or overrunning NULL slots")
    }
    # create new mosquito
    # MPopF$mosy[[MPopIx]]$id          <<- MPopIx
    MPopF$mosy[[MPopIx]]$id          <<- paste0(tMax,"_",MPopIx)
    MPopF$mosy[[MPopIx]]$bDay        <<- tm
    MPopF$mosy[[MPopIx]]$tNow        <<- tm
    MPopF$mosy[[MPopIx]]$tNext       <<- tm
    MPopF$mosy[[MPopIx]]$ix          <<- ix
    MPopF$mosy[[MPopIx]]$state       <<- state
    MPopF$mosy[[MPopIx]]$stateNew    <<- state
    MPopF$mosy[[MPopIx]]$EIP         <<- EIP
    MPopF$mosy[[MPopIx]]$inPointSet  <<- inPointSet
    MPopF$mosy[[MPopIx]]$female      <<- female
    MPopF$mosy[[MPopIx]]$history$ixH <<- ix
    # adjust null index vector
    MPopF$nullIx                     <<- MPopF$nullIx[-1]

  } else {

    # extend the mosquito cohort object
    if(length(MPopM$nullIx) < 2){
      enlargeMosyPop(female = FALSE)
    }
    MPopIx = MPopM$nullIx[1]
    # check for overwrite
    if(!is.null(MPopM$mosy[[MPopIx]]$id)){
      browser("overwriting non-NULL mosquito or overrunning NULL slots")
    }
    # create new mosquito
    # MPopM$mosy[[MPopIx]]$id          <<- MPopIx
    MPopM$mosy[[MPopIx]]$id          <<- paste0(tMax,"_",MPopIx)
    MPopM$mosy[[MPopIx]]$bDay        <<- tm
    MPopM$mosy[[MPopIx]]$tNow        <<- tm
    MPopM$mosy[[MPopIx]]$tNext       <<- tm
    MPopM$mosy[[MPopIx]]$ix          <<- ix
    MPopM$mosy[[MPopIx]]$state       <<- state
    MPopM$mosy[[MPopIx]]$stateNew    <<- state
    MPopM$mosy[[MPopIx]]$inPointSet  <<- inPointSet
    MPopM$mosy[[MPopIx]]$female      <<- female
    MPopM$mosy[[MPopIx]]$history$ixH <<- ix
    # adjust null index vector
    MPopM$nullIx                     <<- MPopM$nullIx[-1]

  }
}

# addCohort: add a cohort of mosquitoes based on ImagoQ
# tNow: current time tick
addCohort <- function(tNow, ...){

  for(ii in 1:LANDSCAPE$nA){ # ii: iterate over aquatic habitats

    for(jj in 1:length(LANDSCAPE$aquaSites[[ii]]$ImagoQ)){ # jj: iterate over ImagoQ

      lambda = LANDSCAPE$aquaSites[[ii]]$ImagoQ[[jj]]$N # adult emergence from site ii; ImagoQ jj
      if(lambda == 0 | LANDSCAPE$aquaSites[[ii]]$ImagoQ[[jj]]$t > tNow){ # skip if no emergence or time incorrect
        next()
      } else { # add lambda mosquitoes to cohort
        # browser("addCohort being called when there are mosy to add")
        for(kk in 1:lambda){
          addMosquito(tm = tNow,ix = ii, female = TRUE, ...)
          addMosquito(tm = tNow,ix = ii, female = FALSE, ...)
        }

      }
      zeroAdult(ixQ = jj, ix = ii) # clear out ImagoQ
    } # jj: end iteration over ImagoQ
  } # ii: end iteration over aquatic habitats
}


###########################################
#  Setup the EggQ and associated functions
###########################################

# addBatch2Q: called from makeBatches() in MBITES-Energetics
addBatch2Q <- function(eggs, ix, tm, dam, sire){
  #. addBatch2Q

  emptyIx = which(sapply(LANDSCAPE$aquaSites[[ix]]$EggQ,function(x){x$N==0}))
  if(length(emptyIx)==0){
    extendEggQ(ix = ix)
    emptyIx = which(sapply(LANDSCAPE$aquaSites[[ix]]$EggQ,function(x){x$N==0}))
  }
  ixQ = emptyIx[1]
  addBatch(ixQ = ixQ, ix = ix, N = eggs, tm = tm, dam, sire)
}

eggBatch <- function(){
  list(
        N       = 0, # Eggs
        tm      = 0, # Time of oviposition
        ix      = 0, # ix of aqua site
        dam     = 0, # id of dam
        sire    = 0  # id of sire
  )
}

# addBatch: add an egg batch to the landscape
addBatch <- function(ixQ, ix, N, tm, dam, sire){
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$N     <<- N
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$tm    <<- tm
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$ix    <<- ix
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$dam   <<- dam
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$sire  <<- sire
}

# zeroBatch: zero out an egg batch
zeroBatch <- function(ixQ, ix){
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$N     <<- 0
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$tm    <<- 0
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$ix    <<- 0
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$dam   <<- 0
  LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$sire  <<- 0
}

# allocEggQ:
allocEggQ <- function(N){
  batches = replicate(n = N,expr = eggBatch(),simplify = FALSE)
  return(batches)
}

#extendEggQ: extend EggQ
extendEggQ <- function(ix,offset = NULL){

  if(is.null(offset)){
    offset = ceiling(length(LANDSCAPE$aquaSites[[ix]]$EggQ)*1.5)
  }

  LANDSCAPE$aquaSites[[ix]]$EggQ <<- c(LANDSCAPE$aquaSites[[ix]]$EggQ,replicate(n = offset,expr = eggBatch(),simplify = FALSE))
}


#################################################################
#
# Numerical tracking of EggQ:
# observe values of state varibles across each site in landscape
#
#################################################################

#' Initialize Text Connection to Track EggQ
#'
#' Generate a text connection to write .csv EggQ values to. This records daily egg laying in each aquatic habitat.
#' The text connection generated must be named \code{EggQCon}; this function also defines a global flag \code{EggQ_TRACK} that is
#' used by \code{oneDay_EL4P()} or \code{oneDay_emerge()} to determine whether or not to write to output.
#'
#' @param directory directory that will be created; files will be put in directory/OUTPUT/..
#' @param fileName name of the file to write to; directory/OUTPUT/fileName.csv
#' @return a text connection
#' @examples
#' trackEggQ_init(directory, fileName)
trackEggQ_init <- function(directory, fileName){
  if(!dir.exists(paste0(directory,"OUTPUT"))){
    dir.create(paste0(directory))
    dir.create(paste0(directory,"OUTPUT"))
  }
  if(file.exists(paste0(directory,"OUTPUT/",fileName))){
    stop("trackEggQ_init cannot init a file that already exists!")
  }
  EggQ_TRACK <<- TRUE
  conEggQ = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  writeLines(text = paste0(c("time",paste0("ix",1:LANDSCAPE$nA)),collapse = ","),con = conEggQ,sep = "\n") # write header
  return(conEggQ)
}

#' Write EggQ Values to .csv
#'
#' Generate daily output of EggQ tracking to .csv file through text connection created with \code{trackEggQ_init()}.
#'
#' @param con text connection to write to
#' @return none
#' @examples
#' trackEggQ(con)
trackEggQ <- function(con){
  eggs = sapply(LANDSCAPE$aquaSites,function(x){
    sum(sapply(x$EggQ,function(xx){xx$N}))
  })
  writeLines(text = paste0(c(tMax,eggs),collapse = ","),con = .GlobalEnv$EggQCon,sep = "\n")
}

#' Import Logged EggQ Data from .csv
#'
#' Generate daily output of EggQ tracking to .csv file through text connection created with \code{trackEggQ_init()}.
#'
#' @param directory directory that OUTPUT/.. is in
#' @param fileName file name in OUTPUT/..
#' @return data frame
#' @examples
#' importimportEggQEL4P(directory, fileName)
importEggQ <- function(directory, fileName){
  read.csv(file = paste0(directory,"OUTPUT/",fileName),header = TRUE)
}
