##################################################################
##################################################################
##
##  M-BITES (Mosquito Bout-based and Individual-based Transmission Ecology Simulation)
##  Version 0.9
##  April 6, 2017
##
##  This version was designed and written by David L. Smith (aka.
##  Dave), Sean Wu, and Hector Sanchez.
##  Please send bug reports, comments, and suggestions to
##  <smitdave@gmail.com> or <slwu89@berkeley.edu>.
##
##  Robert C. Reiner, Jr. (aka Bobby) <bcreiner@uw.edu>, Hector
##  Sanchez Castellanos <sanchez.hmsc@gmail.com> Sean Wu
##  <slwu89@berkeley.edu>, and Amit Verma <amit.verma13@gmail.com>
##  helped with development, debugging and documentation of
##  version 1.0.
##
##  M-BITES (formerly DHM) was conceived of by David Smith, and it was inspired
##  by discussions with many people, including Bobby, Hector,
##  Sean, Amit, Arnaud Le Menach, Nick Ruktanonchai, Samson
##  Kiware, Gerry Killeen, Tom Scott, Ellis McKenzie, Steven W.
##  Lindsay, Willem Takken, Philip Eckhoff, Nicole Achee, Chris
##  Barker, Nakul Chitnis, Justin Cohen, Su Yun Kang, Audrey
##  Lenhart, John Marshall, Phil McCall, Catherine Moyes, Doug
##  Norris, Alex Perkins, Chris Stone, Edward Wenger, and Anne
##  Wilson.
##
##################################################################
##################################################################

#################################################################
#
#   MASH/MBITES
#   MBITES Male Life Cycle Routines
#   MBITES-Mating.R defines female mating behaviors
#   R version
#   Sean Wu
#   January 27, 2017
#
#################################################################


#################################################################
#  MatingQ
#################################################################

# allocMatingQ: allocate a mating queue with nMales slots
allocMatingQ <- function(nMales){
  #. allocMatingQ: allocate a mating queue
  #nMales: number of empty slots in queue
  list(
    id  = rep(NA,nMales)
  )
}

# addMale2Q: add male ID: id to matingQ at site indexed by: ix
addMale2Q <- function(ix, id){

  emptyIx = which(is.na(LANDSCAPE$swarmSites[[ix]]$matingQ$id))
  if(length(emptyIx)==0){
    extendMatingQ(ix = ix)
    emptyIx = which(is.na(LANDSCAPE$swarmSites[[ix]]$matingQ$id))
  }
  qIx = emptyIx[1]
  addMale(qIx, ix, id)

}

# addMale: add id to queue position qIx at site ix
addMale <- function(qIx, ix, id){
  LANDSCAPE$swarmSites[[ix]]$matingQ$id[qIx] <<- id
}

# extendMatingQ: extend the matingQ for site ix
extendMatingQ <- function(ix, offset = NULL){

  if(is.null(offset)){
    offset = ceiling(length(LANDSCAPE$swarmSites[[ix]]$matingQ$id)*1.5)
  }

  LANDSCAPE$swarmSites[[ix]]$matingQ$id <<- c(LANDSCAPE$swarmSites[[ix]]$matingQ$id,rep(NA,offset))

}

# clearMatingQ: clear all matingQs
clearMatingQ <- function(){
  for(ix in 1:LANDSCAPE$nM){
    qSize = length(LANDSCAPE$swarmSites[[ix]]$matingQ$id)
    LANDSCAPE$swarmSites[[ix]]$matingQ$id <<- rep(NA,qSize)
  }
}


#################################################################
#  Male Survival
#################################################################

# SenesceMale:
SenesceMale <- function(M,P){
  #. SenesceMale
  age = M$tNow - M$bDay
  with(P,{
    return((2+sns.b)/(1+sns.b) - exp(sns.a*age)/(sns.b + exp(sns.a*age)))
  })
}

# maleFlightStress
maleFlightStress <- function(M,P){
  #. maleFlightStress
  if(isActive(M)){
    p = getSFp(M,P) # baseline survival
    if(P$TATTER){
      M$damage = M$damage + rTatterSize(P)
      p = p * pTatter(M,P)
    }
    if(P$SENESCE){
      p = p * SenesceMale(M,P)
    }
    if(!rbinom(1,1,p)){
      M$stateNew = "D"
    }
  }
  return(M)
}


################################################################################
#  Generic Bout:
#
#  boutF is any bout function and ... are additional named
#  parameters to be passed to boutF.
#  The generic bout runs necessary updates of timing, state,
#  survival, energetics, and queue checks prior to calling boutF
#  and checks mosquito is alive/active before calling the bout
#
#  This corresponds to the following Gillespie-style algorithm:
#   1. tNow is set to tNext from previous bout
#   2. rMove: movement between point classes (if needed)
#   3. boutFun: run bout function
#   4. run energetics and check if alive
#   5. run landingSpot and check if alive
#   6. run surviveResting/surviveFlight and check if alive
#   7. update tNext
#   8. update state to stateNew which is determined in boutGeneric
#
################################################################################

# boutGeneric updates tNext and stateNew
boutGenericMale <- function(M,P,boutFun,...){

  # update time and state
  M$tNow = M$tNext # update time
  M$state = M$stateNew # update current state
  M = timingExponential(M,P) # update tNext

  # movement
  move = rMove(ix = M$ix,pSet = M$inPointSet,bState = M$state)
  M$ix = move$ixNew
  M$inPointSet = move$pSetNew

  # bout
  M = boutFun(M,P,...)

  # energetics
  M = energetics(M,P) # MBITES-Energetics.R

  # landing spot
  M = landingSpot(M,P) # MBITES-Bouts.R

  # survival
  M = surviveResting(M,P) # MBITES-Survival.R
  M = maleFlightStress(M,P) # MBITES-Swarming.R

  # check queueing
  # queueEstivation()
  # queueMating() only queue if mated = FALSE

  # log history
  if(P$HISTORY){M = historyTrack(M)}

  return(M)
}


#################################################################
#  Swarming Bout :: M
#################################################################

# enterSwarm:
enterSwarm <- function(M,P){
  addMale2Q(M$ix, M$id)
}

# boutMM:
boutMM <- function(M,P){
  if(isAlive(M) && rbinom(1,1,P$maleM.s)){
    enterSwarm(M,P)
    M$stateNew = "S"
  }
  # M=fSwarmSpraying(M)
  return(M)
}


#############################################################################
#   MBITES
#   Male Mosquito
#   Mosquito Bout-based & Individual-based Transmission Ecology Simulation
#############################################################################

# MBITESmalei:
# tMax: gloabl time tick; ie all the mosy start born at time=0 and ticks start 1,2,3,...
# run bouts while tNext (time of next bout/action) < global tMax time tick
# run bouts while mosquito is still alive
MBITESmalei <- function(M,P){
  while(M$tNext < tMax && M$stateNew != "D"){
    M = switch(M$stateNew,
      M = boutGenericMale(M,P,boutMM),
      S = boutGenericMale(M,P,boutS),
      R = boutGenericMale(M,P,boutR)
    )
  }
  return(M)
}

# # MBITESmale:
MBITESmale <- function(P){

  clearMatingQ() # clear matingQs prior to running daily M-BITES

  # iterate over non-null, living mosquitoes
  ixM = which(sapply(MPopM$mosy,function(x){!is.null(x$id) & isAlive(x)}))

  for(ix in ixM){
    MPopM$mosy[[ix]] <<- MBITESmalei(MPopM$mosy[[ix]],P)
  }

}
