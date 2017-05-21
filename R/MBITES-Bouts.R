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
#   Bout routines
#   R version
#   Sean Wu
#   March 6, 2017
#
#################################################################


#################################################################
#
#  HOUSE ENTERING & RESTING BEHAVIOR:
#  At the end of the search bout, attempt bout, or after Egg
#  laying a mosquito has entered the area around a feeding
#  station and either rested or attempted to rest:
#    l) Leave the area
#    r) Reattempt Without Resting;
#    v) Rest on vegetation
#    w) Rest on the Outside wall of a structure
#    i) Rest on the Inside wall of a structure
#
#################################################################

#getWTS: return landing spot weights based on behavior state
getWTS <- function(M,P){
  with(P,{
    switch(M$state,
      F = Fwts,
      B = Fwts,
      R = Rwts,
      L = Lwts,
      O = Owts,
      M = Mwts,
      S = Swts
    )
  })
}

#lspotIx: return integer index of lspot
lspotIx <- function(lspot){
  switch(lspot,
    i = 1L,
    w = 2L,
    v = 3L,
    r = 4L,
    l = 5L
  )
}

# newSpot: return choose landing spot for focal mosquito
newSpot <- function(M,P){
  if(M$inPointSet == "f"){
    probs = P$InAndOut[lspotIx(M$lspot),] * getWTS(M,P)
    return(sample(x = c("i","w","v","r","l"),size = 1,prob = probs))
  } else {
    return("v")
  }
}

#enterHouse: mosquito attempts to enter a house and appropriate events are called
enterHouse <- function(M,P){
  if(runif(1) < LANDSCAPE$feedSites[[M$ix]]$enterHouseP){ #mosquito found a gap
    # M=fEaveTube(M)
    # M=fImproveHome(M)
  } else { #mosquito did not find a gap
    M$lspot = newSpot(M,P)
    M = surviveFlight(M,P)
    if(M$lspot == "i"){
      M = Recall(M,P)
    }
  }

  return(M)
}

# landingSpot:
landingSpot <- function(M,P){
  if(isActive(M)){
    oldSpot = M$lspot
    M$lspot = newSpot(M,P) # choose new lspot
    if(oldSpot != "i" & M$lspot == "i"){
      M = enterHouse(M,P) # enterHouse
    }
    # if(M$lspot == "l" & M$state == "B"){
    #   M$stateNew = "F" # repeat search bout (boutF)
    # }
    # if(M$lspot == "l" & M$state == "O"){
    #   M$stateNew = "L" # repeat search bout (boutL)
    # }
  }
  return(M)
}

############################################################
# Checks to determine eligibility for action
############################################################

isAlive <- function(M){
  if(M$stateNew == "D" || M$state == "D"){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

isActive <- function(M){
  if(M$state == "E"){
    return(FALSE)
  } else {
    return(isAlive(M))
  }
}


############################################################
# History & Bionomics
############################################################

# track history
historyTrack <- function(M){
  M$history$stateH     = c(M$history$stateH,M$state)     # state trajectory
  M$history$timeH      = c(M$history$timeH,M$tNow)         # transition times
  M$history$ixH        = c(M$history$ixH,M$ix)              # sites visited
  M$history$pSetH      = c(M$history$pSetH,M$inPointSet)    # point sets visited
  if(!isAlive(M)){ # return final state transition & calculate bionomics upon death
    M$history$stateH     = c(M$history$stateH,M$stateNew)     # state trajectory
    M$history$timeH      = c(M$history$timeH,M$tNext)         # transition times
    if(M$female){M$bionomics = bionomics(M$history)}
  }
  return(M)
}

# calculate bionomics upon death
bionomics <- function(history){
  with(history,{
    if(tail(stateH,n=1) != "D"){
      stop("mosquito not dead yet")
    }
    bionomics = list()
    if(!is.null(batchH)){
      bionomics$mBatch = mean(batchH) # mean batch size
      bionomics$tBatch = sum(batchH) # total batch size
    } else {
      bionomics$mBatch = 0L
      bionomics$tBatch = 0L
    }
    bionomics$feedAllH =  feedAllH # total number of bloodmeals
    bionomics$feedHumanH = feedHumanH # number of human bloodmeals
    bionomics$bmInt = diff(feedAllT) # all bloodmeal intervals
    bionomics$bmIntH = diff(feedHumanT) # human bloodmeal intervals
    bionomics$lifespan = tail(timeH,n=1) - head(timeH,n=1) # lifespan
    return(bionomics)
  })
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
boutGeneric <- function(M,P,boutFun,...){

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

  # # Intervetions routines
  # M=fIRS(M)
  # M=fAerialSpray(M)

  # survival
  M = surviveResting(M,P) # MBITES-Survival.R
  M = surviveFlight(M,P) # MBITES-Survival.R

  # check queueing
  M = queueEstivation(M,P) # MBITES-Estivate.R

  # log history
  if(P$HISTORY){M = historyTrack(M)}

  return(M)
}


#################################################################
#  Blood Feeding Search Bout :: F
#################################################################

boutF <- function(M,P){
  if(M$lspot != "l" && rbinom(n=1,size=1,prob=P$F.s)){
    M$stateNew = "B"
  }

  return(M)
}


#################################################################
#  Blood Feeding Attempt Bout :: B
#################################################################

boutB <- function(M,P){
  if(rbinom(1,1,P$B.s)){
    M = chooseHost(M) # MBITES-ChooseHost.R
  } else {
    M$hostID = 0
  }

  if(M$hostID > 0){
    M = humanEncounter(M,P) # MBITES-HostEncounter.R
  } else if(M$hostID == -1){
    M = zooEncounter(M,P) # MBITES-HostEncounter.R
  } else if(M$hostID == 0){
    M = nullEncounter(M,P) # MBITES-HostEncounter.R
  } else {
    stop("hostID of mosy id: ",M$id," not a recognized host ID")
  }

  return(M)
}


#################################################################
#  Post-Prandial Resting Bout :: R
#################################################################

boutR <- function(M,P){
  #. boutR: Mosquito resting bout
  if(!rbinom(1,1,P$R.p)){
    M$stateNew = "D"
  } else {
    if(M$female){ # female behavior
      if(!M$mature){ # immature female
        M$stateNew = sample(x = c("F","S"),size = 1)
      } else { # mature female
        M = reFeed(M,P) # MBITES-Energetics.R
      }
    } else { # male behavior
      M$stateNew = "M"
    }
  }

  return(M)
}

#################################################################
#  Egg Laying Search Bout :: L
#################################################################

boutL <- function(M,P){
  if(rbinom(1,1,P$L.s)){
    M$stateNew = "O"
  }

  return(M)
}


#################################################################
#  Egg Laying Attempt Bout :: O
#################################################################

layEggs <- function(M,P){
  if(rbinom(1,1,P$O.s)){
    makeBatches(M) # MBITES-Energetics.R
    M$batch = 0
    M$stateNew = "F"
  }

  return(M)
}

boutO <- function(M,P){
  # M=fOvitrap(M)
  if(isAlive(M)){
    M = layEggs(M,P)
  }
  return(M)
}


#################################################################
#  Sugar Feeding Attempt Bout :: O
#################################################################

boutS <- function(M,P){
  with(P,{
    if(M$female){

      # female behavior
      M$stateNew = sample(x = c("F","D"), size = 1, prob = c(S.s,1-S.s))
      if(isAlive(M)){
        M$energy = 1

        if(!M$mature){
          M$energyPreG = M$energyPreG - preGsugar
          if(M$energyPreG <= 0){
            M$mature = TRUE
          }
        }
      }

    } else {

      # male behavior
      M$stateNew = sample(x = c("R","D"), size = 1, prob = c(S.s,1-S.s))
      if(isAlive(M)){
        M$energy = 1
      }

    }

    # M=fATSB(M)
    return(M)
  })
}


#################################################################
#  Female Mating Bout :: M
#################################################################

boutMF <- function(M,P){
  M = chooseMate(M,P) # MBITES-Mating.R
  # M = fSwarmSpraying(M)
  return(M)
}


#################################################################
#  Estivation Bout :: E
#################################################################

boutE <- function(M,P){
  M$stateNew = "F"
  return(M)
}


#############################################################################
#   MBITES
#   Female Mosquito
#   Mosquito Bout-based & Individual-based Transmission Ecology Simulation
#############################################################################

# MBITESi: run M-BITES model for a single mosquito
# tMax: gloabl time tick; ie all the mosy start born at time=0 and ticks start 1,2,3,...
# run bouts while tNext (time of next bout/action) < global tMax time tick
# run bouts while mosquito is still alive
MBITESi <- function(M,P){

  while(M$tNext < tMax && M$stateNew != "D"){
    M = switch(M$stateNew,
      F = boutGeneric(M,P,boutF),
      B = boutGeneric(M,P,boutB),
      R = boutGeneric(M,P,boutR),
      L = boutGeneric(M,P,boutL),
      O = boutGeneric(M,P,boutO),
      M = boutGeneric(M,P,boutMF),
      S = boutGeneric(M,P,boutS),
      E = boutGeneric(M,P,boutE)
    )
  }
  return(M)
}

# MBITES: run M-BITES model for all mosquitoes
MBITES <- function(P){
  # iterate over non-null, living mosquitoes
  ixM = which(sapply(MPopF$mosy,function(x){!is.null(x$id) & isAlive(x)}))

  for(ix in ixM){
    MPopF$mosy[[ix]] <<- MBITESi(MPopF$mosy[[ix]],P)
  }

}
