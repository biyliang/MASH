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
#   Cohort
#   MBITES-Basic
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################


###########################################
# choose mate routine
###########################################

# chooseMate: focal mosquito M chooses a mate
chooseMate.basic <- function(M,P){

  if(rbinom(1,1,P$M.s)){
    M$sire = 1
    M$mated = TRUE
    M$stateNew = sample(x = c("F","D"), size = 1, prob = c(P$M.s,1-P$M.s))
  }

  return(M)
}

###########################################
# choose host routine
###########################################

chooseHost.basic <- function(M,P){

  if(M$inPointSet != "f"){ #check M is in a feeding site
    stop(paste0("chooseHost error; mosy ",M$id," not in a feeding site"))
  }

  if(runif(1) < P$Q){
    M$hostID = 1
  } else {
    M$hostID = -1
  }

  return(M)
}


###########################################
# Movement function
###########################################

# bMove.nomove: used for MBITESi.basic.nomove; updates inPointSet but not ix
bMove.nomove <- function(M){
  if(!M$state %in% c("F","L","S","M")){
    return(M)
  } else if(M$state == "F"){
    M$inPointSet = "f"
    return(M)
  } else if(M$state == "L"){
    M$inPointSet = "l"
    return(M)
  } else if(M$state == "S"){
    M$inPointSet = "s"
    return(M)
  } else if(M$state == "M"){
    M$inPointSet = "m"
    return(M)
  } else {
    browser("unrecognized behavioral state")
  }
}

# bMove.move: used for MBITESi.basic.move; updates inPointSet and ix
bMove.move <- function(M,MvOb){
  if(!M$state %in% c("F","L","S","M")){
    return(M)
  } else if(M$state == "F"){
    M$ix = switch(M$inPointSet,
        f = bMove.newIx(MvMat = MvOb$F2F[[M$ix]]),
        s = bMove.newIx(MvMat = MvOb$S2F[[M$ix]]),
        m = bMove.newIx(MvMat = MvOb$M2F[[M$ix]]),
        l = bMove.newIx(MvMat = MvOb$L2F[[M$ix]])
      )
    M$inPointSet = "f"
    return(M)
  } else if(M$state == "L"){
    M$ix = switch(M$inPointSet,
        f = bMove.newIx(MvMat = MvOb$F2L[[M$ix]]),
        s = bMove.newIx(MvMat = MvOb$S2L[[M$ix]]),
        m = bMove.newIx(MvMat = MvOb$M2L[[M$ix]]),
        l = bMove.newIx(MvMat = MvOb$L2L[[M$ix]])
      )
    M$inPointSet = "l"
    return(M)
  } else if(M$state == "S"){
    M$ix = switch(M$inPointSet,
        f = bMove.newIx(MvMat = MvOb$F2S[[M$ix]]),
        s = bMove.newIx(MvMat = MvOb$S2S[[M$ix]]),
        m = bMove.newIx(MvMat = MvOb$M2S[[M$ix]]),
        l = bMove.newIx(MvMat = MvOb$L2S[[M$ix]])
      )
    M$inPointSet = "s"
    return(M)
  } else if(M$state == "M"){
    M$ix = switch(M$inPointSet,
        f = bMove.newIx(MvMat = MvOb$F2M[[M$ix]]),
        s = bMove.newIx(MvMat = MvOb$S2M[[M$ix]]),
        m = bMove.newIx(MvMat = MvOb$M2M[[M$ix]]),
        l = bMove.newIx(MvMat = MvOb$L2M[[M$ix]])
      )
    M$inPointSet = "m"
    return(M)
  } else {
    browser("unrecognized behavioral state")
  }
}

# bMove.newIx: select a new index given a MvMat (element in MvOb)
bMove.newIx <- function(MvMat){
  x = runif(1)

  ixNew = with(MvMat,{
    if(x <= PR[1]){ #no movement
      ix
    } else {
      if(x <= PR[1] + PR[2]){ #near movement
        # sample(x = near$id,size = 1,prob = near$pr)
        ixNear = sample(x = 1:length(near$id),size = 1,prob = near$pr)
        near$id[ixNear]
      } else {
        if(x <= sum(PR)){ #around movement
          browser("'around' movement not yet implemented")
        } else { #moveFar movement
          browser("'moveFar' movement not yet implemented")
        }
      }
    }
  })

  return(ixNew)
}

##########################################
# Encounter Host Routines
##########################################

humanEncounter.basic <- function(M,P){
  #. humanEncounter: run encounter routine for human host
  #M: mosquito object
  #P: parameters list
  with(P,{

    if(!rbinom(n=1,size=1,prob=surviveH)){ # does not survive to probe
      M$stateNew = "D"
    } else { # survives to probe
      if(rbinom(n=1,size=1,prob=probeH)){ # undeterred
        # M = probing(M) # probe the host
        if(!rbinom(n=1,size=1,prob=surviveprobeH)){ # does not survive probing
          M$stateNew = "D"
        } else { # survives probing
          if(rbinom(n=1,size=1,prob=feedH)){ # successfully begins feeding
            M = BloodMeal(M,P) # MBITES-Energetics.R
            # M = getInfected(M) # MBITES-HostEncounter.R
            if(HISTORY){M = historyFeed(M)} # track history
            M$stateNew = "R"
          }
        }
      }
    }

  return(M)
  })
}

#zooEncounter: run encounter routine for non-human host
zooEncounter.basic <- function(M,P){
  #. zooEncounter: run encounter routine for non-human host
  #M: mosquito object
  #P: parameters list
  with(P,{

    if(!rbinom(n=1,size=1,prob=surviveL)){ # does not survive encounter
      M$stateNew = "D"
    } else { # survives encounter
      if(rbinom(n=1,size=1,prob=feedL)){ # successful feed
        M=BloodMeal(M,P) # MBITES-Energetics.R
        if(HISTORY){M = historyFeed(M)} # track history
        M$stateNew="R"
      }
    }

  return(M)
  })
}


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

#enterHouse.basic: mosquito attempts to enter a house and appropriate events are called
enterHouse.basic <- function(M,P){
  if(runif(1) < P$enterhouse.mean){ #mosquito found a gap
    # call EAVE.TUBE
  } else { #mosquito did not find a gap
    M$lspot = newSpot(M,P)
    M = surviveFlight(M,P)
    if(M$lspot == "i"){
      M = Recall(M,P)
    }
  }

  return(M)
}

# landingSpot.basic:
landingSpot.basic <- function(M,P){
  if(isActive(M)){
    oldSpot = M$lspot
    M$lspot = newSpot(M,P) # choose new lspot
    if(oldSpot != "i" & M$lspot == "i"){
      M = enterHouse.basic(M,P) # enterHouse
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


##########################################
# Resting Survival
##########################################

# myHaz: get landscape hazards
myHaz.basic <- function(M,P){
  switch(M$inPointSet,
    f = P$feedHaz.mean,
    l = P$aquaHaz.mean,
    m = P$swarmHaz.mean,
    s = P$sugarHaz.mean
  )
}

surviveResting.basic <- function(M,P){
  if(isActive(M)){
    if(!rbinom(1,1,myHaz.basic(M,P))){
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
#   2. bMove: movement between point classes (if needed)
#   3. boutFun: run bout function
#   4. run energetics and check if alive
#   5. run landingSpot and check if alive
#   6. run surviveResting/surviveFlight and check if alive
#   7. update tNext
#   8. update state to stateNew which is determined in boutGeneric
#
################################################################################

# boutBasic updates tNext and stateNew without point movement
boutBasic.nomove <- function(M,P,boutFun,...){

  # update time and state
  M$tNow = M$tNext # update time
  M$state = M$stateNew # update current state
  M = timingExponential(M,P) # update tNext

  # movement
  M = bMove.nomove(M)

  # bout
  M = boutFun(M,P,...)

  # energetics
  M = energetics(M,P) # MBITES-Energetics.R

  # landing spot
  M = landingSpot.basic(M,P) # MBITES-Basic.R

  # survival
  M = surviveResting.basic(M,P) # MBITES-Basic.R
  M = surviveFlight(M,P) # MBITES-Survival.R

  # check queueing
  M = queueEstivation(M,P) # MBITES-Estivate.R

  # log history
  if(P$HISTORY){M = historyTrack(M)}

  return(M)
}

# boutBasic.move updates tNext and stateNew with point movement
boutBasic.move <- function(M,P,MvOb,boutFun,...){

  # update time and state
  M$tNow = M$tNext # update time
  M$state = M$stateNew # update current state
  M = timingExponential(M,P) # update tNext
  curState = M$state
  # movement
  M = bMove.move(M,MvOb)

  # bout
  M = boutFun(M,P,...)

  # energetics
  M = energetics(M,P) # MBITES-Energetics.R

  # landing spot
  M = landingSpot.basic(M,P) # MBITES-Basic.R

  # survival
  M = surviveResting.basic(M,P) # MBITES-Basic.R
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

boutF.basic <- function(M,P){
  if(M$lspot != "l" && rbinom(n=1,size=1,prob=P$F.s)){
    M$stateNew = "B"
  }

  return(M)
}


#################################################################
#  Blood Feeding Attempt Bout :: B
#################################################################

boutB.basic <- function(M,P){
  if(rbinom(1,1,P$B.s)){
    M = chooseHost.basic(M,P) # MBITES-ChooseHost.R
  } else {
    M$hostID = 0
  }

  if(M$hostID > 0){
    M = humanEncounter.basic(M,P) # MBITES-Basic.R
  } else if(M$hostID == -1){
    M = zooEncounter.basic(M,P) # MBITES-Basic.R
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

boutR.basic <- function(M,P){
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

boutL.basic <- function(M,P){
  if(rbinom(1,1,P$L.s)){
    M$stateNew = "O"
  }

  return(M)
}


#################################################################
#  Egg Laying Attempt Bout :: O
#################################################################

layEggs.basic <- function(M,P){
  if(rbinom(1,1,P$O.s)){
    M$batch = 0
    M$stateNew = "F"
  }

  return(M)
}

boutO.basic <- function(M,P){
  if(isAlive(M)){
    M = layEggs.basic(M,P)
  }
  return(M)
}


#################################################################
#  Sugar Feeding Attempt Bout :: O
#################################################################

boutS.basic <- function(M,P){
  with(P,{

    # female behavior
    M$stateNew = sample(x = c("F","D"), size = 1, prob = c(S.s,1-S.s))
    if(isAlive(M)){
      M$energy = 1

      if(!M$mature){
        M$energyPreG = M$energyPreG - preGblood
        if(M$energyPreG <= 0){
          M$mature = TRUE
        }
      }
    }

    return(M)
  })
}


#################################################################
#  Female Mating Bout :: M
#################################################################

boutMF.basic <- function(M,P){
  M = chooseMate.basic(M,P) # MBITES-Mating.R
  return(M)
}


#################################################################
#  Estivation Bout :: E
#################################################################

boutE.basic <- function(M,P){
  M$stateNew = "F"
  return(M)
}


#############################################################################
#   MBITES-Basic
#   Female Mosquito
#   Mosquito Bout-based & Individual-based Transmission Ecology Simulation
#############################################################################

# MBITESi: run M-BITES model for a single mosquito
# tMax: gloabl time tick; ie all the mosy start born at time=0 and ticks start 1,2,3,...
# run bouts while tNext (time of next bout/action) < global tMax time tick
# run bouts while mosquito is still alive
MBITESi.basic.nomove <- function(M,P){

  while(M$stateNew != "D"){
    M = switch(M$stateNew,
      F = boutBasic.nomove(M,P,boutF.basic),
      B = boutBasic.nomove(M,P,boutB.basic),
      R = boutBasic.nomove(M,P,boutR.basic),
      L = boutBasic.nomove(M,P,boutL.basic),
      O = boutBasic.nomove(M,P,boutO.basic),
      M = boutBasic.nomove(M,P,boutMF.basic),
      S = boutBasic.nomove(M,P,boutS.basic),
      E = boutBasic.nomove(M,P,boutE.basic)
    )
  }
  return(M)
}

# MBITESi.basicVec: vectorized version over M
MBITESi.basicVec.nomove <- Vectorize(FUN = MBITESi.basic.nomove,vectorize.args = "M",SIMPLIFY = FALSE)

# MBITESi: run M-BITES model for a single mosquito with movement
# tMax: gloabl time tick; ie all the mosy start born at time=0 and ticks start 1,2,3,...
# run bouts while tNext (time of next bout/action) < global tMax time tick
# run bouts while mosquito is still alive
MBITESi.basic.move <- function(M,P,MvOb){

  while(M$stateNew != "D"){
    M = switch(M$stateNew,
      F = boutBasic.move(M,P,MvOb,boutF.basic),
      B = boutBasic.move(M,P,MvOb,boutB.basic),
      R = boutBasic.move(M,P,MvOb,boutR.basic),
      L = boutBasic.move(M,P,MvOb,boutL.basic),
      O = boutBasic.move(M,P,MvOb,boutO.basic),
      M = boutBasic.move(M,P,MvOb,boutMF.basic),
      S = boutBasic.move(M,P,MvOb,boutS.basic),
      E = boutBasic.move(M,P,MvOb,boutE.basic)
    )
  }
  return(M)
}

# MBITESi.basicVec: vectorized version over M
MBITESi.basicVec.move <- Vectorize(FUN = MBITESi.basic.move,vectorize.args = "M",SIMPLIFY = FALSE)

# MBITES.basic: run MBITES-Basic for a cohort of mosquitoes
# N: size of cohort
# P: M-BITES parameters
MBITES.basic <- function(N, P, parallel = TRUE, move = FALSE, ...){

  mosy = vector(mode="list",length=N)
  for(i in 1:N){
    mosy[[i]] = makeMosquito(id = i,ix = 1,tm = 1,female = TRUE)
  }

  if(parallel){
    if(move){
      out = parallel::pvec(v = mosy,FUN = MBITESi.basicVec.move,P = P,MvOb = MvOb,...)
    } else {
      out = parallel::pvec(v = mosy,FUN = MBITESi.basicVec.nomove,P = P,...)
    }
  } else {
    if(move){
      out = lapply(X = mosy,FUN = MBITESi.basic.move,P=P,MvOb = MvOb,...)
    } else {
      out = lapply(X = mosy,FUN = MBITESi.basic.nomove,P=P,...)
    }
  }

  return(out)
}
