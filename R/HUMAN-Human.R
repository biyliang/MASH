#################################################################
#
#   MASH/MBITES
#   Humans object generation routines
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################

makeHuman <- function(id, hhIx, w, bday = 0, dday = 75*365){
  #. makeHuman: Make a default human object
  list(
    myID      = id,
    hhID      = hhIx,

    Alive     = TRUE,     # alive or dead
    runIt     = "go",     # "go" or "pause" computation

    #################################################
    # General Information
    #################################################
    bDay      = bday,     # birthday
    sex       = NULL,     # male or female
    weight    = 0,        # body weight
    height    = 0,        # height

    #################################################
    # The health event queue &
    # a record of health events averted
    #################################################
    queueN     = 1,        # the number of events in the eventQ
    eventQ = replicate(n=1,expr=maxDeath,simplify=FALSE),
    averted   = NULL,
    History = list(),

    #################################################
    # Health State Objects
    #################################################
    Fever          = NULL,
    SevereDisease  = NULL,
    Nutrition      = NULL,
    Anemia         = NULL,
    Genotype       = NULL,

    #################################################
    # Malaria Parameters
    #################################################
    w = w,        # Personal biting propensity

    #################################################
    # Individual Health Behavior Parameters
    #################################################
    HealthSeeking     = NULL,
    AntiMosy          = NULL,

    # Drugs & Drug Taking OBJECT
    RX                = NULL,

    # The Mobility object
    atRisk = NULL,

    # Pathogens
    Pathogens = list(
      Pf = FALSE, #indicates current infection status
      tmInf = NULL, #vector of times of infection
      tmRec = NULL, #vector of times of recovery
      mosyID = NULL #vector of mosy ID to infect this human
    )
)}

hhAges <- function(N){
  #. hhAges: hhAges
  a = runif(1,20,40)*365
  a =  c(a, rexp(N-1,1/20/365))
  ix = which(a > 60*365)
  while(length(ix > 0)){
    a[ix] = rexp(length(ix), 1/30/365)
    ix = which(a > 60*365)
  }
  round(a)
}

makeHumans <- function(nH, hhSizes, hhIx){
  #. makeHumans: make the humans object
  k = 1
  w = rgamma(nH,1,1)
  humans = vector(mode="list",length=nH)
  for(i in 1:length(hhSizes)){
    ages = hhAges(hhSizes[i])
    for(j in 1:hhSizes[i]){
      humans[[k]] = makeHuman(id = unlist(hhIx)[k],hhIx = i,w = w[k],bday = -ages[j])
      k = k + 1
    }
  }
  return(humans)
}


##########################################
# Event Queue: Management
##########################################

########################################################################
#
#  eventQ functions
#
#  each health event has the form event.name(d, F, PAR)
#     d:    the time now
#
#     PAR:  a parameter or parameter list
#
#     F:    a function of the form F(hIx, d, PAR) that
#           executes the event
#
########################################################################

addEvent2Q <- function(ixH,event){
  #. addEvent2Q: Add single event to human associated with index ixH's event queue
  thisQ = HUMANS[[ixH]]$eventQ
  NN = HUMANS[[ixH]]$queueN
  tt = sapply(1:NN, getTfromQ, evQ=thisQ)
  HUMANS[[ixH]]$queueN <<- NN+1
  ix = which(event$t>tt)
  if(length(ix)==0) { #if event occured at event$t <= tt; advance to front of queue
    HUMANS[[ixH]]$eventQ <<- c(list(event), thisQ)
  } else { #if event will occur at event$t > tt; place at proper position in queue
    ixn = c(1:NN)[-ix]
    HUMANS[[ixH]]$eventQ <<- c(thisQ[ix], list(event), thisQ[ixn])
  }
}

rmFirstEventFromQ <- function(ixH){
  #. rmFirstEventFromQ: Remove first event from focal human ixH's event queue
  HUMANS[[ixH]]$eventQ <<- HUMANS[[ixH]]$eventQ[-1]
  HUMANS[[ixH]]$queueN <<- HUMANS[[ixH]]$queueN - 1
}

oneEvent <- function(ixH, tPause){
  #. oneEvent: run first event in focal human ixH's event queue
  with(HUMANS[[ixH]]$eventQ[[1]],{
    F(ixH = ixH, t = t, PAR = PAR)
  })
  rmFirstEventFromQ(ixH)
}

#Return the time of evQ
getTfromQ <- function(i, evQ){
  return(evQ[[i]]$t)
}


##########################################
# Event Queue: Events
##########################################

add2Q_death <- function(t, PAR = NULL){
  #. add2Q_death: Add death event to event queue
  list(t = t, PAR = PAR, F = eventDeath)
}

eventDeath <- function(ixH, t, PAR=NULL){
  #. eventDeath: death event
  HUMANS[[ixH]]$Alive <<- FALSE
}

#maximum death event (if H alive at t=75*365 kill immediately)
maxDeath <- list(t = 75*365, PAR = list(cause="oldage"), F = eventDeath, tag = "killHuman")


##########################################
# Event Queue: Simulation
##########################################

liveLife <- function(ixH, tPause){
  #. liveLife: Will execute all events in focal human H's eventQ while they are still alive and event time < tPause
  while(HUMANS[[ixH]]$Alive & HUMANS[[ixH]]$event[[1]]$t < tPause){
    oneEvent(ixH, tPause)
    if(HUMANS[[ixH]]$queueN == 0){ #if event queue empties break out of loop
      break()
    }
  }
}

runHumanEventQ <- function(tPause){
  #. runHumanEventQ: Will execute event queues for all humans for specific tPause
  for(ixH in 1:length(HUMANS)){
    liveLife(ixH, tPause)
  }
}
