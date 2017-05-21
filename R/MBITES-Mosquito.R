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
#   Mosquito object
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################


##########################################################
#  Define the mosquito object
##########################################################

makeMosquito <- function(id=0,tm=NULL,ix=NULL,state="M",inPointSet="l",EIP=12,
                         mature=FALSE,energyPreG=0,female=TRUE,history=TRUE){

  M = list(

    # ID and time
    id        = id,        # mosquito id
    bDay      = tm,        # time of emergence
    tNow      = tm,        # time of last event
    tNext     = tm,        # time to next event

    # State and Location
    state      = state,       # {F,B,R,L,O,S,M,E,D}
    stateNew   = state,       # {F,B,R,L,O,S,M,E,D}
    inPointSet = inPointSet,  # class of site {f,l,s,m}
    ix         = ix,          # index of site
    mature     = mature,      # mature

    # Other State Variables
    lspot     = "l",      # landing spot (i: inside wall, w: outside wall, v: outside vegetation, r: feed, l: leave)
    damage    = 0,        # wing tattering
    energy    = 1,        # energy reserves
    female    = female    # sex of mosquito
  )

  if(female){
    # Egg Batch Variables
    M$bmSize = 0         # the size of the blood meal, relative to max
    M$batch  = 0         # female eggs in batch
    M$eggT   = 0         # the minimum time before eggs are mature
    M$eggP   = 0         # the mimimum provision for eggs to mature

    # Variables for Mosquito Mating
    M$mated       = FALSE
    M$sire        = 0
    M$energyPreG  = energyPreG  # pre-gonotrophic energy requirement

    # Infection events
    M$hostID  = 0          # the id of the host: -1::none; 0::not human
    # M$spz     = 0          # presence/absence of sporozoites
    M$EIP     = EIP        # presence/absence of sporozoites
    M$Pf      = list(n=0, t=0, spz=FALSE, PfM = list())
  }

  if(history){
    M$history = makeHistory(M)
    M$bionomics = NULL
  }

  return(M)
}

makeHistory <- function(M){
  if(M$female){
    list(
      stateH     = NULL,  # state trajectory
      timeH      = NULL,      # transition times
      ixH        = M$ix,     # sites visited
      pSetH      = M$inPointSet,    # point sets visited

      feedAllH   = 0,      # number of blood meals
      feedAllT   = NULL,   # times of blood meals
      feedHumanH = 0,      # number of blood meals on human hosts
      feedHumanT = NULL,   # times of blood meals on human hosts
      feedIxH    = NULL,   # ids of all blood hosts

      bmSizeH    = NULL,
      batchH     = NULL
      )
  } else {
    list(
      stateH     = NULL,  # state trajectory
      timeH      = NULL,      # transition times
      ixH        = M$ix,     # sites visited
      pSetH      = M$inPointSet    # point sets visited
      )
  }
}

##########################################################
#  Make a New Mosquito Population
##########################################################

#makeMosquitoesFemale: make a new mosquito population
#N: size of cohort
#female: boolean
#t: bDay/tNow/tNext set to t
#offset: number of null memory slots at end of list
#...: additional named parameters for makeMosquito()
makeMosquitoCohort <- function(N, female, tm = 0, offset = NULL, ...){

  if(is.null(offset)){
    offset = floor(N * 1.5)
  }

  mosy = vector(mode="list",length=N+offset)
  if(!exists("LANDSCAPE")){
    message("LANDSCAPE not defined; initial ix set to 1 for all mosquitoes")
    initIx = rep(1,length=N)
  } else {
    initIx = sample(x = LANDSCAPE$nA,size = N,replace = TRUE)
  }

  for(i in 1:N){ #fill up initial cohort
    # mosy[[i]] = makeMosquito(id = i,ix = initIx[i],tm = tm,female = female, ...)
    mosy[[i]] = makeMosquito(id = paste0(tm,"_",i),ix = initIx[i],tm = tm,female = female, ...)
  }
  for(i in (N+1):(N+offset)){ #prealloc memory slots
    mosy[[i]] = makeMosquito(id = NULL,ix = NULL,female = female, ...)
  }

  nullIx = which(sapply(mosy,FUN = function(x){ #get NULL ix
    is.null(x$id)
  }))

  return(list(mosy=mosy,nullIx=nullIx))

}


##########################################################
#  Manage Mosquito Population
##########################################################

# enlargeMosyPop: enlarge a mosquito population
# offset: number of additional slots to add; if NULL add 1.5 times the current size of list
enlargeMosyPop <- function(female = TRUE, offset = NULL){

  # enlarge the female population
  if(female){
    if(is.null(offset)){
      offset = floor(length(MPopF$mosy) * 1.5)

      N = length(MPopF$mosy)
      for(i in (N+1):(N+offset)){
        MPopF$mosy[[i]] <<- makeMosquito(id = NULL,ix = NULL,female = female)
      }
      MPopF$nullIx <<- which(sapply(MPopF$mosy,FUN = function(x){ #get NULL ix
        is.null(x$id)
      }))

    }
  # enlarge the male population
  } else {
    if(is.null(offset)){
      offset = floor(length(MPopM$mosy) * 1.5)
    }

    N = length(MPopM$mosy)
    for(i in (N+1):(N+offset)){
      MPopM$mosy[[i]] <<- makeMosquito(id = NULL,ix = NULL,female = female)
    }
    MPopM$nullIx <<- which(sapply(MPopM$mosy,FUN = function(x){ #get NULL ix
      is.null(x$id)
    }))

  }

}


# resetMosyPop: after outputting dead mosquito histories to external file, find dead ix, replace and add to nullIx to be reused
resetMosyPop <- function(female){

  # reset the female population
  if(female){
    deadIx = which(sapply(MPopF$mosy,Negate(isAlive)))
    for(ix in deadIx){
      MPopF$mosy[[ix]] <<- makeMosquito(id = NULL,ix = NULL,female = female)
    }
    MPopF$nullIx <<- which(sapply(MPopF$mosy,FUN = function(x){ #get NULL ix
      is.null(x$id)
    }))

  # reset the male population
  } else {
    deadIx = which(sapply(MPopM$mosy,Negate(isAlive)))
    for(ix in deadIx){
      MPopM$mosy[[ix]] <<- makeMosquito(id = NULL,ix = NULL,female = female)
    }
    MPopM$nullIx <<- which(sapply(MPopM$mosy,FUN = function(x){ #get NULL ix
      is.null(x$id)
    }))

  }

}
