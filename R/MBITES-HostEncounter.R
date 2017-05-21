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
#   MBITES generic host encounter routines
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################


##########################################
# Encounter Host Routines
##########################################

humanEncounter <- function(M,P){
  #. humanEncounter: run encounter routine for human host
  #M: mosquito object
  #P: parameters list
  with(P,{

    if(!rbinom(n=1,size=1,prob=surviveH)){ # does not survive to probe
      M$stateNew = "D"
    } else { # survives to probe
      if(rbinom(n=1,size=1,prob=probeH)){ # undeterred
        M = probing(M) # probe the host
        if(!rbinom(n=1,size=1,prob=surviveprobeH)){ # does not survive probing
          M$stateNew = "D"
        } else { # survives probing
          if(rbinom(n=1,size=1,prob=feedH)){ # successfully begins feeding
            M = BloodMeal(M,P) # MBITES-Energetics.R
            M = getInfected(M) # MBITES-HostEncounter.R
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
zooEncounter <- function(M,P){
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

# nullEncounter: encounter for null host
nullEncounter <- function(M,P){
  #. nullEncounter: encounter for null host
  return(M)
}


##########################################
# Probing & Infection
##########################################

# probeHost & infectMosquito are generic functions that are defined in the appropriate PATHOGEN module

# probing: probe a host
probing <- function(M){
  if(M$Pf$n > 0){ # if mosquito has at least one infection
    M = updateMosqInf(M)
    with(M,{
      # probeHost(ixH = hostID, tBite = tNow, ixS = ix, ixM = id, Pf = Pf) # in appropriate PATHOGEN module
      probeHost(tBite = tNow, ixH = hostID, ixS = ix, ixM = id, Pf = Pf)
      return(M)
    })
  } else { # if no infections immediately return
    return(M)
  }
}

#updateMosqInf: update infection status of infected/infectious mosquitoes
updateMosqInf <- function(M){
  for(n in 1:M$Pf$n){
    if(M$tNow > (M$Pf$t[n] + M$EIP)){ #if past EIP has sporozoites; if not past EIP no sporozoites
      M$Pf$spz[n] = TRUE
    }
  }
  return(M)
}

# getInfected: infectMosquito defined in PfSI.R
getInfected <- function(M){
  im = infectMosquito(tBite = M$tNow, ixH = M$hostID, ixS = M$ix, ixM = M$id) # in appropriate PATHOGEN module
  if(im$infected==TRUE){ #if mosy got infected
    #update Pf list of mosy
    M$Pf$n                = M$Pf$n+1
    M$Pf$t[M$Pf$n]        = M$tNow
    M$Pf$spz[M$Pf$n]      = FALSE # however sporozoites are not produced by this clonal variant until passes EIP
    M$Pf$PfM[[M$Pf$n]]    = im$PfM
  }
  return(M)
}


##########################################
# History tracking
##########################################

historyFeed <- function(M){

  if(M$hostID > 0){ # human host
    M$history$feedAllH   = M$history$feedAllH + 1           # number of blood meals
    M$history$feedAllT   = c(M$history$feedAllT,M$tNow)     # times of blood meals
    M$history$feedHumanH = M$history$feedHumanH + 1         # number of blood meals on human hosts
    M$history$feedHumanT = c(M$history$feedHumanT,M$tNow)   # times of blood meals on human hosts
    M$history$feedIxH    = c(M$history$feedIxH,M$hostID)    # ids of all blood hosts
    M$history$bmSizeH    = c(M$history$bmSizeH,M$bmSize)    # size of blood meal
    M$history$batchH     = c(M$history$batchH,M$batch)      # size of egg batch
  } else { # animal host
    M$history$feedAllH   = M$history$feedAllH + 1         # number of blood meals
    M$history$feedAllT   = c(M$history$feedAllT,M$tNow)   # times of blood meals
    M$history$feedIxH    = c(M$history$feedIxH,M$hostID)  # ids of all blood hosts
    M$history$bmSizeH    = c(M$history$bmSizeH,M$bmSize)  # size of blood meal
    M$history$batchH     = c(M$history$batchH,M$batch)    # size of egg batch
  }

  return(M)
}
