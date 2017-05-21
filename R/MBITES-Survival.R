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
#   Survival
#   Senescence, flight and resting survival
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################


##########################################
# Flight Survival
##########################################

# getSFp: get baseline survival probability
getSFp <- function(M,P){
  with(P,{
    switch(M$state,
      F = F.p,
      B = B.p,
      R = R.p,
      L = L.p,
      O = O.p,
      M = M.p,
      S = S.p,
      E = 1
    )
  })
}

# surviveFlight:
surviveFlight <- function(M,P){
  if(isActive(M)){
    p = getSFp(M,P) # baseline survival
    if(P$TATTER){
      M$damage = M$damage + rTatterSize(P)
      p = p * pTatter(M,P)
    }
    if(P$SENESCE){
      p = p * pSenesce(M,P)
    }
    if(!rbinom(1,1,p)){
      M$stateNew = "D"
    }
  }
  return(M)
}


##########################################
# Resting Survival
##########################################

# myHaz: get landscape hazards
myHaz <- function(M){
  switch(M$inPointSet,
    f = LANDSCAPE$feedSites[[M$ix]]$haz,
    l = LANDSCAPE$aquaSites[[M$ix]]$haz,
    m = LANDSCAPE$swarmSites[[M$ix]]$haz,
    s = LANDSCAPE$sugarSites[[M$ix]]$haz
  )
}

surviveResting <- function(M,P){
  if(isActive(M)){
    if(!rbinom(1,1,myHaz(M))){
      M$stateNew = "D"
    }
  }

  return(M)
}


#############################################
#  Wing Tattering
#############################################

# rTatterSize: zero-inflated per-bout additive wing damage from tattering
rTatterSize <- function(P){
  with(P,{
    if(runif(1) > ttsz.p){
      return(0)
    } else {
      return(rbeta(1,ttsz.a,ttsz.b))
    }
  })
}

# pTatter: probability of death due to tattering
pTatter <- function(M,P){
  with(P,{
    return((2+ttr.b)/(1+ttr.b) - exp(M$damage*ttr.a)/(ttr.b + exp(M$damage*ttr.a)))
  })
}


#############################################
#  Senescence
#############################################

# pSenesce:
pSenesce <- function(M,P){
  age = M$tNow - M$bDay
  with(P,{
    return((2+sns.b)/(1+sns.b) - exp(sns.a*age)/(sns.b + exp(sns.a*age)))
  })
}
