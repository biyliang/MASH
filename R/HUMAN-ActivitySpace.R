#################################################################
#
#   MASH/MBITES
#   Activity Space generation and simulation routines
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################


###########################################
# initialize myTimeAtRisk object
###########################################

#init_myTimeAtRiskIx: initialize myTimeAtRisk object for single human
#hIx: id of human
#nDaily: average number of other sites visited
#move: if FALSE humans do not leave their home sites
init_myTimeAtRiskIx <- function(hIx, nDaily, move){

  if(move){ #humans move

    Nplaces = 1 + rpois(1,3)  #number of places visited
    Nplaces = min(Nplaces, LANDSCAPE$nF-1) #for very small landscapes Nplaces > nF is possible
    myPAR = list(
      nDaily = nDaily,
      Nplaces = Nplaces,
      p = rbeta(1,100,6), #parameter used to weight proportion of time spent at home
      loc = sample(x = (1:LANDSCAPE$nF)[-HUMANS[[hIx]]$hhID],size = Nplaces) #vector of other sites on the landscape to visit
    )

  } else { #humans do not move

    myPAR = list(
      nDaily = nDaily,
      Nplaces = 0,  #number of places visited
      p = 1, #parameter used to weight proportion of time spent at home
      loc = NULL
    )

  }

  HUMANS[[hIx]]$atRisk <<- myPAR
}



#init_myTimeAtRisk: initialize myTimeAtRisk for all humans
#nDaily: average number of other sites visited
#move: if FALSE humans do not leave their home sites
init_myTimeAtRisk <- function(nDaily = 1.4, move = TRUE){
  for(i in 1:LANDSCAPE$nH){
    init_myTimeAtRiskIx(hIx = i,nDaily = nDaily,move = move)
  }
}


###########################################
# simulate activity space
###########################################

#activitySpaceIx: simulate activity space for single human
#hIx: id of human
activitySpaceIx <- function(hIx){

  with(HUMANS[[hIx]],{
    with(atRisk,{

      #add risk at home site
      pD = rbeta(1,100,100*(1-p)/p) #proportion of time at home site today
      add_riskListIx(fIx = hhID,hIx = myID,pTm = pD,w = w) #add home site risk to LANDSCAPE

      #add risk from visited sites
      nD = min(rpois(1,nDaily),Nplaces) #draw random number of other sites visited
      if(nD > 0){
        fD = sample(x = Nplaces, size = nD, replace = FALSE) #sample sites in loc vector to visits
        for(i in 1:nD){
          add_riskListIx(fIx = loc[fD[i]],hIx = myID,pTm = (1-pD)/nD,w = w)
        }
      }

    })
  })

}

#activitySpace: simulate activity space for all humans
activitySpace <- function(){

  for(ix in 1:LANDSCAPE$nF){ #zero out all riskLists in LANDSCAPE
    LANDSCAPE$feedSites[[ix]]$riskList$N   <<- (LANDSCAPE$feedSites[[ix]]$riskList$N)*0
    LANDSCAPE$feedSites[[ix]]$riskList$who <<- (LANDSCAPE$feedSites[[ix]]$riskList$who)*0
    LANDSCAPE$feedSites[[ix]]$riskList$pTm <<- (LANDSCAPE$feedSites[[ix]]$riskList$pTm)*0
    LANDSCAPE$feedSites[[ix]]$riskList$w   <<- (LANDSCAPE$feedSites[[ix]]$riskList$w)*0
  }

  for(jx in 1:LANDSCAPE$nH){ #run the activity space for each human
    activitySpaceIx(hIx = jx)
  }
}
