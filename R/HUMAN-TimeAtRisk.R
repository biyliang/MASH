#################################################################
#
#   MASH/MBITES
#   Risklist generation and simulation routines
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################


########################################
# Initialize riskList for sites
########################################

#init_riskListIx: sets up riskList object for a single site
init_riskListIx <- function(maxH = 10){
  list(
    maxH  = maxH, #maximum number of humans in this atRiskQ
    N     = 0,
    who   = rep(0, maxH), #who: ix of hosts at risk
    pTm   = rep(0, maxH), #pTm: person time at risk
    w     = rep(0, maxH), #w: biting weight on hosts
    nO = 1, #nO: number of other hosts
    other = list(
      zoo = list(
        w = rgamma(n = 1,shape = 1,rate = 1), #biting weight on livestock
        typeID = -1 #livestock host id
      )
    )
  )
}

#init_riskList: set up daily riskList for all sites
init_riskList <- function(){

  for(i in 1:LANDSCAPE$nF){
    LANDSCAPE$feedSites[[i]]$riskList <<- init_riskListIx()
  }

}


########################################
# Manage risk object size
########################################

#lengthen_riskListIx: lengthen a single daily risk list by single human slot
lengthen_riskListIx <- function(ix){
  with(LANDSCAPE$feedSites[[ix]],{
    riskList$who  = c(riskList$who,0)
    riskList$pTm  = c(riskList$pTm,0)
    riskList$w    = c(riskList$w,0)
    riskList$maxH = riskList$maxH + 1
    riskList
  })
}

#lengthen_riskListIx: lengthen all daily risk lists by a single human slot
lengthen_riskList <- function(){
  for(ix in 1:LANDSCAPE$nF){
    LANDSCAPE$feedSites[[ix]]$riskList <<- lengthen_riskListIx(ix = ix)
  }
}

#add_riskListIx: add risk information for a single human for single site
#fIx: index of feeding site
#hIx: index of human
#pTm: time at risk
#w: biting weight
add_riskListIx <- function(fIx, hIx, pTm, w){

  nH = LANDSCAPE$feedSites[[fIx]]$riskList$N + 1
  if(nH > LANDSCAPE$feedSites[[fIx]]$riskList$maxH){
    lengthen_riskList()
  }
  LANDSCAPE$feedSites[[fIx]]$riskList$N       <<- nH   #update number of humans
  LANDSCAPE$feedSites[[fIx]]$riskList$pTm[nH] <<- pTm  #update time at risk
  LANDSCAPE$feedSites[[fIx]]$riskList$w[nH]   <<- w    #update biting weights
  LANDSCAPE$feedSites[[fIx]]$riskList$who[nH] <<- hIx  #update who is at risk
}
