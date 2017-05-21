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
#   Choosehost routines
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################


#chooseHost: mosquito host choosing routine
#M: focal mosquito
chooseHost <- function(M){

  if(M$inPointSet != "f"){ #check M is in a feeding site
    stop(paste0("chooseHost error; mosy ",M$id," inPointSet: ",M$inPointSet," , not in a feeding site"))
  }

  with(LANDSCAPE$feedSites[[M$ix]]$riskList,{

    pTmIx = pTm
    whoIx = who
    wIx = w

    #non-human hosts
    if(nO > 0){
      for(i in 1:nO){
        pTmIx = c(pTmIx,1)
        whoIx = c(whoIx,other[[i]]$typeID)
        wIx = c(wIx,other[[i]]$w)
      }
    }

    M$hostID = sample(x = whoIx,size = 1,prob = wIx*pTmIx) #select a host
    return(M)
  })

}
