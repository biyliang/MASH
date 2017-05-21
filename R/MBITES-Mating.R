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
#   MBITES Female Mating Routines
#   MBITES-Swarming.R defines male behaviors
#   R version
#   Sean Wu
#   January 27, 2017
#
#################################################################

# chooseMate: focal mosquito M chooses a mate
chooseMate <- function(M,P){

  # if matingQ empty fly away
  if(all(is.na(LANDSCAPE$swarmSites[[M$ix]]$matingQ$id))){
    M$lspot = "l"
    return(M)
  }

  if(rbinom(1,1,P$M.s)){
    ixM = sample(x = na.omit(LANDSCAPE$swarmSites[[M$ix]]$matingQ$id), size = 1)
    M$sire = ixM
    M$mated = TRUE
    M$stateNew = sample(x = c("F","D"), size = 1, prob = c(P$M.s,1-P$M.s))
  }

  return(M)
}
