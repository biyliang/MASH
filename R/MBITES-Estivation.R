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
#   MBITES Estivation Routines
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################

##################################################################
# ESTIVATION
#
#   Note that survival through estivation
#   is determined in the queue estivation function
#   and if a mosquito would not have survived,
#   it is immediately killed.
#
#   The "estivation bout" simply changes
#   bState to "F"
#
##################################################################

queueEstivation <- function(M,P){
  with(P,{

    if(ESTIVATION && isAlive(M)){
      if(rbinom(n=1,size=1,prob=prEstivate(M$tNow,Emax,Eb))){
        if(rbinom(n=1,size=1,prob=E.p)){
          M$stateNew = "E"
          M$tNext = wakeUpTime(eEndm,eEndV)
        } else {
          M$stateNew = "D"
        }
      }
    }

    return(M)
  })
}

prEstivate <- function(tNow, Emax, Eb){
  pmax(0, cos(2*pi*(tNow-Emax)/365)-Eb)
}

wakeUpTime <- function(eEndm, eEndV){
  rnorm(n = 1, mean = eEndm, sd = eEndV)
}
