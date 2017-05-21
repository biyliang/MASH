#################################################################
#
#   MASH/MBITES
#   PfPedigree
#   R version
#   Sean Wu
#   March 13, 2017
#
#################################################################


############################################################################
# Generic PfTransmission management functions
############################################################################

#' Initialize PfTransmission to Track Pathogen Transmission Networks
#'
#' Generate a text connection to write .csv transmission chains. This connection \strong{must} be named \code{PfTransmissionCon}.
#' See \code{\link{trackPfTransmission}}.
#'
#' @param directory directory that will be created; files will be put in directory/OUTPUT/..
#' @param fileName name of the file to write to; directory/OUTPUT/fileName.csv
#' @return a text connection
#' @examples
#' PfTransmission_init(directory, fileName)
PfTransmission_init <- function(directory, fileName){
  if(!dir.exists(paste0(directory,"OUTPUT"))){
    dir.create(paste0(directory))
    dir.create(paste0(directory,"OUTPUT"))
  }
  if(file.exists(paste0(directory,"OUTPUT/",fileName))){
    stop("PfTransmission_init cannot init a file that already exists!")
  }
  PfTransmission_TRACK <<- TRUE
  conPf = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  writeLines(text = paste0(c("time","M2H","H2M","ixH","ixS","ixM","PfID"),collapse = ","),con = conPf,sep = "\n") # write header
  return(conPf)
}

#' Add New Transmisson Event to PfTransmission
#'
#' Add a new transmission event PfTransmission .csv for data logging. If the global flag \code{PfTransmission_TRACK} is set to TRUE, this
#' function will be called during \code{\link{infectiousBite_PfSI}} (vector to human transmission) and \code{\link{infectMosquito_PfSI}} (human to vector transmission).
#'
#' @param M2H logical: direction of transmission
#' @param tBite time of infectious bite (mosquito to human)
#' @param ixH index of human host
#' @param ixS site of infectious bite
#' @param ixM index of mosquito vector
#' @param PfM Pf object from infectious mosquito
#' @return Write the following Pf pedigree to file
#' * tStart: time of start of infection
#' * tBite: time of infectious bite
#' * tMosy: time mosquito initially infected
#' * ixH: index of human host
#' * ixS: site of vector to human transmission
#' * ixM: index of mosquito vector
#' * damID: female gametocyte ID
#' * sireID: male gametocyte ID
#' * PfID: the Pf ID; every infection that makes it to bloodstream stage is considered a new clonal variant
#' @md
trackPfTransmission <- function(M2H, tBite, ixH, ixS, ixM, PfM){
  if(M2H){ # vector to human
    txtOut = paste0(c(tBite, 1L, 0L, ixH, ixS, ixM, PfM$pfid),collapse = ",")
    writeLines(text = txtOut,con = .GlobalEnv$PfTransmissionCon,sep = "\n")
  } else { # human to vector
    txtOut = paste0(c(tBite, 0L, 1L, ixH, ixS, ixM, PfM$pfid),collapse = ",")
    writeLines(text = txtOut,con = .GlobalEnv$PfTransmissionCon,sep = "\n")
  }
}

#' Import PfTransmission Data from .csv
#'
#' Import PfTransmission data written to .csv from \code{\link{trackPfTransmission}}.
#'
#' @param directory directory of output; files are in directory/OUTPUT/..
#' @param fileName name of the file; directory/OUTPUT/fileName.csv
#' @return a data frame
#' @examples
#' importPfTransmission(directory, fileName)
importPfTransmission <- function(directory, fileName){
  read.csv(file = paste0(directory,"OUTPUT/",fileName),header = TRUE)
}

#################################################################
#
#   Generic definitions of variables:
#
#   tStart  : when parasites emerge from the liver
#   tBite   : time of infectious bite
#   tMosy   : time mosquito vector initially infected
#   ixH     : index of human host
#   ixM     : index of mosquito vector
#   ixS     : site of infectious bite
#   PfID    : ID of new Pf infection
#
#   To initialize PfPedigree we must initialize the following:
#   initPfPedigree(n, PfPedigree_XX)
#   makePfPedigree = makePfPedigree_XX
#   makePfM = makePfM_XX
#   getPfParent = getPfParent_XX (here XX is one of SI, States, etc)
#
#################################################################

############################################################################
# Generic PfPedigree management functions
############################################################################

#' Initialize PfPedigree to Generation of New Clonal Variants
#'
#' Generate a text connection to write .csv pedigrees. This connection \strong{must} be named \code{PfPedigreeCon}.
#'
#' @param directory directory that will be created; files will be put in directory/OUTPUT/..
#' @param fileName name of the file to write to; directory/OUTPUT/fileName.csv
#' @return a text connection
#' @examples
#' PfPedigree_init(directory, fileName)
PfPedigree_init <- function(directory, fileName){
  if(!dir.exists(paste0(directory,"OUTPUT"))){
    dir.create(paste0(directory))
    dir.create(paste0(directory,"OUTPUT"))
  }
  if(file.exists(paste0(directory,"OUTPUT/",fileName))){
    stop("PfPedigree_init cannot init a file that already exists!")
  }
  PfPedigree_TRACK <<- TRUE
  conPf = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
  writeLines(text = paste0(c("tStart","tBite","tMosy","ixH","ixS","ixM","damID","sireID","PfID"),collapse = ","),con = conPf,sep = "\n") # write header
  return(conPf)
}

#' Add new Mosquito to Human Pf Infection to Pedigree
#'
#' Add a new bloodstream stage infection (new clonal variant) to the PfPedigree .csv for data logging.
#'
#' @param tStart time of start of infection
#' @param tBite time of infectious bite (mosquito to human)
#' @param ixH index of human host
#' @param ixM index of mosquito vector
#' @param ixS site of infectious bite
#' @param PfM Pf object from infectious mosquito
#' @return Write the following Pf pedigree to file
#' * tStart: time of start of infection
#' * tBite: time of infectious bite
#' * tMosy: time mosquito initially infected
#' * ixH: index of human host
#' * ixS: site of vector to human transmission
#' * ixM: index of mosquito vector
#' * damID: female gametocyte ID
#' * sireID: male gametocyte ID
#' * PfID: the Pf ID; every infection that makes it to bloodstream stage is considered a new clonal variant
#' @md
addPf2Pedigree <- function(tStart, tBite, ixH, ixM, ixS, PfM){
  with(PfM,{
    txtOut = paste0(c(tStart,tBite,tm,ixH,ixS,ixM,damID,sireID,pfid),collapse = ",")
    writeLines(text = txtOut,con = .GlobalEnv$PfPedigreeCon,sep = "\n")
  })
}

#' Import PfPedigree Data from .csv
#'
#' Import PfPedigree data written to .csv from \code{\link{addPf2Pedigree}}.
#'
#' @param directory directory of output; files are in directory/OUTPUT/..
#' @param fileName name of the file; directory/OUTPUT/fileName.csv
#' @return a data frame
#' @examples
#' importPfTransmission(directory, fileName)
importPfPedigree <- function(directory, fileName){
  read.csv(file = paste0(directory,"OUTPUT/",fileName),header = TRUE)
}


############################################################################
# getPfParent_XX functions:
#
# for PfSI module this will only return a single pfid but in general
# the damID and sireID could be a vector of male and female gametocyte IDs
# ingested by the mosquito during host to vector transmission during a
# blood meal.
#
############################################################################

#' Return a pfid During Human to Mosquito Transmission
#'
#' For PfSI module this will only return a single pfid but in general the damID and sireID could be a vector of male and female gametocyte IDs
#' ingested by the mosquito during host to vector transmission during a blood meal.
#'
#' @param ixH index of human
#' @return a single pfid (integer)
#' @examples
#' getPfParent_SI(ixH)
getPfParent_SI <- function(ixH){
  HUMANS[[ixH]]$Pathogens$Pf$pfid
}


############################################################################
# PfPedigree Modules:
#
# makePfM_XX: make a PfM object, at a minimum must include pfid, damID, and sireID
# PfPedigree_XX: the pedigree slot; contains all information for that pedigree
# makePfPedigree_XX: during infectious bite, make a new pedigree for that variant
#
############################################################################


######################
# PfPedigree: full
######################

#' Generate the Pf Object for Mosquito Stage of Parasite Life Cycle
#'
#' During successful human to mosquito transmission makePfM() generates the PfM object that is passed from human to mosquito.
#' This object contains all of the Pf information necessary for the module.
#'
#' @param ixH index of human
#' @param tBite time of bite
#' @param ixS index of site
#' @return a list of two elements
#' * infected: TRUE
#' * PfM: The PfM object
#'    * tm: time of human to vector transmission
#'    * ixS: site of human to vector transmission
#'    * ixH: site of human to vector transmission
#'    * damID: female gametocyte IDs
#'    * sireID: male gametocyte IDs
#'    * pfid: Pf ID from the infected human
#' @md
#' @examples
#' makePfM(ixH, tBite, ixS)
makePfM <- function(ixH, tBite, ixS){
  damID  = getPfParent(ixH)
  sireID = getPfParent(ixH)
  pfid = getPfParent(ixH)
  list(
    infected = TRUE,
    PfM = list(tm=tBite, ixS=ixS, ixH=ixH, damID=damID, sireID=sireID, pfid=pfid)
  )
}

#' Generate the Pf Object for Simulated Biting
#'
#' This generates the Pf object for \code{\link{simbite_PfSI}}. It does not have the element \code{infected} in the returned list
#' because this is for vector to human transmission; therefore \code{spz} must be set to \code{TRUE}.
#'
#' @param ixH index of human
#' @param tBite time of bite
#' @param ixS NULL; index of site
#' @return a list of two elements
#' * spz: TRUE
#' * PfM: The PfM object
#'    * tm: time of human to vector transmission
#'    * ixS: NULL
#'    * ixH: site of human to vector transmission
#'    * damID: NULL
#'    * sireID: NULL
#'    * pfid: NULL
#' @md
#' @examples
#' makePfM(ixH, tBite, ixS)
makePf0 <- function(ixH, tBite, ixS = NULL){
  PfM = list()
  PfM[[1]] = list(tm=tBite, ixS=ixS, ixH=ixH, damID=NULL, sireID=NULL, pfid=NULL)
  list(
    spz = TRUE,
    PfM = PfM
  )
}
