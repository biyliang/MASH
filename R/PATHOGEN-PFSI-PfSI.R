########################################
#  PFSI
#  Sean Wu
#  April 5, 2017
########################################

###################################################################
# From infectious bite to infection
###################################################################

# probeHost_PfSI: probeHost called from probing(); defined in MBITES-HostEncounter.R
# infect a human
probeHost_PfSI <- function(tBite, ixH, ixS, ixM, Pf){
  if(any(Pf$spz)){ # sample a clonal variant if multiple
    PfClonalVar = which(Pf$spz)
    PfIx = sample(x = PfClonalVar, size = 1)
    infectiousBite_PfSI(tBite = tBite, ixH = ixH, ixS = ixS, ixM = ixM, PfM = Pf$PfM[[PfIx]])
  }
}

#' Infectious Bite (Vector to Human Transmission)
#'
#' This is the routine that handles a bite from an infectious vector on a human for PfSI module. It calls \code{link{add2Q_startPfSI}}
#' to begin queuing human infection events if tranmission of sporozoites into the bloodstream is successful (probability \code{b}).
#' Depending on the global flags \code{PfPedigree_TRACK} and \code{PfTransmission_TRACK} it may call \code{\link{addPf2Pedigree}} or \code{\link{trackPfTransmission}} to log data.
#'
#'
#' @param ixH index of human
#' @return a single pfid (integer)
#' @examples
#' getPfParent_SI(ixH)
infectiousBite_PfSI <- function(tBite, ixH, ixS, ixM, PfM){
  if(rbinom(1,1,HUMANS[[ixH]]$Pathogens$Pf$b)){
    tInfStart = tBite + ttInfectionPf() # when does latent -> infected occur
    if(PfPedigree_TRACK){
      addPf2Pedigree(tStart = tInfStart, tBite = tBite, ixH = ixH, ixS = ixS, ixM = ixM, PfM = PfM)
    }
    if(PfTransmission_TRACK){
      trackPfTransmission(M2H = TRUE, tBite = tBite, ixH = ixH, ixS = ixS, ixM = ixM, PfM = PfM)
    }
    add2Q_startPfSI(ixH, tInfStart, PfM$pfid) #FIX THIS LATER; NEW CLONAL VARIANT = NEW ID
  }
}

###################################################################
# Simulated bites
###################################################################

#' \code{add2Q_XX} Family: Add Simulated Bite to EventQ
#'
#' \code{add2Q_simbitePfSI} adds a simulated bite to human ixH's eventQ.
#' For more information on how this family of functions works, see \url{http://smitdave.github.io/MASH-Development/docHumans.html#add2q_xx}
#'
#' @param ixH index of human
#' @param t time the event will occur
#' @param PAR NULL; the function will generate a default pathogen object with \code{\link{makePf0}}
#' @examples
#' add2Q_simbitePfSI(ixH, t, PAR = NULL)
add2Q_simbitePfSI <- function(ixH, t, PAR = NULL){
  Pf0 = makePf0(ixH = ixH, tBite = t)
  addEvent2Q(ixH, event_simbitePfSI(t = t, PAR = Pf0))
}

#' \code{event_XX} Family: Simulated Bite Event in EventQ
#'
#' This creates an event for the simulated bite in a human's eventQ and is called by \code{\link{add2Q_simbitePfSI}}.
#' For more information on how this family of functions works, see \url{http://smitdave.github.io/MASH-Development/docHumans.html#add2q_xx}
#'
#' @param ixH index of human
#' @param t time the event will occur
#' @param PAR passed from Pf0
#' @examples
#' getPfParent_SI(ixH)
event_simbitePfSI <- function(t, PAR){
  list(t=t, PAR = PAR, F=simbite_PfSI, tag="simbite_PfSI")
}

#' \code{XX} Family: Simulated Bite
#'
#' This function runs a simulated bite on a human. It is called by \code{\link{event_simbitePfSI}} and simulates \code{\link{probeHost_PfSI}}.
#' For more information on how this family of functions works, see \url{http://smitdave.github.io/MASH-Development/docHumans.html#add2q_xx}
#'
#' @param ixH index of human
#' @param t time the event will occur
#' @param PAR passed from Pf0
#' @examples
#' simbite_PfSI(ixH, t, PAR)
simbite_PfSI <- function(ixH, t, PAR){
  with(PAR,{
    probeHost_PfSI(tBite = t, ixH = ixH, ixS = PfM$ixS, ixM = NULL, Pf = PAR)
  })
}

###################################################################
# Queue up infectious bites
###################################################################

add2Q_startPfSI <- function(ixH, t, pfid){
  addEvent2Q(ixH, event_startPfSI(t, pfid))
}

event_startPfSI <- function(t, pfid){
  list(t=t, PAR=list(pfid=pfid), F=infectHuman_PfSI, tag="infectHuman_PfSI")
}

infectHuman_PfSI <- function(ixH, t, PAR){
  #Infect
  with(PAR,{
    if(HUMANS[[ixH]]$Pathogens$Pf$infected == FALSE & HUMANS[[ixH]]$Pathogens$Pf$chemoprophylaxis == FALSE){
        PfSIHistory(ixH, t, "I")
        HUMANS[[ixH]]$Pathogens$Pf$infected <<- TRUE
        HUMANS[[ixH]]$Pathogens$Pf$t0 <<- t
        HUMANS[[ixH]]$Pathogens$Pf$pfid <<- pfid
        if(rbinom(1,1,FeverPf)){
          add2Q_feverPfSI(ixH, t)
        }
        add2Q_endPfSI(ixH, t, pfid)
      }
  })

}


###################################################################
# End an infection
###################################################################

add2Q_endPfSI <- function(ixH, t, pfid){
  addEvent2Q(ixH, event_endPfSI(t, pfid))
}

event_endPfSI <- function(t,pfid){
  tE = t+ttClearPf()
  list(t=tE, PAR=pfid, F = endPfSI, tag = "endPfSI")
}

endPfSI <- function(ixH, t, PAR=NULL){
  # Clear
  if(HUMANS[[ixH]]$Pathogens$Pf$infected == TRUE){
    PfSIHistory(ixH, t, "S")
    HUMANS[[ixH]]$Pathogens$Pf$infected <<- FALSE
  }
}


###################################################################
# Fever
###################################################################

add2Q_feverPfSI <- function(ixH, t){
  addEvent2Q(ixH, event_feverPfSI(t))
}

event_feverPfSI <- function(t){
  ttF = t + ttFeverPf()
  list(t=ttF, PAR=NULL, F=fever_PfSI, tag = "fever_PfSI")
}

fever_PfSI <- function(ixH, t, PAR=NULL){
  #Fever
  PfSIHistory(ixH, t, "F")
  if(rbinom(1,1,TreatPf)){
    add2Q_treatPfSI(ixH, t)
  }
}


###################################################################
# Treatment
###################################################################

add2Q_treatPfSI <- function(ixH, t){
  addEvent2Q(ixH, event_treatPfSI(t))
}

event_treatPfSI <- function(t){
  ttT = t+ttTreatPf()
  list(t=ttT, PAR=NULL, F=treat_PfSI, tag = "treat_PfSI")
}

treat_PfSI <- function(ixH, t, PAR){
  # Treat
  if(HUMANS[[ixH]]$Pathogens$Pf$infected == TRUE){
    HUMANS[[ixH]]$Pathogens$Pf$infected <<- FALSE
    PfSIHistory(ixH, t, "S")
  }

  HUMANS[[ixH]]$Pathogens$Pf$chemoprophylaxis <<- TRUE
  PfSIHistory(ixH, t, "P")
  # Initiate a period of protection from chemoprophlaxis
  add2Q_endprophylaxisPfSI(ixH, t)
}


###################################################################
# End of Chemoprophylaxis
###################################################################

add2Q_endprophylaxisPfSI <- function(ixH, t){
  addEvent2Q(ixH, event_endprophylaxisPfSI(t))
}

event_endprophylaxisPfSI <- function(t){
  ttS = t + ttSusceptiblePf()
  list(t=ttS, PAR=NULL, F=endprophylaxis_PfSI, tag = "endprophylaxis_PfSI")
}

endprophylaxis_PfSI <- function(ixH, t, PAR){
  # End Prophylaxis
  PfSIHistory(ixH, t, "S")
  HUMANS[[ixH]]$Pathogens$Pf$chemoprophylaxis <<- FALSE
}


###################################################################
# HUMAN PE vaccination functions
###################################################################

add2Q_pevaccinatePfSI <- function(ixH, t){
  addEvent2Q(ixH, event_pevaccinatePfSI(t))
}

event_pevaccinatePfSI <- function(t){
  list(t = t, PAR = NULL, F = pevaccinate_PfSI, tag = "pevaccinate_PfSI")
}

pevaccinate_PfSI <- function(ixH, t, PAR){
  if(rbinom(1,1,PEProtectPf)){
    HUMANS[[ixH]]$Pathogens$Pf$b <<- Pf_b * (1-peBlockPf)
    add2Q_pewanePfSI(ixH, t)
  }
}

add2Q_pewanePfSI <- function(ixH, t){
  ttw = t + ttPEWanePf()
  addEvent2Q(ixH, event_pewanePfSI(ttw))
}

event_pewanePfSI <- function(t){
  list(t=t, PAR=NULL, F=pewane_PfSI, tag="pewane_PfSI")
}

pewane_PfSI <- function(ixH, t, PAR){
  HUMANS[[ixH]]$Pathogens$Pf$b <<- Pf_b
}


###################################################################
# HUMAN GS vaccination functions
###################################################################

add2Q_gsvaccinatePfSI <- function(ixH, t){
  addEvent2Q(ixH, event_gsvaccinatePfSI(t))
}

event_gsvaccinatePfSI <- function(t){
  list(t=t, PAR=NULL, F=gsvaccinate_PfSI, tag = "gsvaccinate_PfSI")
}

gsvaccinate_PfSI <- function(ixH, t, PAR){
  if(rbinom(1,1,GSProtectPf)){
    HUMANS[[ixH]]$Pathogens$Pf$c <<- Pf_c*(1-gsBlockPf)
    add2Q_gswanePfSI(ixH, t)
  }
}

add2Q_gswanePfSI <- function(ixH, t){
  ttw = t + ttGSWane()
  addEvent2Q(ixH, event_gswanePfSI(ttw))
}

event_gswanePfSI <- function(t){
  list(t=t, PAR=NULL, F=gswane_PfSI, tag = "gswane_PfSI")
}

gswane_PfSI <- function(ixH, t, PAR){
 HUMANS[[ixH]]$Pathogens$Pf$c <<- Pf_c
}


###################################################################
# PfSI Diagnostics
###################################################################

rdtTest_PfSI <- function(ixH){
 ifelse(HUMANS[[ixH]]$Pathogens$Pf$infected==TRUE,rbinom(1,1,rdtSensPf), rbinom(1,1,rdtSpecPf))
}

lmTest_PfSI <- function(ixH){
 ifelse(HUMANS[[ixH]]$Pathogens$Pf$infected==TRUE,rbinom(1,1,lmSensPf), rbinom(1,1,lmSpecPf))
}


###################################################################
# PfSI History
###################################################################

PfSIHistory <- function(ixH, t, event){
  if(KeepPfHistory == TRUE){
    HUMANS[[ixH]]$Pathogens$Pf$eventT <<- c(HUMANS[[ixH]]$Pathogens$Pf$eventT,t)
    HUMANS[[ixH]]$Pathogens$Pf$events <<- c(HUMANS[[ixH]]$Pathogens$Pf$events,event)
  }
}

###################################################################
# From bloodstream infection to infect the mosquito
###################################################################

#' Infectious Bite for Human to Vector Transmission
#'
#' This function handles human to vector transmission for the PfSI module and is called from \code{getInfected()} defined in MBITES-HostEncounter.R.
#' If the human has an active infection and there is a successful bloodstream to mosquito transfer of gametocytes a \code{PfM} object is made by \code{makePfM()}
#' and is passed to the mosquito. If global flag \code{\link{PfTransmission_TRACK}} is true it will call \code{\link{trackPfTransmission}} to track human to vector transmission.
#'
#' @param ixH index of human
#' @param t time of bite (local time of mosquito)
#' @param ixS index of feeding site
#' @return a named list
#' * infected: TRUE or FALSE depending on human infection status and transmisson efficiency
#' * PfM: if transmisson successful, the pathogen object passed to the mosquito
#' @md
infectMosquito_PfSI <- function(tBite, ixH, ixS, ixM){
  with(HUMANS[[ixH]]$Pathogens$Pf,{
    if(infected==TRUE & rbinom(1,1,HUMANS[[ixH]]$Pathogens$Pf$c)){
      infObj = makePfM(ixH, tBite, ixS)
      if(PfTransmission_TRACK){
        trackPfTransmission(M2H = FALSE, tBite = tBite, ixH = ixH, ixS = ixS, ixM = ixM, PfM = infObj$PfM)
      }
      return(infObj)
    } else {
      infObj = list(infected = FALSE)
      return(infObj)
    }
  })
}

#' Pathogen Object for Humans in PfSI Module
#'
#' This function generates the Pf object for humans Pathogen slot \code{HUMANS[[ixH]]$Pathogens$Pf <<- pathOBJ_PfSI()}
#' It is called during initialization by \code{PFSI.INIT()}.
#'
#' @param b transmission efficiency: infected mosquito to human
#' @param c transmission efficiency: infected human to mosquito
#' @return a named list
#' * infected: TRUE or FALSE
#' * chemoprophylaxis: is human protected by chemoprophylaxis?
#' @md
pathOBJ_PfSI <- function(b=Pf_b, c=Pf_c){
  list(
    infected = FALSE,
    chemoprophylaxis = FALSE,
    b = b,
    c = c,
    pfid = NULL,
    t0 = NULL,
    eventT = -1,
    events = "init"
  )
}
