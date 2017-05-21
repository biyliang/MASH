########################################
#  PFSI Setup
#  Sean Wu
#  April 5, 2017
########################################

#' Initialize PfSI Module Parameters
#'
#' Generate a list of parameters and define functions in the global environment for PfSI. All arguments have default values which are listed below before the definition.
#'
#' @param Pf_c 0.15; transmission efficiency: infected human to mosquito
#' @param Pf_b 0.55; transmission efficiency: infected mosquito to human
#' @param DurationPf 200; duration of infection (How many days does the infection last?)
#' @param LatentPf 10; latency (How many days after the infectious bite does the infection start?)
#' @param FeverPf 0.3; probability of fever
#' @param mnFeverPf 10; mean of Gaussian-distributed timing of fever incident (relative to start of the infection)
#' @param vrFeverPf 0.1; standard deviation of Gaussian-distributed timing of fever incident (relative to start of the infection)
#' @param TreatPf 0.5; probability of treatment
#' @param mnTreatPf 1; timing of treatment (relative to start of the fever)
#' @param mnChemoprophylaxisPf 32; Prophylaxis, time to susceptibility (duration of protection due to chemoprophylaxis)
#' @param PEProtectPf 0.99; proportion protected by PE vaccination (probability vaccination successful)
#' @param peBlockPf 1; proportion of infections blocked by PE vaccination
#' @param mnPEPf 270; mean duration of protection from PE vaccination
#' @param vrPEPf 50; standard deviation of protection from PE vaccination
#' @param GSProtectPf 1; proportion protected by GS vaccination (probability vaccination successful)
#' @param gsBlockPf 0.9; proportion of infections blocked by GS vaccination
#' @param mnGSPf 180; mean duration of protection from GS vaccination
#' @param vrGSPf 20; standard deviation of protection from GS vaccination
#' @param rdtSensPf 0.9; RDT sensitivity
#' @param rdtSpecPf 0.1; RDT specificity
#' @param lmSensPf 0.9; Light Microscopy sensitivity
#' @param lmSpecPf 0.1; Light Microscopy specificity
#' @param KeepPfHistory TRUE; record PfSI events in HUMANS[[ixH]]$Pathogens$Pf$eventT?
#' @param NOISY FALSE; print events to console?
#' @return define functions and parameters in global environment
#' @examples
#' PFSI.SETUP()
PFSI.SETUP <- function(

    ########################################
    #  Parameters
    ########################################

    # Transmission efficiency
    Pf_c   = 0.15,
    Pf_b   = 0.55,

    # Duration of Infection
    # (How many days does the infection last?)
    DurationPf = 200,

    # Latency:
    # (How many days after the infectious
    #  bite does the infection start?)
    LatentPf = 10,

    # Probability of Fever
    FeverPf = 0.3,

    # Timing of Fever Incident
    # (relative to start of the infection)
    mnFeverPf = 10,
    vrFeverPf = .1,

    # Probability of Treating
    TreatPf = .5,

    # Timing of Treatment
    # (relative to start of the fever)
    mnTreatPf = 1,

    # Prophylaxis, time to susceptibility
    mnChemoprophylaxisPf = 32,

    # Proportion Protected by PE Vaccination
    PEProtectPf = .99,

    # Proportion of infections Blocked
    peBlockPf = 1,
    mnPEPf = 270,
    vrPEPf = 50,

    # Proportion Protected by GS Vaccination
    GSProtectPf = 1,

    # Proportion of infections Blocked
    gsBlockPf = .9,
    mnGSPf = 180,
    vrGSPf = 20,

    #  RDT Probe and Accuracy
    rdtSensPf  = .9,
    rdtSpecPf  = .1,

    #  Light Microscopy Probe and Sensitivity
    lmSensPf  = .9,
    lmSpecPf  = .1,

    # auxiliary
    KeepPfHistory = TRUE,
    NOISY = FALSE

  ){

  # Define functions
  probeHost <<- probeHost_PfSI
  infectMosquito <<- infectMosquito_PfSI

  # Transmission efficiency
  Pf_c   <<- Pf_c
  Pf_b   <<- Pf_b

  # Duration of Infection
  # (How many days does the infection last?)
  DurationPf <<- DurationPf
  ttClearPf <<- function(){rexp(1, 1/DurationPf)}

  # Latency:
  # (How many days after the infectious
  #  bite does the infection start?)
  LatentPf <<- LatentPf
  ttInfectionPf <<- function(){LatentPf}

  # Probability of Fever
  FeverPf <<- FeverPf

  # Timing of Fever Incident
  # (relative to start of the infection)
  mnFeverPf <<- mnFeverPf
  vrFeverPf <<- vrFeverPf
  ttFeverPf <<- function(){mnFeverPf}

  # Probability of Treating
  TreatPf <<- TreatPf

  # Timing of Treatment
  # (relative to start of the fever)
  mnTreatPf <<- mnTreatPf
  ttTreatPf <<- function(){mnTreatPf}

  # Prophylaxis, time to susceptibility
  mnChemoprophylaxisPf <<- mnChemoprophylaxisPf
  ttSusceptiblePf <<- function(){mnChemoprophylaxisPf}

  # Proportion Protected by PE Vaccination
  PEProtectPf <<- PEProtectPf

  # Proportion of infections Blocked
  peBlockPf <<- peBlockPf
  mnPEPf <<- mnPEPf
  vrPEPf <<- vrPEPf
  ttPEWanePf <<- function(){
    rnorm(1, mnPEPf, vrPEPf)
  }

  # Proportion Protected by GS Vaccination
  GSProtectPf <<- GSProtectPf

  # Proportion of infections Blocked
  gsBlockPf <<- gsBlockPf
  mnGSPf <<- mnGSPf
  vrGSPf <<- vrGSPf
  ttGSWanePf <<- function(){
    rnorm(1, PfmeanGSProtect, PfvarGSprotect)
  }

  #  RDT Probe and Accuracy
  rdtSensPf <<- rdtSensPf
  rdtSpecPf <<- rdtSpecPf
  rdtTest <<- rdtTest_PfSI

  #  Light Microscopy Probe and Sensitivity
  lmSensPf <<- lmSensPf
  lmSpecPf <<- lmSpecPf
  lmTest <<- lmTest_PfSI

  # Pedigree
  getPfParent <<- getPfParent_SI

  # auxiliary
  KeepPfHistory <<- KeepPfHistory
  NOISY <<- NOISY

  # tracking flags
  PfPedigree_TRACK <<- FALSE
  PfTransmission_TRACK <<- FALSE

  # placeholder objects for simulated biting
  Pf0 <<- list(spz = TRUE)
  Pf0$PfM[[1]] <<- list(tm=0, ixS=0, ixH=0, damID=0, sireID=0, pfid=0)

}

#' Initialize PfSI Human Infections
#'
#' Initialize PfSI module infections in the human population (requires \code{HUMANS} to be defined in the global environment)
#' \code{PFSI.SETUP()} should be called prior to initializing human infections.
#'
#' @param PfPR initial parasite prevalence in the human population
#' @return modify \code{HUMANS} object in global environment
#' @examples
#' PFSI.INIT()
PFSI.INIT <- function(PfPR){
  if(exists(x = "LANDSCAPE",where = .GlobalEnv)){
    if(length(HUMANS)!=LANDSCAPE$nH){
      stop("HUMANS object not equal to LANDSCAPE$nH")
    }
  }
  if(!exists(x = "tStart",where = .GlobalEnv)){
    print(paste0("tStart not defined globally; setting initial infection times to 0"))
    tStart = 0
  }
  PfID <<- 1
  for(ixH in 1:length(HUMANS)){
    HUMANS[[ixH]]$Pathogens$Pf <<- pathOBJ_PfSI()
    if(runif(1) < PfPR){
      infectHuman_PfSI(ixH = ixH,t = tStart,PAR = list(pfid = PfID))
      PfID <<- PfID + 1
    } else {
      PfSIHistory(ixH = ixH,t = tStart,event = "S")
    }
  }
}
