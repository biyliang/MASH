#################################################################
#
#   MASH/MBITES
#   Example Simulation and Analysis MASH-MBPT
#   R version
#   Sean Wu
#   April 10, 2017
#
#################################################################

######################################################
# Initialize Ross-MacDonald parameters and EL4P
######################################################

rm(list=ls())
library(MASH.MBPT)

# global parameters
P = MBITES.PAR(maxBatch = 35)

# set up LANDSCAPE
LANDSCAPE = makeLandscape(nF = 25,nA = 25,nS = 5,nM = 5,aquaMod = "el4p",pointGen = pointsPoisson,hhSize = 10,hhMin = 2,aquaSD = 0.05)
staticLandscapePlot()
MvOb = exactAll(LANDSCAPE)

# test MBITES-Basic
basicN = 5e3
basicOut = MBITES.basic(N = basicN,P = P,parallel = TRUE,move = TRUE)

cohortBionomics = cohortBionomics(mosyPop = basicOut,ix = 1:basicN)
aquaEq = aquaIx_equilibrium(mosyPop = basicOut)

eqM = 50

Aq.PAR = makePAR_EL4P(nA = LANDSCAPE$nA,nH = LANDSCAPE$nH,aquaEq = aquaEq,M = floor(eqM*aquaEq),R0 = 4,par = P,summary = cohortBionomics$summary)
# Aq.Equilibrium = setupAquaPop_EL4PsamplePoints(PAR = Aq.PAR,tol = 100,plot = TRUE)
Aq.Equilibrium = setupAquaPop_EL4Pexact(PAR = Aq.PAR,tol = 100,plot = TRUE)

# set LANDSCAPE aquatic populations to equilibrium values
for(ix in 1:Aq.PAR$nA){
  # probably need to have the eggs go into the EggQs as well
  LANDSCAPE$aquaSites[[ix]]$EL4P = Aq.Equilibrium$EL4P_pops[[ix]]
  LANDSCAPE$aquaSites[[ix]]$alpha = Aq.Equilibrium$PAR$alpha[[ix]]
  LANDSCAPE$aquaSites[[ix]]$psi = Aq.Equilibrium$PAR$psiHat[[ix]]
  LANDSCAPE$aquaSites[[ix]]$p = Aq.Equilibrium$PAR$p
}

######################################################
# Initialize HUMANS and mosquito objects
######################################################

# generate humans object
HUMANS = makeHumans(nH = LANDSCAPE$nH,hhSizes = LANDSCAPE$hhSizes,hhIx = LANDSCAPE$hhIx)

# initialize activity space
init_riskList() #initialize riskList for all sites
init_myTimeAtRisk() #initialize myTimeAtRisk parameters for all humans
activitySpace() #run daily activity space

# generate mosquito populations
tStart = 1
MPopF = makeMosquitoCohort(N = eqM,female = TRUE,tm = tStart-1, state = "M", EIP = 1, mature = TRUE, offset = 6e4)
MPopM = makeMosquitoCohort(N = eqM,female = FALSE,tm = tStart-1, state = "M", offset = 6e4)

# initialize PfSI module
PFSI.SETUP()
PFSI.INIT(PfPR = 0.5)


##########################################
# Run MASH
##########################################

# PLEASE SET TO YOUR OWN DESIRED FOLDER TO HOLD OUTPUT/.. (this must end in /)
out = "/Users/slwu89/Desktop/mash.out/"

clearOutput(directory = out) # THIS FUNCTION ERASES ALL FILES IN OUTPUT/ FOLDER; USE WITH CAUTION

# vector population tracking
el4pCon = trackEL4P_init(directory = out,fileName = "el4p.csv")
adultCon = trackAdults_init(directory = out,fileName = "adults.csv")
EggQCon = trackEggQ_init(directory = out,fileName = "egg.csv")

# pathogen transmission tracking
PfPedigreeCon = PfPedigree_init(directory = out,fileName = "pfpedigree.csv") # this MUST be called PfPedigreeCon
PfTransmissionCon = PfTransmission_init(directory = out,fileName = "pftransmission.csv") # this MUST be called PfTransmissionCon

# test MASH
tEnd = tStart+(361*2)
for(tMax in tStart:tEnd){

  print(paste0("tMax: ",tMax))

  # human activity space simulation
  activitySpace()

  # "EL4P" Aquatic Ecology
  oneDay_EL4P(tNow = tMax)
  addCohort(tNow = tMax)

  # male MBITES
  MBITESmale(P)

  # female MBITES
  MBITES(P)

  # human events
  runHumanEventQ(tPause = tMax)

  # track all output
  trackEL4P(con = el4pCon)
  trackAdults(con = adultCon)
  # log mosquito data and clear out vector every 10 days
  if((tMax > tStart) & (tMax %% 10 == 0)){
    print(paste0("logging mosquito bionomics and histories"))
    trackBionomics(directory = out,fileName = paste0("bionomics",tMax,".json"))
    trackHistory(directory = out,fileName = paste0("historyF",tMax,".json"))
    trackHistoryM(directory = out,fileName = paste0("historyM",tMax,".json"))
    resetMosyPop(female = TRUE)
    resetMosyPop(female = FALSE)
  }
  # log human event data when simulation completes
  if(tMax == tEnd){
    writeHumanEvent_PfSI(directory = out,fileName = "humanPfSI.json")
  }

}


# close all output connections
close(el4pCon)
close(adultCon)
close(EggQCon)

close(PfPedigreeCon)
close(PfTransmissionCon)


##########################################
# Analyze MASH
##########################################

# import .json data
bionomics = importBionomics(directory = out) # import mosquito bionomics
history = importHistory(directory = out) # import mosquito event histories
humanHistory = importHumanEvent_PfSI(directory = out) # import human event histories

# import .csv data
adults = importAdults(directory = out,fileName = "adults.csv") # import mosquito adult densities
el4p = importEL4P(directory = out,fileName = "el4p.csv") # import mosquito aquatic stage densities
egg = importEggQ(directory = out,fileName = "egg.csv") # import mosquito egg laying densities

PfTransmission = importPfTransmission(directory = out,fileName = "pftransmission.csv") # import pf transmission histories
PfPedigree = importPfPedigree(directory = out,fileName = "pfpedigree.csv") # import pf pedigree histories (currently placeholder)

# analyze aquatic dynamics
plotEL4P(el4p = el4p,egg = egg)

# analyze adult dynamics
plotAdults(adults)

# analyze mosquito cohorts
cohortIx = getCohortIndices.history(history = history,sites = FALSE)
cohortT = getStateTraj.history(ix = cohortIx[[1]],history = history,female = TRUE)
plotStateTrajectory(stateT = cohortT,logX = F,normY = T, lwd = 1.5)
cohortTrajectory.history(ix = cohortIx[[1]],history = history,cex = 2)

cohortBionomics.history(bionomics = bionomics,ix = cohortIx[[1]])

# analyze human infection
plotPfsiTrajectory(history = humanHistory)
