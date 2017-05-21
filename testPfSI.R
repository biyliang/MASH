#################################################################
#
#   MASH/MBITES
#   Example Simulation and Analysis of PfSI Module
#   R version
#   Sean Wu
#   April 10, 2017
#
#################################################################

######################################################
# Initialize HUMANS object
######################################################

rm(list=ls())
library(MASH.MBPT)

nH = 1e3 # nH must be divisible by nS for now; fix later.
nS = 100

hhSizes = replicate(n = nS,expr = nH/nS,simplify = TRUE)

hhIx = sapply(X = 1:nS,FUN = function(x,hhSizes){
  ((hhSizes[x]*(x-1))+1) : (hhSizes[x]*x)
},hhSizes = hhSizes)

HUMANS = makeHumans(nH = nH,hhSizes = hhSizes ,hhIx = hhIx)


######################################################
# Initialize PfSI Module
######################################################

PFSI.SETUP(NOISY = FALSE,mnChemoprophylaxisPf = 100)
PFSI.INIT(PfPR = 0)


######################################################
# Simulate PfSI Module
######################################################

# Queue infections
tMax = 8*365
for(i in 1:nH){
  t=0
  print(paste0("queueing for human ixH: ",i))
  while(t < tMax){
    t = t + rexp(1, 1/20)
    add2Q_simbitePfSI(ixH = i,t = t)
  }
}

# Vaccinate
tMax = 5*365
for(i in 1:round(nH/1.5)){
  print(paste0("vaccinating for human ixH: ",i))
  add2Q_pevaccinatePfSI(ixH = i,t = 730)
  add2Q_treatPfSI(ixH = i,t = 730.1)
}

# Simulating
tMax = 15*365
for(i in 1:nH){
  print(paste0("simulating human ixH: ",i))
  liveLife(ixH = i,tPause = tMax)
}

# PLEASE SET TO YOUR OWN DESIRED FOLDER TO HOLD OUTPUT/.. (this must end in /)
out = "/Users/slwu89/Desktop/mash.out/"
writeHumanEvent_PfSI(directory = out,fileName = "humanPfSI.json")


######################################################
# Analysis and Visualization
######################################################

humanHistories = importHumanEvent_PfSI(directory = out)
history = humanHistories

plotPfsiTrajectory(history = history)
