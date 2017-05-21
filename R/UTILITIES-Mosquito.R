#################################################################
#
#   MASH/MBITES
#   Mosquito Profile
#   R version
#   Sean Wu
#   January 24, 2017
#
#################################################################


##########################################
# Auxiliary Functions
##########################################

# extractPointTraj: extract xy trajectory from a mosquito history object
extractPointTraj <- function(history){

  xyH = with(history,{

    pointHist = switch(pSetH[1],
                       f = LANDSCAPE$feedSites[[ixH[1]]]$siteXY,
                       l = LANDSCAPE$aquaSites[[ixH[1]]]$siteXY,
                       s = LANDSCAPE$sugarSites[[ixH[1]]]$siteXY,
                       m = LANDSCAPE$swarmSites[[ixH[1]]]$siteXY
    )

    for(i in 2:length(ixH)){

      if((ixH[i] != ixH[(i-1)]) | (pSetH[i] != pSetH[(i-1)])){

        pointI = switch(pSetH[i],
                        f = LANDSCAPE$feedSites[[ixH[i]]]$siteXY,
                        l = LANDSCAPE$aquaSites[[ixH[i]]]$siteXY,
                        s = LANDSCAPE$sugarSites[[ixH[i]]]$siteXY,
                        m = LANDSCAPE$swarmSites[[ixH[i]]]$siteXY
        )
        pointHist = rbind(pointHist,pointI)
      }

    }

    return(unname(pointHist))
  })

  return(xyH)
}

# extractStateTraj: extract state trajectory from a mosquito history object
extractStateTraj <- function(history){
  list(state=history$stateH,time=history$timeH)
}


##########################################
# Mosquito Landscape Trajectories
##########################################

# oneTrajectory: plot a single mosquito trajectory
# history: history object
# traceCol: hue used for the trajectory
# bgcol: background color of plot the trajectory will be overlaid onto
oneTrajectory <- function(history, traceCol = "#FF80FFFF", bgCol = "grey20", offset = 0.05, cex = 0.75){

  feedXY = t(sapply(LANDSCAPE$feedSites,function(x){x$siteXY}))
  aquaXY = t(sapply(LANDSCAPE$aquaSites,function(x){x$siteXY}))
  sugarXY = t(sapply(LANDSCAPE$sugarSites,function(x){x$siteXY}))
  mateXY = t(sapply(LANDSCAPE$swarmSites,function(x){x$siteXY}))

  #set up plotting options
  if(!exists(x = "defaultPar",where = .GlobalEnv)){
    # if does not exist, save default graphical parameters
    .GlobalEnv$defaultPar = par()
  }
  par(bg=bgCol,mar=c(0,0,0,0),mgp=c(0,0,0))

  #colors
  setup_col = ggCol(n=4)

  #set up empty grid
  plot(1,type="n",axes=F,frame.plot=F,ann=F,xlim=c(LANDSCAPE$xLim[1]-offset,LANDSCAPE$xLim[2]+offset),
       ylim=c(LANDSCAPE$yLim[1]-offset,LANDSCAPE$yLim[2]+offset))
  # plot(seq(from=(LANDSCAPE$xLim[1]-offset),to=(LANDSCAPE$xLim[2]+offset),by=offset),
  #      seq(from=(LANDSCAPE$yLim[1]-offset),to=(LANDSCAPE$yLim[2]+offset),by=offset)
  #      ,type="n",axes=FALSE,frame.plot=FALSE,ann=FALSE)

  points(x = feedXY, pch = 15, col = setup_col[1], cex = cex) #draw feeding sites
  points(x = aquaXY, pch = 16, col = setup_col[2], cex = cex) #draw aquatic habitats
  points(x = sugarXY, pch = 17, col = setup_col[3], cex = cex) #draw sugar sources
  points(x = mateXY, pch = 18, col = setup_col[4], cex = cex) #draw mating sites

  legend(x = "topleft", legend = c("Feeding Site","Aquatic Habitat","Sugar Source","Mating Site"),pch = 15:18,
         col = setup_col, bty = "n", text.col = "grey80")


  xyH = with(history,{

    pointHist = switch(pSetH[1],
                       f = LANDSCAPE$feedSites[[ixH[1]]]$siteXY,
                       l = LANDSCAPE$aquaSites[[ixH[1]]]$siteXY,
                       s = LANDSCAPE$sugarSites[[ixH[1]]]$siteXY,
                       m = LANDSCAPE$swarmSites[[ixH[1]]]$siteXY
    )

    for(i in 2:length(ixH)){

      if((ixH[i] != ixH[(i-1)]) | (pSetH[i] != pSetH[(i-1)])){

        pointI = switch(pSetH[i],
                        f = LANDSCAPE$feedSites[[ixH[i]]]$siteXY,
                        l = LANDSCAPE$aquaSites[[ixH[i]]]$siteXY,
                        s = LANDSCAPE$sugarSites[[ixH[i]]]$siteXY,
                        m = LANDSCAPE$swarmSites[[ixH[i]]]$siteXY
        )
        pointHist = rbind(pointHist,pointI)
      }

    }

    return(unname(pointHist))
  })

  colVec = colorRampPalette(colors=c(bgCol,traceCol),alpha=0.8)(nrow(xyH)) # 1st element is brightest, fades to bgCol

  for(i in 1:(nrow(xyH)-1)){
    segments(x0 = xyH[i,1], y0 = xyH[i,2], x1 = xyH[(i+1),1], y1 = xyH[(i+1),2], col = colVec[i+1])
  }

  par(bg=.GlobalEnv$defaultPar$bg,mar=.GlobalEnv$defaultPar$mar,mgp=.GlobalEnv$defaultPar$mgp)

}

# cohortTrajectory: plot all mosquito trajectories for a cohort
# ix: vector of mosquito ID belonging to a single cohort (use getCohortIndices() to extract this)
# MPop: mosquito population object
# traceCol: hue used for the trajectory
# bgcol: background color of plot the trajectory will be overlaid onto
# jitterSD: if != 0 add rnorm(0,jitterSD) amount of jitter to trajectories
cohortTrajectory <- function(ix, MPop, traceCol = "#FF80FFFF", bgCol = "grey20", offset = 0.05, cex = 0.75, alpha = 0.5){

  feedXY = t(sapply(LANDSCAPE$feedSites,function(x){x$siteXY}))
  aquaXY = t(sapply(LANDSCAPE$aquaSites,function(x){x$siteXY}))
  sugarXY = t(sapply(LANDSCAPE$sugarSites,function(x){x$siteXY}))
  mateXY = t(sapply(LANDSCAPE$swarmSites,function(x){x$siteXY}))

  #set up plotting options
  .GlobalEnv$defaultPar = par()
  par(bg=bgCol,mar=c(0,0,0,0),mgp=c(0,0,0))

  #colors
  setup_col = ggCol(n=4)

  #set up empty grid
  plot(1,type="n",axes=F,frame.plot=F,ann=F,xlim=c(LANDSCAPE$xLim[1]-offset,LANDSCAPE$xLim[2]+offset),
       ylim=c(LANDSCAPE$yLim[1]-offset,LANDSCAPE$yLim[2]+offset))
  # plot(seq(from=(LANDSCAPE$xLim[1]-offset),to=(LANDSCAPE$xLim[2]+offset),by=offset),
  #      seq(from=(LANDSCAPE$yLim[1]-offset),to=(LANDSCAPE$yLim[2]+offset),by=offset)
  #      ,type="n",axes=FALSE,frame.plot=FALSE,ann=FALSE)

  points(x = feedXY, pch = 15, col = setup_col[1], cex = cex) #draw feeding sites
  points(x = aquaXY, pch = 16, col = setup_col[2], cex = cex) #draw aquatic habitats
  points(x = sugarXY, pch = 17, col = setup_col[3], cex = cex) #draw sugar sources
  points(x = mateXY, pch = 18, col = setup_col[4], cex = cex) #draw mating sites

  legend(x = "topleft", legend = c("Feeding Site","Aquatic Habitat","Sugar Source","Mating Site"),pch = 15:18,
         col = setup_col, bty = "n", text.col = "grey80")

  # extract all xy trajectories from that cohort
  for(i in 1:length(ix)){
    xyH = extractPointTraj(MPop$mosy[[ix[i]]]$history)
    colVec = colorRampPalette(colors=c(bgCol,traceCol),alpha=alpha)(nrow(xyH)) # 1st element is brightest, fades to bgCol
    for(j in 1:(nrow(xyH)-1)){
      segments(x0 = xyH[j,1], y0 = xyH[j,2], x1 = xyH[(j+1),1], y1 = xyH[(j+1),2], col = colVec[j+1])
    }
  }

  suppressWarnings(par(.GlobalEnv$defaultPar))
}

# cohortTrajectory: plot all mosquito trajectories for a cohort from .json history
# ix: vector of mosquito ID belonging to a single cohort (use getCohortIndices() to extract this)
# MPop: mosquito population object
# traceCol: hue used for the trajectory
# bgcol: background color of plot the trajectory will be overlaid onto
# jitterSD: if != 0 add rnorm(0,jitterSD) amount of jitter to trajectories
# NOTE: requires the LANDSCAPE the cohort was run on to be defined in global environment
cohortTrajectory.history <- function(ix, history, traceCol = "#FF80FFFF", bgCol = "grey20", offset = 0.05, cex = 0.75, alpha = 0.5){

  feedXY = t(sapply(LANDSCAPE$feedSites,function(x){x$siteXY}))
  aquaXY = t(sapply(LANDSCAPE$aquaSites,function(x){x$siteXY}))
  sugarXY = t(sapply(LANDSCAPE$sugarSites,function(x){x$siteXY}))
  mateXY = t(sapply(LANDSCAPE$swarmSites,function(x){x$siteXY}))

  #set up plotting options
  .GlobalEnv$defaultPar = par()
  par(bg=bgCol,mar=c(0,0,0,0),mgp=c(0,0,0))

  #colors
  setup_col = ggCol(n=4)

  #set up empty grid
  plot(1,type="n",axes=F,frame.plot=F,ann=F,xlim=c(LANDSCAPE$xLim[1]-offset,LANDSCAPE$xLim[2]+offset),
       ylim=c(LANDSCAPE$yLim[1]-offset,LANDSCAPE$yLim[2]+offset))
  # plot(seq(from=(LANDSCAPE$xLim[1]-offset),to=(LANDSCAPE$xLim[2]+offset),by=offset),
  #      seq(from=(LANDSCAPE$yLim[1]-offset),to=(LANDSCAPE$yLim[2]+offset),by=offset)
  #      ,type="n",axes=FALSE,frame.plot=FALSE,ann=FALSE)

  points(x = feedXY, pch = 15, col = setup_col[1], cex = cex) #draw feeding sites
  points(x = aquaXY, pch = 16, col = setup_col[2], cex = cex) #draw aquatic habitats
  points(x = sugarXY, pch = 17, col = setup_col[3], cex = cex) #draw sugar sources
  points(x = mateXY, pch = 18, col = setup_col[4], cex = cex) #draw mating sites

  legend(x = "topleft", legend = c("Feeding Site","Aquatic Habitat","Sugar Source","Mating Site"),pch = 15:18,
         col = setup_col, bty = "n", text.col = "grey80")

  # extract all xy trajectories from that cohort
  for(i in 1:length(ix)){
    xyH = extractPointTraj(history[[ix[i]]])
    colVec = colorRampPalette(colors=c(bgCol,traceCol),alpha=alpha)(nrow(xyH)) # 1st element is brightest, fades to bgCol
    for(j in 1:(nrow(xyH)-1)){
      segments(x0 = xyH[j,1], y0 = xyH[j,2], x1 = xyH[(j+1),1], y1 = xyH[(j+1),2], col = colVec[j+1])
    }
  }

  suppressWarnings(par(.GlobalEnv$defaultPar))
}


##########################################
# Mosquito State Space Trajectories
##########################################

# cohortStateTrajectory: plot a mosquito cohort through state space
# stateT: state trajectory of a cohort (use getStateTraj() to extract this)
# logX: log-transform time axis? (note: 1 will be added to prevent numerical inconsistencies)
# normY: plot normalized proportions or counts?
# ...: additional named parameters passed to lines()
plotStateTrajectory <- function(stateT,logX = TRUE, normY = TRUE, ...){

  times = sapply(stateT,function(x){x$time})
  timesLabel = times
  maxY = max(unlist(stateT[[1]][-1]))

  if(logX){
    if(0 %in% times){
      times = log(times+1)
    } else {
      times = log(times)
    }
    xTitle = "Time (log)"
  } else {
    xTitle = "Time"
  }

  if(normY){
    yLimU = 1
    yTitle = "Proportion"
  } else {
    yLimU = maxY
    yTitle = "Count"
  }

  stateSpace = names(stateT[[1]])[-1]
  stateCol = ggCol(n = length(stateSpace))

  # set up empty plot with correct parameters
  # plot(1, type = "n", xaxt = "n", xlim = c(0,max(times)),ylim = c(0,yLimU), xlab = xTitle, ylab = yTitle, main = "Cohort State Space Trajectory")
  # axis(side = 1, at = times, labels = c(0,rep(x = NA,times = length(times)-2),round(max(timesLabel),digits = 3)))
  plot(1, type = "n", xaxt = "n", xlim = c(min(times),max(times)),ylim = c(0,yLimU), xlab = xTitle, ylab = yTitle, main = "Cohort State Space Trajectory")
  axis(side = 1, at = times, labels = c(min(times),rep(x = NA,times = length(times)-2),round(max(timesLabel),digits = 3)))
  grid()

  # plot state trajectories
  for(st in 1:length(stateSpace)){
    stateY = sapply(stateT,function(x){x[[stateSpace[st]]]})
    if(normY){
      stateY = stateY / maxY
    }
    lines(x = times,y = stateY,col = stateCol[st], ...)
  }

  legend("topright",legend = stateSpace,col = stateCol,bty = "n",pch = 16)

}

# getStateTraj.history: get state trajectory for a single cohort from .json history
# ix: vector of mosquito ID belonging to a single cohort (use getCohortIndices() to extract this)
# MPop: mosquito population object
getStateTraj.history <- function(ix, history, female = TRUE){

  #extract the trajectory of a cohort
  cohortTraj = lapply(history[ix],function(x){extractStateTraj(x)})
  cohortN = length(cohortTraj)

  # create time bin for each discrete event
  eventTimes = unname(sort(unlist(sapply(cohortTraj,function(x){x$time}))))
  timeBins = unname(c(min(eventTimes),Filter(f = function(y){y > min(eventTimes)},x = eventTimes)))

  # create empty time series (each element is slice of occupancy vector at that point in time)
  timeSeries = lapply(X = timeBins, oneTime, female = female) # bins

  initState = unique(sapply(cohortTraj,function(x){x$state[1]})) # initial state
  if(length(initState)>1){ # sanity check
    stop("more than one initial state")
  }

  # set initial occupancy vector
  switch(initState,
         "F" = for(i in 1:length(timeSeries)){timeSeries[[i]]$F <- cohortN},
         "B" = for(i in 1:length(timeSeries)){timeSeries[[i]]$B <- cohortN},
         "R" = for(i in 1:length(timeSeries)){timeSeries[[i]]$R <- cohortN},
         "L" = for(i in 1:length(timeSeries)){timeSeries[[i]]$L <- cohortN},
         "O" = for(i in 1:length(timeSeries)){timeSeries[[i]]$O <- cohortN},
         "M" = for(i in 1:length(timeSeries)){timeSeries[[i]]$M <- cohortN},
         "S" = for(i in 1:length(timeSeries)){timeSeries[[i]]$S <- cohortN},
         "D" = stop("initial state should not be D")
  )

  for(ixM in 1:cohortN){ # iterate over mosquitoes

    time = cohortTraj[[ixM]]$time
    state = cohortTraj[[ixM]]$state
    for(ixT in 2:length(time)){
      tIter = which(timeBins >= time[ixT]) # times over which to propagate current state change
      for(i in tIter){
        timeSeries[[i]][[state[ixT-1]]] = timeSeries[[i]][[state[ixT-1]]] - 1 # mosy ixM no longer in state[ixT-1] in all times ixT
        timeSeries[[i]][[state[ixT]]] = timeSeries[[i]][[state[ixT]]] + 1 # mosy ixM in state[ixT] is propagated forward
      }
    }

  }

  return(timeSeries)
}


# getStateTraj: get state trajectory for a single cohort
# ix: vector of mosquito ID belonging to a single cohort (use getCohortIndices() to extract this)
# MPop: mosquito population object
getStateTraj <- function(ix, mosyPop, female = TRUE){

  #extract the trajectory of a cohort
  cohortTraj = lapply(mosyPop[ix],function(x){extractStateTraj(x$history)})
  cohortN = length(cohortTraj)

  # create time bin for each discrete event
  eventTimes = sort(unlist(sapply(cohortTraj,function(x){x$time})))
  timeBins = c(min(eventTimes),Filter(f = function(y){y > min(eventTimes)},x = eventTimes))

  timeSeries = lapply(X = timeBins, oneTime, female = female) # bins

  initState = unique(sapply(cohortTraj,function(x){x$state[1]})) # initial state

  if(length(initState)>1){
    stop("more than one initial state")
  }

  # set initial occupancy vector
  switch(initState,
         "F" = for(i in 1:length(timeSeries)){timeSeries[[i]]$F <- cohortN},
         "B" = for(i in 1:length(timeSeries)){timeSeries[[i]]$B <- cohortN},
         "R" = for(i in 1:length(timeSeries)){timeSeries[[i]]$R <- cohortN},
         "L" = for(i in 1:length(timeSeries)){timeSeries[[i]]$L <- cohortN},
         "O" = for(i in 1:length(timeSeries)){timeSeries[[i]]$O <- cohortN},
         "M" = for(i in 1:length(timeSeries)){timeSeries[[i]]$M <- cohortN},
         "S" = for(i in 1:length(timeSeries)){timeSeries[[i]]$S <- cohortN},
         "D" = stop("initial state should not be D")
  )

  for(ixM in 1:cohortN){ # iterate over mosquitoes

    time = cohortTraj[[ixM]]$time
    state = cohortTraj[[ixM]]$state
    for(ixT in 2:length(time)){
      tIter = which(timeBins >= time[ixT]) # times over which to propagate current state change
      for(i in tIter){
        timeSeries[[i]][[state[ixT-1]]] = timeSeries[[i]][[state[ixT-1]]] - 1 # mosy ixM no longer in state[ixT-1] in all times ixT
        timeSeries[[i]][[state[ixT]]] = timeSeries[[i]][[state[ixT]]] + 1 # mosy ixM in state[ixT] is propagated forward
      }
    }

  }

  return(timeSeries)
}

# single time slice of a cohort's occupancy vector
oneTime <- function(time, female){
  if(female){
    list(
      time = time,
      F = 0,
      B = 0,
      R = 0,
      L = 0,
      O = 0,
      M = 0,
      S = 0,
      D = 0
    )
  } else {
    list(
      time = time,
      R = 0,
      M = 0,
      S = 0,
      D = 0
    )
  }
}


####################################################################################
#
# Retrieve mosquito indices by date of emergence
# or by date and site of emergence for a cohort
#
####################################################################################

# getCohortIndices: return named list where each element is the cohort of mosquitoes with a shared day of emergence
# sites: partition by sites?
getCohortIndices <- function(MPop, sites = FALSE){

  emergeDays = unlist(lapply(MPop$mosy,function(x){x$bDay}))
  emergeDaysUnique = unique(emergeDays) # all unique days of emergence
  cohortDay = lapply(emergeDaysUnique,function(x){which(emergeDays == x)}) # list of ids for each day of emergence

  if(sites){

    sites = sapply(LANDSCAPE$aquaSites,function(x){x$ix}) # vector of site ID
    iter = expand.grid(ixS=sites,ixD=emergeDaysUnique) # create iterator object

    cohortDaySite = parallel::mcmapply(FUN = function(ixS,ixD){
      filtM = Filter(f = function(xx){xx$bDay == ixD & xx$history$ixH[1] == ixS},x = MPop$mosy) # filter MPop to those born on ixD and at site ixS
      return(sapply(X = filtM,function(x){x$id}))
    },ixS=iter$ixS,ixD=iter$ixD)

    names(cohortDaySite) = paste0(iter$ixS,":",iter$ixD) # ixS : ixD
    return(cohortDaySite)

  } else {

    names(cohortDay) = emergeDaysUnique
    return(cohortDay)

  }
}

# getCohortIndices.history: return named list where each element is the cohort of mosquitoes with a shared day of emergence from imported .json histories
# sites: partition by sites?
getCohortIndices.history <- function(history, sites = FALSE){

  emergeDays = unlist(lapply(history,function(x){x$timeH[1]}))
  emergeDaysUnique = unique(emergeDays) # all unique days of emergence

  if(sites){

    print("warning: running 'sites' currently requires having the same LANDSCAPE used in the simulation defined in the global environment")
    sites = sapply(LANDSCAPE$aquaSites,function(x){x$ix}) # vector of site ID
    iter = expand.grid(ixS=sites,ixD=emergeDaysUnique) # create iterator object

    emergeSites = unlist(lapply(history,function(x){x$ixH[1]}))
    ixM = 1:length(history)
    names(ixM) = names(history)

    # pull out mosquito ix in history and tagged id from M-BITES by date and site of emergence
    cohortDaySite = parallel::mcmapply(FUN = function(ixS,ixD,emergeSites,emergeDays,ixM){

      ixC = which(emergeDays == ixD & emergeSites == ixS)
      names(ixC) = names(ixM[ixC])
      return(ixC)

    },ixS=iter$ixS,ixD=iter$ixD,MoreArgs = list(emergeSites = emergeSites, emergeDays = emergeDays, ixM = ixM))

    names(cohortDaySite) = paste0(iter$ixS,":",iter$ixD) # ixS : ixD (site index : day index)
    return(cohortDaySite)

  } else {

    cohortDay = lapply(emergeDaysUnique,function(x){which(emergeDays == x)}) # list of ids for each day of emergence
    names(cohortDay) = emergeDaysUnique
    return(cohortDay)

  }

}


##########################################
# Mosquito Cohort Summary Bionomics
#
# Derive basic bionomics from a cohort
# for RM calculations
#
# Correspondence with original DHM-Basic:
# S = stateH
# T = timeH
# f = bmInt
# a = bmIntH
# Hfeed.N = feedHumanH
# feed.N = feedAllH
# libfa.span = lifespan
# m.eg = mBatch
# t.eg = tBatch
#
##########################################

# cohortBionomics: calculate bionomics for a single cohort
# ix: index of a cohort
cohortBionomics <- function(mosyPop, ix){
  feedAllC = feedHumanC = lifespanC = mBatchC = tBatchC = rep(NA,length=length(ix))
  bmIntC = bmIntHC = rep(NA,length=length(ix)*3)

  i = j = k = 1
  for(ixM in ix){
    # calculate summary bionomics for mosy ixM
    bmIntLen = length(mosyPop[[ixM]]$bionomics$bmInt)
    if(bmIntLen>0){
      bmIntC[j:(j+bmIntLen-1)] = mosyPop[[ixM]]$bionomics$bmInt
      j = j + bmIntLen
    }
    bmIntHLen = length(mosyPop[[ixM]]$bionomics$bmIntH)
    if(bmIntHLen>0){
      bmIntHC[k:(k+bmIntHLen-1)] = mosyPop[[ixM]]$bionomics$bmIntH
      k = k + bmIntHLen
    }
    feedAllC[i] = mosyPop[[ixM]]$bionomics$feedAllH
    feedHumanC[i] = mosyPop[[ixM]]$bionomics$feedHumanH
    lifespanC[i] = mosyPop[[ixM]]$bionomics$lifespan
    mBatchC[i] = mosyPop[[ixM]]$bionomics$mBatch
    tBatchC[i] = mosyPop[[ixM]]$bionomics$tBatch
    i = i + 1
  }

  bmIntC = Filter(Negate(is.na),x = bmIntC)
  bmIntHC = Filter(Negate(is.na),x = bmIntHC)

  # return summary statistics for cohort
  return(list(
    data = list(feedAllC=feedAllC, feedHumanC=feedHumanC, lifespanC=lifespanC, mBatchC=mBatchC, tBatchC=tBatchC),
    summary = list(feedAllC=mean(feedAllC), feedHumanC=mean(feedHumanC), lifespanC=mean(lifespanC), mBatchC=mean(mBatchC), tBatchC=mean(tBatchC))
  ))
}

cohortBionomics.history <- function(bionomics, ix){
  feedAllC = feedHumanC = lifespanC = mBatchC = tBatchC = rep(NA,length=length(ix))
  bmIntC = bmIntHC = rep(NA,length=length(ix)*3)

  i = j = k = 1
  for(ixM in ix){
    # calculate summary bionomics for mosy ixM
    bmIntLen = length(bionomics[[ixM]]$bmInt)
    if(bmIntLen>0){
      bmIntC[j:(j+bmIntLen-1)] = bionomics[[ixM]]$bmInt
      j = j + bmIntLen
    }
    bmIntHLen = length(bionomics[[ixM]]$bmIntH)
    if(bmIntHLen>0){
      bmIntHC[k:(k+bmIntHLen-1)] = bionomics[[ixM]]$bmIntH
      k = k + bmIntHLen
    }
    feedAllC[i] = bionomics[[ixM]]$feedAllH
    feedHumanC[i] = bionomics[[ixM]]$feedHumanH
    lifespanC[i] = bionomics[[ixM]]$lifespan
    mBatchC[i] = bionomics[[ixM]]$mBatch
    tBatchC[i] = bionomics[[ixM]]$tBatch
    i = i + 1
  }

  bmIntC = Filter(Negate(is.na),x = bmIntC)
  bmIntHC = Filter(Negate(is.na),x = bmIntHC)

  # return summary statistics for cohort
  return(list(
    data = list(feedAllC=feedAllC, feedHumanC=feedHumanC, lifespanC=lifespanC, mBatchC=mBatchC, tBatchC=tBatchC),
    summary = list(feedAllC=mean(feedAllC), feedHumanC=mean(feedHumanC), lifespanC=mean(lifespanC), mBatchC=mean(mBatchC), tBatchC=mean(tBatchC))
  ))
}

# cohortBionomics.basic: simplification of cohortBionomics for MBITES-Basic
# MPop: output of MBITES.basic
cohortBionomics.basic <- function(MPop){

  ix = length(MPop)

  feedAllC = feedHumanC = lifespanC = mBatchC = tBatchC = rep(NA,length=length(ix))
  bmIntC = bmIntHC = rep(NA,length=length(ix)*3)

  for(ixM in ix){
    # calculate summary bionomics for mosy ixM
    bmIntLen = length(MPop[[ixM]]$bionomics$bmInt)
    if(bmIntLen>0){
      bmIntC[j:(j+bmIntLen-1)] = MPop$mosy[[ixM]]$bionomics$bmInt
      j = j + bmIntLen
    }
    bmIntHLen = length(MPop[[ixM]]$bionomics$bmIntH)
    if(bmIntHLen>0){
      bmIntHC[k:(k+bmIntHLen-1)] = MPop$mosy[[ixM]]$bionomics$bmIntH
      k = k + bmIntHLen
    }
    feedAllC[i] = MPop[[ixM]]$bionomics$feedAllH
    feedHumanC[i] = MPop[[ixM]]$bionomics$feedHumanH
    lifespanC[i] = MPop[[ixM]]$bionomics$lifespan
    mBatchC[i] = MPop[[ixM]]$bionomics$mBatch
    tBatchC[i] = MPop[[ixM]]$bionomics$tBatch
    i = i + 1
  }

  bmIntC = Filter(Negate(is.na),x = bmIntC)
  bmIntHC = Filter(Negate(is.na),x = bmIntHC)

  # return summary statistics for cohort
  return(list(
    data = list(feedAllC=feedAllC, feedHumanC=feedHumanC, lifespanC=lifespanC, mBatchC=mBatchC, tBatchC=tBatchC),
    summary = list(feedAllC=mean(feedAllC), feedHumanC=mean(feedHumanC), lifespanC=mean(lifespanC), mBatchC=mean(mBatchC), tBatchC=mean(tBatchC))
  ))
}


##############################################################################################################################
# Functions to calculate empirical equilibrium distribution of time at aquatic habitats
##############################################################################################################################

# aquaIx_equilibrium: calculate distribution of time at different aquatic habitats during oviposition from MBITESi.basic.move
aquaIx_equilibrium <- function(mosyPop){

  aquaIx = sapply(LANDSCAPE$aquaSites,function(x){x$ix}) # vector of aquatic sites

  aquaVecAll = parallel::mclapply(X = mosyPop,FUN = aquaIx_oneMosy, aquaIx = aquaIx)
  aquaVec = Reduce(f = "+",x = aquaVecAll)

  return(aquaVec/sum(aquaVec))
}

# aquaIx_oneMosy: calculate count vector of visits to each aquatic habitat:
# function needs to be vectorized to run over the mosyPop list.
aquaIx_oneMosy <- function(oneMosy, aquaIx){
  aquaVec = rep(0,length = length(aquaIx))

  with(oneMosy$history,{
    if(!"O" %in% stateH){
      return(aquaVec)
    } else {
      oviIx = which(stateH == "O")
      for(ix in oviIx){
        aquaVec[ixH[oviIx]] = aquaVec[ixH[oviIx]] + 1
      }
      return(aquaVec)
    }
  })

}


##########################################
# Plot population dynamics
##########################################

#' Plot Imported Adult Dynamics .csv Data
#'
#' This is a utility to plot logged adult dynamics data after it has been imported by \code{\link{importAdults}}.
#'
#' @param egg a data frame (output of \code{\link{importAdults}})
#' @return plot
#' @examples
#' plotAdults(adults)
plotAdults <- function(adults){

  par(mfrow=c(1,2))

  # female parameters
  fStates = colnames(adults)[grep("f$",x = colnames(adults))]
  fStatesLegend = substr(x = fStates,start = 1,stop = 1)
  fCol = ggCol(n = length(fStates),alpha = 0.75)
  fColSmooth = ggCol(n = length(fStates))
  names(fColSmooth) = fStates

  # plot female dynamics
  matplot(x = adults[,fStates],type = "l",col = fCol,lty = 1:length(fStates),
          ylab="Adult Densities",xlab="Time (Days)",main = "Females")
  legend("topleft",legend = fStatesLegend,col = fCol,lty = 1:length(fStates),bty="n")

  # smoothed female dynamics
  for(ix in fStates){
    stateSmooth = ksmooth(x = adults$time,y = adults[,ix],bandwidth = max(adults$time)/20)
    lines(x = stateSmooth$x,y = stateSmooth$y,col = fColSmooth[ix],lwd = 1.5)
  }

  # male parameters
  mStates = colnames(adults)[grep("m$",x = colnames(adults))]
  mStatesLegend = substr(x = mStates,start = 1,stop = 1)
  mCol = ggCol(n = length(mStates),alpha = 0.75)
  mColSmooth = ggCol(n = length(mStates))
  names(mColSmooth) = mStates

  # plot male dynamics
  matplot(x = adults[,mStates],type = "l",col = mCol,lty = 1:length(mStates),
          ylab="Adult Densities",xlab="Time (Days)",main = "Males")
  legend("topleft",legend = mStatesLegend,col = mCol,lty = 1:length(mStates),bty="n")

  # smoothed male dynamics
  for(ix in mStates){
    stateSmooth = ksmooth(x = adults$time,y = adults[,ix],bandwidth = max(adults$time)/20)
    lines(x = stateSmooth$x,y = stateSmooth$y,col = mColSmooth[ix],lwd = 1.5)
  }

  par(mfrow=c(1,1))

}
