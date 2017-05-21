#################################################################
#
#   MASH/MBITES
#   Human Profile
#   R version
#   Sean Wu
#   January 24, 2017
#
#################################################################


##########################################
# Generic Event History Management Functions
##########################################

#' Write eventT and events from HUMANS to .json
#'
#' Write HUMANS histories from PfSI module to .json file. This should typically be called at the end of a model run, assuming demography is not implemented.
#'
#' @param directory directory; files will be put in directory/OUTPUT/..
#' @param fileName name of the file to write to; directory/OUTPUT/fileName.csv
#' @param overWrite allow overwriting of files in directory/OUTPUT/..? Use with caution. Generally it is reccomended to use \code{\link{clearOutput}} to clear the OUTPUT/.. folder after moving data to a permanant location and prior to running the simulation.
#' @return nothing
#' @examples
#' writeHumanEvent(directory, fileName overWrite = FALSE)
writeHumanEvent_PfSI <- function (directory, fileName, overWrite = FALSE)
{
  if (!dir.exists(paste0(directory, "OUTPUT"))) {
    dir.create(paste0(directory, "OUTPUT"), recursive = T)
  }
  if (!overWrite & file.exists(paste0(directory, "OUTPUT/", fileName))) {
    stop("writeHumanEvent_PfSI cannot write to a file that already exists!")
  }
  con = file(description = paste0(directory, "OUTPUT/", fileName),
             open = "wt")
  humanHistories = lapply(HUMANS, function(x) {
    list(eventT = x$Pathogens$Pf$eventT, events = x$Pathogens$Pf$events)
  })
  humanID = sapply(HUMANS, function(x) {
    x$myID
  })
  names(humanHistories) = paste0("human", humanID)
  writeLines(text = jsonlite::toJSON(x = humanHistories, pretty = TRUE),
             con = con)
  close(con)
}

# writeHumanEvent_PfSI <- function(directory, fileName){
#   if(!dir.exists(paste0(directory,"OUTPUT"))){
#     dir.create(paste0(directory))
#     dir.create(paste0(directory,"OUTPUT"))
#   }
#   if(file.exists(paste0(directory,"OUTPUT/",fileName))){
#     stop("writeHumanEvent_PfSI cannot write to a file that already exists!")
#   }
#
#   con = file(description = paste0(directory,"OUTPUT/",fileName),open = "wt")
#
#   humanHistories = lapply(HUMANS,function(x){
#     list(
#       eventT = x$Pathogens$Pf$eventT,
#       events = x$Pathogens$Pf$events
#     )
#   })
#   humanID = sapply(HUMANS,function(x){x$myID})
#   names(humanHistories) = paste0("human",humanID)
#   writeLines(text = jsonlite::toJSON(x = humanHistories,pretty = TRUE),con = con)
#   close(con)
# }



#' Reimport Human eventT and events from .json
#'
#' Import human histories exported to .json from PfSI module.
#'
#' @param directory directory; files are in directory/OUTPUT/..
#' @return list of histories
#' * eventT: times events occur
#' * events: events
#' @md
#' @examples
#' writeHumanEvent(directory, fileName)
importHumanEvent_PfSI <- function(directory){
  dirFiles = system(command = paste0("ls ",directory,"OUTPUT/"),intern = TRUE)
  humanFile = grep("humanPfSI.json",dirFiles)
  humanHistories = jsonlite::fromJSON(txt = paste0(directory,"OUTPUT/",dirFiles[humanFile]))
  return(humanHistories)
}

#' Extract PfSI History from Human Object
#'
#' Essentially go through the process of exporting via \code{\link{writeHumanEvent_PfSI}} and importing via \code{\link{importHumanEvent_PfSI}}
#' without exporting data to .json to return formatted human event histories.
#'
#' @param HUMANS the \code{HUMANS} object
#' @return list of histories
#' * eventT: times events occur
#' * events: events
#' @md
#' @examples
#' writeHumanEvent(directory, fileName)
convertHumanEvent_PfSI <- function(HUMANS){
  humanHistories = lapply(HUMANS,function(x){
    list(
      eventT = x$Pathogens$Pf$eventT,
      events = x$Pathogens$Pf$events
    )
  })
  humanID = sapply(HUMANS,function(x){x$myID})
  names(humanHistories) = paste0("human",humanID)
  return(humanHistories)
}


##########################################
# PfSI Trajectories
##########################################

#' Plot a Single Human Trajectory for PfSI Module
#'
#' Plots a single trajectory through state space for PfSI module.
#'
#' @param ixH index of human
#' @param oneHistory a single human histoy \code{history[[ixH]]}
#' @param tMax maximum time event observed in all histories
#' @return plot
#' @examples
#' plotPfsiOneTrajectory
plotPfsiOneTrajectory <- function(ixH, oneHistory, tMax){

  times = oneHistory$eventT[-1]
  events = oneHistory$events[-1]

  # fever
  ixF = which(events == "F")
  if(length(ixF)>0){
    points(x = times[ixF],y = rep(x = ixH,times = length(ixF)),pch=17,cex=0.5,col="red")
  }

  # prophylaxis
  ixP = which(events == "P")
  if(length(ixP)>0){
    points(x = times[ixP],y = rep(x = ixH,times = length(ixP)),pch=16,cex=0.5,col="blue")
  }

  events = events[-ixF]
  times = times[-ixF]

  # state trajectory
  if(length(events) > 0){

    for(i in 1:(length(events))){

      if(i == length(events)){ # final trajectory
        if(events[i] == "I"){
          segments(x0 = times[i],y0 = ixH,x1 = tMax,y1 = ixH,col = "red",lwd = 2)
        }
        if(events[i] == "S"){
          segments(x0 = times[i],y0 = ixH,x1 = tMax,y1 = ixH,col = "grey50",lwd = 2)
        }
        if(events[i] == "P"){
          segments(x0 = times[i],y0 = ixH,x1 = tMax,y1 = ixH,col = "blue",lwd = 2)
        }

      } else { # interior trajectory
        if(events[i] == "I"){
          segments(x0 = times[i],y0 = ixH,x1 = times[i+1],y1 = ixH,col = "red",lwd = 2)
        }
        if(events[i] == "S"){
          segments(x0 = times[i],y0 = ixH,x1 = times[i+1],y1 = ixH,col = "grey50",lwd = 2)
        }
        if(events[i] == "P"){
          segments(x0 = times[i],y0 = ixH,x1 = times[i+1],y1 = ixH,col = "blue",lwd = 2)
        }
      }

    } # end for
  } # end state trajectory

}

#' Plot Human Trajectories from PfSI Module Histories
#'
#' Plots all human trajectories through state space for PfSI module.
#'
#' @param history imported human histories from \code{\link{importHumanEvent_PfSI}}.
#' @return plot
#' @examples
#' plotPfsiOneTrajectory
plotPfsiTrajectory <- function(history){

  tMax = max(sapply(X = history,FUN = function(x){max(x$eventT)}))

  plot(1,type="n",xaxt="n",yaxt="n",ylab="Humans",xlab="Time (Years)",xlim=c(0,tMax),ylim=c(0,length(history)))
  ttMax = tMax/365
  axis(side = 1,at = c(0:ttMax)*365,labels = c(0:ttMax))

  for(ixH in 1:length(history)){
    plotPfsiOneTrajectory(ixH = ixH,oneHistory = history[[ixH]],tMax = tMax)
  }

}
