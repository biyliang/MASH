#################################################################
#
#   MASH/MBITES
#   Visualization and Analysis of logged data for AQUATIC ECOLOGY modules
#   R version
#   Sean Wu
#   April 6, 2017
#
#################################################################


#################################################################
# EL4P Basic
#################################################################


#' Plot Imported EL4P .csv Data
#'
#' This is a utility to plot logged EL4P and EggQ data after it has been imported by importEL4P
#'
#' @param el4p a data frame (output of \code{\link{importEL4P}})
#' @param egg a data frame (output of \code{\link{importEggQ}}; may be set to NULL)
#' @return plot
#' @examples
#' formatEL4P(el4p)
plotEL4P <- function(el4p, egg = NULL){

  par(mfrow=c(2,3))

  N = ncol(el4p)-1 # number of sites

  col = ggCol(n = N,alpha = 0.8)

  if(!is.null(egg)){
    matplot(x = egg[,-which(names(egg)=="time")],type = "l",col = col,lty = 1,
            ylab="Eggs",xlab="Time (Days)")
    eggFlat = as.vector(as.matrix(egg[,-which(names(egg)=="time")]))
    eggTime = rep(x = egg$time,times = ncol(egg)-1)
    nwSmooth = ksmooth(x = eggTime,y=eggFlat,bandwidth = max(egg$time)/100)
    lines(x = nwSmooth$x,y = nwSmooth$y,col = "grey20",lwd=2.5,lty=1)
  }

  # L1 stage
  L1 = el4p[el4p$labels=="L1",-1]
  L1$time = 1:nrow(L1)
  matplot(x = L1[,-which(names(L1)=="time")],type = "l",col = col,lty = 1,
          ylab="Larval instar stage 1 (L1)",xlab="Time (Days)")
  L1flat = as.vector(as.matrix(L1[,-which(names(L1)=="time")]))
  L1time = rep(x = L1$time,times = ncol(L1)-1)
  nwSmooth = ksmooth(x = L1time,y=L1flat,bandwidth = max(L1$time)/100)
  lines(x = nwSmooth$x,y = nwSmooth$y,col = "grey20",lwd=2.5)

  # L2 stage
  L2 = el4p[el4p$labels=="L2",-1]
  L2$time = 1:nrow(L2)
  matplot(x = L2[,-which(names(L2)=="time")],type = "l",col = col,lty = 1,
          ylab="Larval instar stage 2 (L2)",xlab="Time (Days)")
  L2flat = as.vector(as.matrix(L2[,-which(names(L2)=="time")]))
  L2time = rep(x = L2$time,times = ncol(L2)-1)
  nwSmooth = ksmooth(x = L2time,y=L2flat,bandwidth = max(L2$time)/100)
  lines(x = nwSmooth$x,y = nwSmooth$y,col = "grey20",lwd=2.5)

  # L3 stage
  L3 = el4p[el4p$labels=="L3",-1]
  L3$time = 1:nrow(L3)
  matplot(x = L3[,-which(names(L3)=="time")],type = "l",col = col,lty = 1,
          ylab="Larval instar stage 3 (L3)",xlab="Time (Days)")
  L3flat = as.vector(as.matrix(L3[,-which(names(L3)=="time")]))
  L3time = rep(x = L3$time,times = ncol(L3)-1)
  nwSmooth = ksmooth(x = L3time,y=L3flat,bandwidth = max(L3$time)/100)
  lines(x = nwSmooth$x,y = nwSmooth$y,col = "grey20",lwd=2.5)

  # L4 stage
  L4 = el4p[el4p$labels=="L4",-1]
  L4$time = 1:nrow(L4)
  matplot(x = L4[,-which(names(L4)=="time")],type = "l",col = col,lty = 1,
          ylab="Larval instar stage 4 (L4)",xlab="Time (Days)")
  L4flat = as.vector(as.matrix(L4[,-which(names(L4)=="time")]))
  L4time = rep(x = L4$time,times = ncol(L4)-1)
  nwSmooth = ksmooth(x = L4time,y=L4flat,bandwidth = max(L4$time)/100)
  lines(x = nwSmooth$x,y = nwSmooth$y,col = "grey20",lwd=2.5)

  # P stage
  P = el4p[el4p$labels=="P",-1]
  P$time = 1:nrow(P)
  matplot(x = P[,-which(names(P)=="time")],type = "l",col = col,lty = 1,
          ylab="Pupae (P)",xlab="Time (Days)")
  Pflat = as.vector(as.matrix(P[,-which(names(P)=="time")]))
  Ptime = rep(x = P$time,times = ncol(P)-1)
  nwSmooth = ksmooth(x = Ptime,y=Pflat,bandwidth = max(P$time)/100)
  lines(x = nwSmooth$x,y = nwSmooth$y,col = "grey20",lwd=2.5)

  par(mfrow=c(1,1))

}
