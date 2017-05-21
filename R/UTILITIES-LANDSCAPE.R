#################################################################
#
#   MASH/MBITES
#   Landscape/Movement Profile
#   R version
#   Sean Wu
#   January 24, 2017
#
#################################################################


##########################################
# Auxiliary Functions
##########################################

#helper function to sample equally from color space
ggCol <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

ggColOffset <- function(n,offset){
  hues = seq(15, 375, length = n + 1 + offset)
  hcl(h = hues,l = 65,c = 100)[(1+offset):(offset+n)]
}

#helper function to brighten or darken colors
colLuminosity <- function(color,factor,bright,alpha=NULL){

  if(!is.logical(bright)){ #sanity check
    stop("i don't know if you wan't me to make your color brighter or darker!")
  }

  col = col2rgb(color,alpha=FALSE) #convert to rgba color space
  if(bright){
    col = col*factor
  } else {
    col = col/factor
  }
  if(!is.null(alpha)){ #adjust alpha if specified
    rbind(col,alpha = 255 * alpha)
  }
  col = rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  return(col)
}


##########################################
# Static Landscape Plot
##########################################

#setup_plot will print the static landscape surface
#offset controls the offset on maximum/minimum XY coordinates due to plotting error
#bgcol is the background color of the plot; default is #f2ffe6
staticLandscapePlot <- function(offset=0.05,bgCol="grey20",cex=0.75){

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

  points(x = feedXY, pch = 15, col = setup_col[1], cex = cex) #draw feeding sites
  points(x = aquaXY, pch = 16, col = setup_col[2], cex = cex) #draw aquatic habitats
  points(x = sugarXY, pch = 17, col = setup_col[3], cex = cex) #draw sugar sources
  points(x = mateXY, pch = 18, col = setup_col[4], cex = cex) #draw mating sites

  legend(x = "topleft", legend = c("Feeding Site","Aquatic Habitat","Sugar Source","Mating Site"),pch = 15:18,
         col = setup_col, bty = "n", text.col = "grey80")

  par(bg=.GlobalEnv$defaultPar$bg,mar=.GlobalEnv$defaultPar$mar,mgp=.GlobalEnv$defaultPar$mgp)
}


##########################################
# MvOb Kernels
##########################################

#showKernelsPlot: generate plot of movement kernels
#S: MvOb of starting sites
#D: MvOb of destination sites
#N: number of random starting sites
#dMax: maximum distance for distance bins
#dMesh: granularity of distance bins
#mtl: optional title
showKernelsPlot <- function(S,D,N=50, dMesh = 0.01, dMax=NULL, mtl = NULL){

  if(exists("defaultPar",.GlobalEnv)){ #reset plot parameters
    suppressWarnings(par(defaultPar))
  }

  M1 = powerKernel(S,D) #Markov transition matrix
  M2 = distanceMat(S,D) #distance matrix
  M1n = M1/sum(M1) #normalized transition matrix

  if(is.null(dMax)){
    dMax = max(M2) + max(M2)*0.05
  }

  xCDF = function(x, M2, M1){ #return movement CDF by distance bins
    sum(M1[which(M2<=x, arr.ind = TRUE)])
  }

  x=seq(0, dMax, by = dMesh) #distance bins
  cdfT = sapply(X = x,FUN = xCDF, M2=M2, M1=M1n) #generate movement CDF
  d = length(S)
  cdf = matrix(0,N,length(x))
  for(i in 1:N){
    j = sample(x = 1:d,size = 1) #pick a random starting site
    cdf[i,] = sapply(x, xCDF, M2=M2[j,], M1=M1[j,])
  }

  xmx = max(cdfT[1], diff(cdfT), cdf[2,]-cdf[1,], cdf[3,]-cdf[2,])

  par(mfrow=c(1,2))
  colors = viridis(n=N)
  plot(x, cdfT, type = "l", xlab = "Distance", ylab = "CDF", ylim = c(0,1), main = mtl)
  grid()
  cdf = cdf[order(apply(cdf,MARGIN = 1,which.max)),] #sort rows by probability of long-range movement
  for(i in 1:N){
    lines(x, cdf[i,], col = colors[i])
  }
  lines(x,cdfT, lwd =3)

  #need to plot the max of the y values as ylin.
  plot(x, c(cdfT[1],diff(cdfT)), xlab = "Distance", ylab = "PDF", type = "l", ylim = c(0, 1), main = mtl)
  grid()
  for(i in 1:N){
    lines(x, c(cdf[i,1],diff(cdf[i,])), col = colors[i])
  }
  lines(x,c(cdfT[1],diff(cdfT)), lwd =3)
  par(mfrow=c(1,1))
}


##########################################
# Aquatic Ecology
##########################################


emergeSeasonsPlot <- function(){
  pcol = viridis(n = LANDSCAPE$nA)
  ylim = range(sapply(LANDSCAPE$aquaSites,function(x){x$season}))
  plot(1:365,type="n",ylim = ylim, ylab = expression(lambda[i]), xlab = "Time (Days)")
  grid()
  for(ix in 1:LANDSCAPE$nA){
    lines(1:365, LANDSCAPE$aquaSites[[ix]]$season, col=pcol[ix])
  }
}
