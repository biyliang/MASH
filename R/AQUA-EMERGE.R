#################################################################
#
#   MASH/MBITES
#   Aquatic Ecology
#   Emergence
#   R Version
#   Sean Wu
#   February 2, 2017
#
#################################################################


#oneDay_emerge: add emerging adults to ImagoQ
#tNow: current daily time tick of model
oneDay_emerge <- function(tNow){

  # log EggQ
  if(EggQ_TRACK){
    trackEggQ(con = .GlobalEnv$EggQCon)
  }

  lambdaExact = sapply(LANDSCAPE$aquaSites,function(x){x$season[floor(tNow)%%365+1]}) #exact emergence rates

  nSites = length(lambdaExact)
  lambdaEmerge = rpois(n = nSites,lambda = lambdaExact)

  #add emerging adults to ImagoQ
  for(ix in 1:LANDSCAPE$nA){
    addAdults2Q(lambda = lambdaEmerge[ix],t = tNow,ix = ix,dam = 0,sire = 0)
  }

}


#makeAquaPop_emerge: generate seasonally forced emergence
#LANDSCAPE: the landscape object
#lambda: net yearly average number of emerging females per human per day (over entire landscape)
#offset: seasonal offset
makeAquaPop_emerge <- function(lambda, offset = 0){

  w = rgamma(n = LANDSCAPE$nA,1,1)
  K = lambda*w / sum(w)
  offset = rep(offset,length=LANDSCAPE$nA)

  for(ix in 1:LANDSCAPE$nA){
    LANDSCAPE$aquaSites[[ix]]$season <<- K[ix]*(1+sin(2*pi*(c(1:365)-offset[ix])/365))
  }
}
