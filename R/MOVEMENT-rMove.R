#################################################################
#
#   MASH/MBITES
#   Movement object generation routines
#   R version
#   Sean Wu
#   January 25, 2017
#
#################################################################

##########################################
# Auxiliary Functions
##########################################
distanceMat <- function(S,D){
  #. distMat: Calculate Euclidean distance matrix
  #S: starting object (list)
  #D: destination object (list)
  Sxy = t(sapply(S,function(x){x$siteXY}))
  Dxy = t(sapply(D,function(x){x$siteXY}))
  dMat = matrix(0,nrow = length(S),ncol = length(D))

  for(i in 1:nrow(Sxy)){
    for(j in 1:nrow(Dxy)){
      dMat[i,j] = sqrt((Sxy[i,1]-Dxy[j,1])^2 + (Sxy[i,2]-Dxy[j,2])^2)
    }
  }
  return(dMat)
}

######################################################################
#  Compute Probability Vectors and Matrices
#  note: the *S2D family of functions are generic.
######################################################################
powerKernel <- function(S,D,sigma=3,eps=0.1,beta=0){
  #. powerKernel: Compute Markov transition matrix between S and D
  dW = sapply(D,function(x){x$w})
  dS2D = distanceMat(S,D)
  S2D = matrix(0,nrow=length(S),ncol=length(D))

  for(i in 1:length(S)){
    allP = dW^(-beta*dS2D[i,]) * (eps + dS2D[i,])^-sigma
    S2D[i,] = allP / sum(allP)
  }

  return(S2D)
}

prSort <- function(id,pr){
  #.prSort: Return sorted list of ids and probabilities
  ot = order(pr,decreasing = TRUE)
  return(list(
    id = id[ot],
    pr = pr[ot]/sum(pr)
  ))
}


############################################################
#  Movement Object :: Exact
#  exactMvOb uses exact probabilities
#  note: Best when ||f|| and ||l|| are small
#  the exactMvOb will be an list of lists
############################################################
exactMvOb <- function(S,D,sigma=3,eps=0.1,beta=0){
  #. exactMvOb: Calculate the movement object between S and D
  ixS = sapply(S,function(x){x$ix}) #id of site in S
  prS2D = powerKernel(S,D,sigma,eps,beta) #movement matrix between S to D
  MvOb = vector(mode="list",length = length(S)) #empty movement object
  nD = length(D) #number of D sites

  if(identical(S,D)){ # movement within same class of site
    for(i in 1:length(S)){
      id = (1:nD)[-i]
      pr = prS2D[i,-i]
      sortedPr = prSort(id,pr)
      MvOb[[i]] = list(
        ix = ixS[i],
        PR = c(0,1,0), # PR = c(pr[i],1-pr[i],0),
        near = sortedPr,
        around = NULL,
        moveFar = NULL
      )
    }
  } else { # movement between classes of sites
    for(i in 1:length(S)){
      sortedPr = prSort(1:nD,prS2D[i,])
      MvOb[[i]] = list(
        ix = ixS[i],
        PR = c(0,1,0), # PR = c(pr[i],1-pr[i],0),
        near = sortedPr,
        around = NULL,
        moveFar = NULL
      )
    }
  }

  return(MvOb)
}

exactAll <- function(landscape,sigma=3,eps=0.1,beta=0){
  #. exactAll: Calculate all MvOb for site classes
  list(

    #move to feeding site
    F2F = exactMvOb(S = landscape$feedSites,D = landscape$feedSites,sigma,eps,beta),
    S2F = exactMvOb(S = landscape$sugarSites,D = landscape$feedSites,sigma,eps,beta),
    M2F = exactMvOb(S = landscape$swarmSites,D = landscape$feedSites,sigma,eps,beta),
    L2F = exactMvOb(S = landscape$aquaSites,D = landscape$feedSites,sigma,eps,beta),

    #move to aquatic habitat
    F2L = exactMvOb(S = landscape$feedSites,D = landscape$aquaSites,sigma,eps,beta),
    S2L = exactMvOb(S = landscape$sugarSites,D = landscape$aquaSites,sigma,eps,beta),
    M2L = exactMvOb(S = landscape$swarmSites,D = landscape$aquaSites,sigma,eps,beta),
    L2L = exactMvOb(S = landscape$aquaSites,D = landscape$aquaSites,sigma,eps,beta),

    # move to sugar site
    F2S = exactMvOb(S = landscape$feedSites,D = landscape$sugarSites,sigma,eps,beta),
    S2S = exactMvOb(S = landscape$sugarSites,D = landscape$sugarSites,sigma,eps,beta),
    M2S = exactMvOb(S = landscape$swarmSites,D = landscape$sugarSites,sigma,eps,beta),
    L2S = exactMvOb(S = landscape$aquaSites,D = landscape$sugarSites,sigma,eps,beta),

    # move to swarming site
    F2M = exactMvOb(S = landscape$feedSites,D = landscape$swarmSites,sigma,eps,beta),
    S2M = exactMvOb(S = landscape$sugarSites,D = landscape$swarmSites,sigma,eps,beta),
    M2M = exactMvOb(S = landscape$swarmSites,D = landscape$swarmSites,sigma,eps,beta),
    L2M = exactMvOb(S = landscape$aquaSites,D = landscape$swarmSites,sigma,eps,beta)

  )

}


###########################################
# Movement function
###########################################


rMove <- function(ix, pSet, bState){
  #. rMove: Based on current location and biological stage, move a mosquito
  #ix: location index of mosquito
  #pSet: point set index of mosquito {f,l}
  #bState: behavioral state of the mosquito
  #MvMat: the movement object

  if(!bState %in% c("F","L","S","M")){
    return(list(ixNew=ix,pSetNew=pSet))
  } else if(bState == "F"){
    if(pSet == "f"){MvMat = MvOb$F2F[[ix]]}
    if(pSet == "s"){MvMat = MvOb$S2F[[ix]]}
    if(pSet == "m"){MvMat = MvOb$M2F[[ix]]}
    if(pSet == "l"){MvMat = MvOb$L2F[[ix]]}
    pSetNew = "f"
  } else if(bState == "L"){
    if(pSet == "f"){MvMat = MvOb$F2L[[ix]]}
    if(pSet == "s"){MvMat = MvOb$S2L[[ix]]}
    if(pSet == "m"){MvMat = MvOb$M2L[[ix]]}
    if(pSet == "l"){MvMat = MvOb$L2L[[ix]]}
    pSetNew = "l"
  } else if(bState == "S"){
    if(pSet == "f"){MvMat = MvOb$F2S[[ix]]}
    if(pSet == "s"){MvMat = MvOb$S2S[[ix]]}
    if(pSet == "m"){MvMat = MvOb$M2S[[ix]]}
    if(pSet == "l"){MvMat = MvOb$L2S[[ix]]}
    pSetNew = "s"
  } else if(bState == "M"){
    if(pSet == "f"){MvMat = MvOb$F2M[[ix]]}
    if(pSet == "s"){MvMat = MvOb$S2M[[ix]]}
    if(pSet == "m"){MvMat = MvOb$M2M[[ix]]}
    if(pSet == "l"){MvMat = MvOb$L2M[[ix]]}
    pSetNew = "m"
  } else {
    browser("unrecognized behavioral state")
  }

  x = runif(1)

  ixNew = with(MvMat,{
    if(x <= PR[1]){ #no movement
      ix
    } else {
      if(x <= PR[1] + PR[2]){ #near movement
        ixNear = sample(x = 1:length(near$id),size = 1,prob = near$pr)
        near$id[ixNear]
      } else {
        if(x <= sum(PR)){ #around movement
          browser("'around' movement not yet implemented")
        } else { #moveFar movement
          browser("'moveFar' movement not yet implemented")
        }
      }
    }
  })

  return(list(ixNew=ixNew,pSetNew=pSetNew))
}
