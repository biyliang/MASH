#################################################################
#
#   MASH/MBITES
#   Aquatic Ecology
#   EL4P Fitting Utilities
#   R Version
#   Sean Wu
#   March 22, 2017
#
#################################################################

# getLambda: return average lambda required to sustain given R0
# R0: desired value of R0 to calibrate against
# summary: output from cohortBionomics() or MBITES-BASIC NEED TO WRITE
# nH: number of humans
# P: M-BITES parameter list
getLambda <- function(R0, summary, nH, P){
  with(c(P,summary),{
    P = exp(-EIP/lifespanC) # probability to survive EIP
    S = feedAllC # stability index
    (R0 * nH * (S^2) * r) / (b * c * P)
  })
}

#################################################################
#
# makePAR_EL4P: make PAR for EL4P
#
# arguments:
# nA: number of aquatic habitats
# nH: number of humans
# R0: desired R0 at equilibrium
# M: vector of densities of mosquitoes at each site ix
# aquaEq: vector of probability of oviposition from aquaIx.equilibrium(...)
# par: parameters from MBITES-Parameters
# summary: cohort summary bionomics from cohortBionomics(...)
# p: 1/p is expected time until advance to next stage
# P: density independent survival from egg to adult
# a: shape parameter for gamma distributed weights
# b: scale parameter for gamma distributed weights
#
# computed values:
# G: mean total lifetime egg production of an adult female
# nu: mean batch size for one oviposition
# lambda: value of adult female emergence to sustain R0
#
#################################################################
makePAR_EL4P <- function(nA, nH, R0, M, aquaEq, par, summary, p = 0.9, P = 0.8, a = 1, b = 1){

  PAR = list(
      M = M,
      aquaEq = aquaEq,
      nA = nA,
      nH = nH,
      lifespan = summary$lifespanC,
      G = summary$tBatchC,
      nu = summary$mBatchC,
      p = p,
      P = P,
      a = a,
      b = b,
      lambda = getLambda(R0, summary, nH, par)
    )

  return(PAR)
}

#################################################################
#
# setupAquaPop_EL4P: fit values of psi for each aquatic habitat
#
# arguments:
# nA: number of aquatic habitats
# nH: number of humans
# PAR: parameters from makePAR_EL4P
#
# computed values:
# W: weights for K
# K: carrying capacities
# pp: mean of initial values for alpha
#
#
#
#################################################################

setupAquaPop_EL4P <- function(PAR, tol = 0.1, plot = FALSE){
  with(PAR,{

    # calculate initial parameter values
    W = rgamma(n = nA,shape = a,scale = b)
    K = (lambda*W) / sum(W)
    PAR$K = K # attach K to PAR
    pp = -log(P^((1-p)/5))
    alpha = abs(rnorm(n = nA,mean = pp,sd = 0.004))
    PAR$alpha = alpha # attach alpha to PAR
    psi = alpha/K
    PAR$psiInit = psi # attach psi to PAR

    # EL4P populations
    EL4P_pops = replicate(n = nA,expr = EL4P(),simplify = FALSE)

    # fit psi
    rng = range(K)
    AquaPops = meshK_EL4P(lK = rng[1],uK = rng[2],EL4P_pops = EL4P_pops,PAR = PAR,plot = plot,tol = tol)
    PAR$psiOptim = AquaPops$psiHat # attach psi fitted by optimize(...) to PAR
    PAR$meshK = AquaPops$meshK # attach sampled meshK to PAR

    if(plot){plotPsi(PAR = PAR)}
    if(plot){par(mfrow=c(1,2))}
    cf = psi2K_cf(AquaPops$meshK,AquaPops$psiHat,plot)
    if(plot){
      psi2K_plot(lmFit = lm(1/K2psi(AquaPops$meshK,cf)~AquaPops$meshK+0),K = AquaPops$meshK,psi = K2psi(AquaPops$meshK,cf))
      par(mfrow=c(1,1))
    }
    psiHat = K2psi(AquaPops$meshK,cf)
    PAR$psiHat = psiHat # attach psi fitted via linear regression to PAR

    # run all EL4P aquatic populations to equilibrium values
    EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = psiHat,pop = AquaPops$EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
                                   mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)
    return(
      list(EL4P_pops=EL4P_pops,
           PAR=PAR
      )
    )
  })
}

# sample from arbitrary mesh of points for K; above function is for exact calculations
# gridN: number of points to sample
setupAquaPop_EL4P_samplePoints <- function(PAR, gridN, tol = 0.1, plot = FALSE){
  with(PAR,{

    # calculate initial parameter values
    W = rgamma(n = nA,shape = a,scale = b)
    K = (lambda*W) / sum(W)
    PAR$K = K # attach K to PAR
    pp = -log(P^((1-p)/5))
    alpha = abs(rnorm(n = nA,mean = pp,sd = 0.004))
    PAR$alpha = alpha # attach alpha to PAR
    psi = alpha/K
    PAR$psiInit = psi # attach psi to PAR

    # EL4P populations
    EL4P_pops = replicate(n = nA,expr = EL4P(),simplify = FALSE)

    # fit psi
    rng = range(K)
    AquaPops = meshK_EL4P(lK = rng[1],uK = rng[2],EL4P_pops = EL4P_pops,PAR = PAR,plot = plot,tol = tol)
    PAR$psiOptim = AquaPops$psiHat # attach psi fitted by optimize(...) to PAR
    PAR$meshK = AquaPops$meshK # attach sampled meshK to PAR

    if(plot){plotPsi(PAR = PAR)}

    cf = psi2K_cf(AquaPops$meshK,AquaPops$psiHat,plot)
    psiHat = K2psi(AquaPops$meshK,cf)
    PAR$psiHat = psiHat # attach psi fitted via linear regression to PAR

    # run all EL4P aquatic populations to equilibrium values
    EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = psiHat,pop = AquaPops$EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
                                   mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)
    return(
      list(EL4P_pops=EL4P_pops,
           PAR=PAR
      )
    )
  })
}


#################################################################
#
# EL4P Fitting Routines
# set values of psi so that lambda = K at (p,G)
# G :: lifetime egg production, per adult (tBatchC from MBITES-BASIC)
#
#################################################################



# meshK_EL4P: generate mesh of K values then fit psi to mesh; interpolate inbetween discrete sampled values of K
# lK: lower bound of K
# uK: upper bound of K
# alpha: vector of alpha parameters
# EL4P_pops: list of EL4P populations
# PAR: passed from Aq.PAR
# tol:
meshK_EL4P <- function(lK, uK, EL4P_pops, PAR, tMax = 500, plot = FALSE, tol = 0.1){

  # sample K on mesh in log space; transform to normal space
  meshK = exp(seq(log(lK),log(uK),length.out=length(EL4P_pops)))
  psiHat = vector(mode="numeric",length=length(EL4P_pops))

  # fit EL4P: set values of psi so lambda = K at (p,G)
  for(ix in 1:length(EL4P_pops)){
    print(paste0("fitting psi for site ix: ",ix))
    psiMin = optimize(f = psiFit,interval = c(0,10), ix=ix, ixEL4P=EL4P_pops[[ix]], ixKmesh=meshK[ix], PAR=PAR)
    psiHat[ix] = psiMin$minimum
  }

  # run all EL4P aquatic populations to equilibrium values
  print(paste0("run aquatic populations to equilibrium values"))
  EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = psiHat,pop = EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
                                 mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)

  return(
    list(
      EL4P_pops = EL4P_pops, # aquatic populations at equilibrium values
      psiHat = psiHat, # fitted values of psi
      meshK = meshK # sampled values of K used to fit psi
    )
  )
}

# psiFit: fit psi for a single aquatic habitat;
# this function does not directly fit psi itself; but it is the loss function to be passed to optimize(...)
# x: to be passed to optim(...)
# ix: the index of this habitat
# ixEL4P: the element of EL4P_pops[[ix]] for this habitat
# ixKmesh: the value of meshK for habitat ix
# ixAlpha: the value of alpha for habitat ix
# PAR: passed from Aq.PAR
psiFit <- function(x, ix, ixEL4P, ixKmesh, PAR){
  psi = abs(x)
  EL4P = runOne_GEL4P(ix = ix,psi = psi,pop = ixEL4P,PAR = PAR)
  return((EL4P$lambda - ixKmesh)^2) # return squared error of empirical lambda around value of meshK
}

# runOne_GEL4P: run single site of EL4P dynamics for psiFit
# ix: index of site
# psi: value of psi for site
# alpha: value of alpha for site
# pop: EL4P population
# PAR: Aq.PAR
# tMax: time to run aquatic dynamics
runOne_GEL4P <- function(ix, psi, pop, PAR, tMax = 150){
  for(i in 1:30){
    pop = oneDay_GEL4P(ix = ix,psi = psi,pop = pop,PAR = PAR)
  }
  lambdaH = vector(mode="numeric",length = tMax+1)
  lambdaH[1] = pop$lambda # record values of lambda
  for(i in 1:tMax){
    pop = oneDay_GEL4P(ix = ix,psi = psi,pop = pop,PAR = PAR)
    PAR$M[ix] = ((exp(-1/PAR$lifespan))*PAR$M[ix]) + pop$lambda # simulate adult population dynamics
    lambdaH[i+1] = pop$lambda
  }
  # plot(lambdaH, type = "l", main = paste0("lambdaH for site: ",ix), col = "red", lwd = 2)
  return(pop)
}

# oneDay_GEL4P: run single step of EL4P dynamics for site ix
# ix: index of site
# psi: value of psi for site
# alpha: value of alpha for site
# pop: EL4P population
# PAR: Aq.PAR
oneDay_GEL4P <- function(ix, psi, pop, PAR){
  with(PAR,{
    L1o = pop$L1; L2o=pop$L2; L3o=pop$L3; L4o=pop$L4
    D   = sum(L1o+L2o+L3o+L4o)
    s1  = exp(-alpha[ix])
    s2  = exp(-(alpha[ix]+ psi*D))
    pop$lambda = s1*pop$P
    pop$P  = s2*p*L4o
    pop$L4 = s2*(p*L3o + (1-p)*L4o)
    pop$L3 = s2*(p*L2o + (1-p)*L3o)
    pop$L2 = s2*(p*L1o + (1-p)*L2o)
    pop$L1 = pop$eggs + s2*(1-p)*L1o
    pop$eggs = M[ix]*aquaEq[ix]*(G/lifespan)
    return(pop)
  })
}

# run2Eq_GEL4P: run a single EL4P population to equilibrium
# ix: index of site
# psi: value of psi for site
# alpha: value of alpha for site
# pop: EL4P population
# PAR: Aq.PAR
# tMax: time to run aquatic dynamics
# tol: tolerance in lambda fluctuations until equilibrium is assumed
run2Eq_GEL4P <- function(ix, psi, pop, PAR, tMax = 800, tol = 0.1){
  pop = burnin_GEL4P(ix, psi, pop, PAR, tMax) # run aquatic populations through burnin
  pop = G2K_GEL4P(ix, psi, pop, PAR, tMax)
  TMAX = 100
  dx = checkDX_GEL4P(ix, psi, pop, PAR, TMAX) # check variance of lambda
  while(var(dx) > tol){ # run until lambda stabilizes
    # print(paste0("current variance of lambda: ",var(dx),", TMAX:",TMAX))
    pop = G2K_GEL4P(ix, psi, pop, PAR, TMAX)
    dx = checkDX_GEL4P(ix, psi, pop, PAR, TMAX)
    TMAX = TMAX + 100
  }
  return(pop)
}

# burnin_GEL4P: run a single population for given amount of time
burnin_GEL4P <- function(ix, psi, pop, PAR, tMax = 800){
  for(i in 1:tMax){
    pop = oneDay_GEL4P(ix = ix, psi = psi, pop = pop, PAR = PAR)
  }
  return(pop)
}

# G2K_GEL4P: run a single population for a given amount of time with simulated adult dynamics
G2K_GEL4P <- function(ix, psi, pop, PAR, tMax = 800){
  PAR$M[ix] = pop$lambda + 1
  for(i in 1:tMax){
    pop = oneDay_GEL4P(ix = ix, psi = psi, pop = pop, PAR = PAR)
    PAR$M[ix] = ((exp(-1/PAR$lifespan))*PAR$M[ix]) + pop$lambda # simulate adult population dynamics
  }
  return(pop)
}

# checkDX_GEL4P: run a single population for a given amount of time with simulated adult dynamics and output vector of lambda
checkDX_GEL4P <- function(ix, psi, pop, PAR, tMax = 800){
  M = pop$lambda # init adult population
  lambdaH = vector(mode="numeric",length=tMax+1)
  lambdaH[1] = pop$lambda
  for(i in 1:100){ # simulate adults for 100 days
    M = ((exp(-1/PAR$lifespan))*PAR$M[ix]) + pop$lambda # simulate adult population dynamics
  }
  PAR$M[ix] = M
  for(i in 1:tMax){ # run full simulation
    pop = oneDay_GEL4P(ix = ix, psi = psi, pop = pop, PAR = PAR)
    PAR$M[ix] = ((exp(-1/PAR$lifespan))*PAR$M[ix]) + pop$lambda # simulate adult population dynamics
    lambdaH[i+1] = pop$lambda
  }
  return(lambdaH)
}

# psi2K_cf: run linear regression of K on psi
psi2K_cf <- function(meshK, psiHat, plot = FALSE){
  psiInv = 1/psiHat
  psi2K = lm(psiInv~meshK+0)
  cf = coef(psi2K)
  # if(plot){
  #   plot(meshK, 1/psiHat)
  #   abline(psi2K)
  #   points(meshK, cf*meshK, col = "red", pch = 3)
  #   points(meshK, 1/K2psi(meshK,cf), col = "purple")
  # }
  if(plot){
    psi2K_plot(lmFit = psi2K,K = meshK,psi = psiHat)
  }
  return(cf)
}

#
psi2K_plot <- function(lmFit,K,psi){
  pCol = ggCol(n = 1,alpha = 0.8)
  plot(K,1/psi,type="p",pch=16,col=pCol)
  abline(lmFit)
  points(K, coef(lmFit)*K, col = "red", pch = 3)
  points(K, 1/K2psi(K = K,cf = coef(lmFit)),col = "purple")
}

# K2psi: convert K and regression coefficient into value of psi
K2psi <- function(K,cf){
  1/(cf*K)
}



#################################################################
# EL4P Visualizations
#################################################################


# plotPsi: visualize distributions of initial parameters
plotPsi <- function(PAR){

  with(PAR,{
    par(mfrow=c(2,2))
    col = ggCol(n = 4)

    hist(K,freq = F,main="K") # K
    lines(density(K,from=0),lwd=5,col=col[1])

    hist(alpha,freq=F,main=expression(alpha),cex.main=2,xlab=expression(alpha)) # alpha
    lines(density(alpha,from=0),lwd=5,col=col[2])

    psiDensity = density(psiInit,from=0) # psi
    hist(psiInit,freq=F,main=expression('initial (uniftted)'~psi),cex.main=2,ylim=c(0,max(psiDensity$y)),xlab=expression(psi))
    lines(psiDensity,lwd=5,col=col[3])

    plot(meshK,main="K mesh (sample in log-space)",type="p",pch=16,col=col[4],cex=2,xaxt="n",xlab="") # meshK
    grid()
    axis(1, at = seq(from=1,to=length(meshK),length.out = 10),labels = round(meshK[seq(from=1,to=length(meshK),length.out = 10)]))

    par(mfrow=c(1,1))
  })

}