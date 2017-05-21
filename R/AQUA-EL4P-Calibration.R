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


#' Calculate Lambda (Equilibrium Rate of Emergence)
#'
#' Calculate parameter \code{lambda} (daily adult female emergence over entire landscape) required to sustain \emph{Plasmodium falciparum} transmission at given value of \emph{R0} based on
#' classical Ross-MacDonald assumptions.
#'
#' @param R0 desired intensity of pathogen transmission
#' @param summary output from \code{\link{cohortBionomics}} which can be calculated from a run of MASH, or, less laboriously, from \code{\link{MBITES.basic}}
#' @param nH number of humans on landscape (ie; number of humans in the MICRO patch)
#' @param P list of parameters from \code{\link{MBITES.PAR}}
#' @return numeric value of lambda
#' @examples
#' getLambda(R0, summary, nH, P)
getLambda <- function(R0, summary, nH, P){
  with(c(P,summary),{
    P = exp(-EIP/lifespanC) # probability to survive EIP
    S = feedAllC # stability index
    # (R0 * nH * (S^2) * r) / (b * c * P)
    (R0*nH*r) / ((S^2)*b*c*P)
  })
}

#' Calculate Parameters for EL4P Aquatic Ecology Module
#'
#' Calculate parameters needed for EL4P fitting.
#'
#' @param nA number of aquatic habitats
#' @param nH number of humans
#' @param R0 desired pathogen transmission level at equilibrium
#' @param M vector of densities of mosquitoes at each aquatic habitat
#' @param aquaEq vector of probability of oviposition at each aquatic habitat (see \code{\link{aquaIx_equilibrium}})
#' @param par parameter list from \code{\link{MBITES.PAR}}
#' @param summary output from \code{\link{cohortBionomics}} which can be calculated from a run of MASH, or, less laboriously, from \code{\link{MBITES.basic}}.
#' @param p expected fraction of cohort that advances to next life stage (1/p is expected time spent in each stage E,L1,L2,L3,L4,P)
#' @param P density independent survival from egg to adult
#' @param a shape parameter for gamma distributed weights on K
#' @param b scale parameter for gamma distributed weights on K
#' @return named list of fitted parameters
#' * M: M
#' * aquaEq: aquaEq
#' * nA: nA
#' * nH: nH
#' * lifespan: mean length of mosquito lifespan
#' * G: mean total lifetime egg production of an adult female
#' * nu: mean batch size for one oviposition
#' * p: p
#' * P: P
#' * a: a
#' * b: b
#' * lambda: daily female emergence required to sustain transmission at input R0 (see \code{\link{getLambda}})
#' @md
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

# setupAquaPop_EL4P <- function(PAR, tol = 0.1, plot = FALSE){
#   with(PAR,{
#
#     # calculate initial parameter values
#     W = rgamma(n = nA,shape = a,scale = b)
#     K = (lambda*W) / sum(W)
#     PAR$K = K # attach K to PAR
#     pp = -log(P^((1-p)/5))
#     alpha = abs(rnorm(n = nA,mean = pp,sd = 0.004))
#     PAR$alpha = alpha # attach alpha to PAR
#     psi = alpha/K
#     PAR$psiInit = psi # attach psi to PAR
#
#     # EL4P populations
#     EL4P_pops = replicate(n = nA,expr = EL4P(),simplify = FALSE)
#
#     # fit psi
#     rng = range(K)
#     AquaPops = meshK_EL4P(lK = rng[1],uK = rng[2],EL4P_pops = EL4P_pops,PAR = PAR,plot = plot,tol = tol)
#     PAR$psiOptim = AquaPops$psiHat # attach psi fitted by optimize(...) to PAR
#     PAR$meshK = AquaPops$meshK # attach sampled meshK to PAR
#
#     if(plot){plotPsi(PAR = PAR)}
#     if(plot){par(mfrow=c(1,2))}
#     cf = psi2K_cf(AquaPops$meshK,AquaPops$psiHat,plot)
#     if(plot){
#       psi2K_plot(lmFit = lm(1/K2psi(AquaPops$meshK,cf)~AquaPops$meshK+0),K = AquaPops$meshK,psi = K2psi(AquaPops$meshK,cf))
#       par(mfrow=c(1,1))
#     }
#     psiHat = K2psi(AquaPops$meshK,cf)
#     PAR$psiHat = psiHat # attach psi fitted via linear regression to PAR
#
#     # run all EL4P aquatic populations to equilibrium values
#     EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = psiHat,pop = AquaPops$EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
#                                    mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)
#     return(
#       list(EL4P_pops=EL4P_pops,
#            PAR=PAR
#       )
#     )
#   })
# }
#
# setupAquaPop_EL4PsamplePoints <- function(PAR, gridN, tol = 0.1, plot = FALSE){
#   with(PAR,{
#
#     # calculate initial parameter values
#     W = rgamma(n = nA,shape = a,scale = b)
#     K = (lambda*W) / sum(W)
#     pp = -log(P^((1-p)/5))
#     alpha = abs(rnorm(n = nA,mean = pp,sd = 0.004))
#     psi = alpha/K
#
#     PAR$alpha = alpha # attach alpha to PAR
#     PAR$K = K # attach K to PAR
#     PAR$psiInit = psi # attach psi to PAR
#
#     if(plot){
#       plotPsi(PAR = PAR)
#       par(mfrow=c(1,3))
#       invisible(psi2K_cf(meshK = K,psiHat = psi,plot = T,main="Prior to Fitting EL4P"))
#     }
#
#     # generate aquatic populations
#     EL4P_pops = replicate(n = nA,expr = EL4P(),simplify = FALSE)
#
#     # fit psi to mesh of K values
#
#     # sample K on mesh in log space; transform to normal space
#     rng = range(K)
#     AquaPops = meshK_EL4P(lK = rng[1],uK = rng[2],EL4P_pops = EL4P_pops,PAR = PAR,plot = plot,tol = tol)
#     PAR$psiOptim = AquaPops$psiHat # attach psi fitted by optimize(...) to PAR
#     PAR$meshK = AquaPops$meshK # attach sampled meshK to PAR
#
#     # regression of psi inverse on K
#     cf = psi2K_cf(meshK = AquaPops$meshK,psiHat = PAR$psiOptim,plot = plot, main = "After Optimization")
#     if(plot){
#       psi2K_plot(lmFit = lm(1/K2psi(AquaPops$meshK,cf)~K+0),K = AquaPops$meshK,psi = K2psi(AquaPops$meshK,cf), main = "After Regression")
#       par(mfrow=c(1,1))
#     }
#
#     # attach psi fitted via linear regression to PAR for final fitted psi
#     PAR$psiHat = K2psi(K = K,cf = cf)
#
#     # run all EL4P aquatic populations to equilibrium values on final fitted psi values
#     EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = PAR$psiHat,pop = AquaPops$EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
#                                    mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)
#
#     return(
#       list(EL4P_pops = EL4P_pops, PAR = PAR)
#     )
#
#   })
# }

#' Fit EL4P Aquatic Ecology Model on Exact Landscape
#'
#' Fit the EL4P Aquatic Ecology module on the exact LANDSCAPE to match a desired level of daily emergence at equilibrium.
#' This means that the site-specific density dependent mortality parameter \code{psi} will be fit
#' to K for each aquatic habitat on the landscape, and then each site will be run to equilibrium. Compare with \code{\link{setupAquaPop_EL4PsamplePoints}} which will
#' fit \code{psi} based on a sampling grid of K.
#'
#' @param PAR named list of parameters calculated from \code{\link{makePAR_EL4P}}
#' @param tol target minimum variance in lambda for aquatic populations equilibrium
#' @param tol tolerance of fluctuations in variance of lambda until convergence to equilibrium is assumed
#' @param plot produce diagnostic plots?
#' @return named list of fitted parameters
#' * EL4P_pops: fitted EL4P aquatic populations
#' * PAR: input PAR with fitted values of \code{psi} and \code{alpha} appended
#' @md
#' @examples
#' setupAquaPop_EL4Pexact(PAR, tol = 0.1, plot = FALSE)
setupAquaPop_EL4Pexact <- function(PAR, tol = 0.1, plot = FALSE){
  with(PAR,{

    # calculate initial parameter values
    W = rgamma(n = nA,shape = a,scale = b)
    K = (lambda*W) / sum(W)
    pp = -log(P^((1-p)/5))
    alpha = abs(rnorm(n = nA,mean = pp,sd = 0.004))
    psi = alpha/K

    PAR$alpha = alpha # attach alpha to PAR
    PAR$K = K # attach K to PAR
    PAR$psiInit = psi # attach psi to PAR

    if(plot){
      plotPsi(PAR = PAR)
      par(mfrow=c(1,3))
      invisible(psi2K_cf(meshK = K,psiHat = psi,plot = TRUE,main="Prior to Fitting EL4P"))
    }

    # generate aquatic populations
    EL4P_pops = replicate(n = nA,expr = EL4P(),simplify = FALSE)

    # fit psi to exact landscape
    psiOptim = vector(mode="numeric",length=length(EL4P_pops))

    # fit EL4P: set values of psi so lambda = K at (p,G)
    for(ix in 1:length(EL4P_pops)){
      print(paste0("fitting psi for site ix: ",ix, ", target K: ",K[ix]))
      psiMin = optimize(f = psiFit,interval = c(0,10), ix=ix, ixEL4P=EL4P_pops[[ix]], ixKmesh=K[ix], PAR=PAR)
      psiOptim[ix] = psiMin$minimum
    }
    PAR$psiOptim = psiOptim # attach psi fitted by optimize(...) to PAR

    # run all EL4P aquatic populations to equilibrium values
    print(paste0("run aquatic populations to equilibrium values"))
    EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = psiOptim,pop = EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
                                   mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)

    # regression of psi inverse on K
    cf = psi2K_cf(meshK = K,psiHat = PAR$psiOptim,plot = plot, main = "After Optimization")
    if(plot){
      psi2K_plot(lmFit = lm(1/K2psi(K,cf)~K+0),K = K,psi = K2psi(K,cf), main = "After Regression")
      par(mfrow=c(1,1))
    }

    # attach psi fitted via linear regression to PAR for final fitted psi
    PAR$psiHat = K2psi(K = K,cf = cf)

    # run all EL4P aquatic populations to equilibrium values on final fitted psi values
    EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = PAR$psiHat,pop = EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
                                   mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)

    return(
      list(EL4P_pops = EL4P_pops, PAR = PAR)
    )
  })
}


#' Fit EL4P Aquatic Ecology Model on Sample Grid
#'
#' Fit the EL4P Aquatic Ecology module on a sampling grid of K values to match a desired level of daily emergence at equilibrium.
#' This means that the site-specific density dependent mortality parameter \code{psi} will be fit
#' to K based on a sampling grid of values for K in log-space. If \code{plot = TRUE}, the linear regression of \code{psi} against logged values of K
#' should show exact linear dependence, indicating fitted \code{psi} will produce desired level of lambda at equilibrium. Compare with \code{\link{setupAquaPop_EL4Pexact}} which will
#' fit \code{psi} based on an exact LANDSCAPE. This will also return the coefficients of a linear regression of K on psi (see \code{\link{psi2K_cf}}) to give the functional relationship between K and psi.
#'
#' @param PAR named list of parameters calculated from \code{\link{makePAR_EL4P}}
#' @param gridN number of points to sample in K
#' @param tol target minimum variance in lambda for aquatic populations equilibrium
#' @param plot produce diagnostic plots?
#' @return named list of fitted parameters
#' * EL4P_pops: fitted EL4P aquatic populations
#' * PAR: input PAR with fitted values of \code{psi} and \code{alpha} appended
#' @md
#' @examples
#' setupAquaPop_EL4PsamplePoints(PAR, gridN, tol = 0.1, plot = FALSE)
setupAquaPop_EL4PsamplePoints <- function(PAR, gridN, tol = 0.1, plot = FALSE){
  with(PAR,{

    # calculate initial parameter values
    W = rgamma(n = gridN,shape = a,scale = b)
    K = (lambda*W) / sum(W)
    pp = -log(P^((1-p)/5))
    alpha = abs(rnorm(n = gridN,mean = pp,sd = 0.004))
    psi = alpha/K

    PAR$alpha = alpha # attach alpha to PAR
    PAR$K = K # attach K to PAR
    PAR$psiInit = psi # attach psi to PAR

    if(plot){
      plotPsi(PAR = PAR)
      par(mfrow=c(1,3))
      invisible(psi2K_cf(meshK = K,psiHat = psi,plot = TRUE,main="Prior to Fitting EL4P"))
    }

    # generate aquatic populations
    EL4P_pops = replicate(n = gridN,expr = EL4P(),simplify = FALSE)

    # fit psi to mesh of K values

    # sample K on mesh in log space; transform to normal space
    rng = range(K)
    AquaPops = meshK_EL4P(lK = rng[1],uK = rng[2],EL4P_pops = EL4P_pops, PAR = PAR, tol = tol)
    PAR$psiOptim = AquaPops$psiHat # attach psi fitted by optimize(...) to PAR
    PAR$meshK = AquaPops$meshK # attach sampled meshK to PAR

    # regression of psi inverse on K
    cf = psi2K_cf(meshK = AquaPops$meshK,psiHat = PAR$psiOptim,plot = plot, main = "After Optimization")
    PAR$cf = cf # attach coefficients to convert K to psi to PAR
    if(plot){
      psi2K_plot(lmFit = lm(1/K2psi(AquaPops$meshK,cf)~K+0),K = AquaPops$meshK,psi = K2psi(AquaPops$meshK,cf), main = "After Regression")
      par(mfrow=c(1,1))
    }

    # attach psi fitted via linear regression to PAR for final fitted psi
    PAR$psiHat = K2psi(K = K,cf = cf)

    # run all EL4P aquatic populations to equilibrium values on final fitted psi values
    EL4P_pops = parallel::mcmapply(FUN = run2Eq_GEL4P,ix = 1:length(EL4P_pops),psi = PAR$psiHat,pop = AquaPops$EL4P_pops, MoreArgs = list(tol = tol, PAR = PAR),
                                   mc.cores = parallel::detectCores()-2L,SIMPLIFY = FALSE)

    return(
      list(EL4P_pops = EL4P_pops, PAR = PAR)
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

#' Fit Psi on Mesh of K Values
#'
#' Generate a mesh of K values in log space and fit psi to mesh via numerical optization, then run input EL4P populations to equilibrium.
#' This sets values of psi so that lambda \eqn{\lambda=K} at given parameter values.
#'
#' @param lK lower bound of K mesh
#' @param uK upper bound of K mesh
#' @param alpha vector of density independent mortality
#' @param EL4P_pops list of aquatic populations
#' @param PAR list of parameters from \code{\link{makePAR_EL4P}}
#' @param tMax maximum time to to run aquatic populations to equilibrium
#' @param tol target minimum variance in lambda for aquatic populations equilibrium
#' @return named list of output
#' * EL4P_pops: aquatic populations run to equilibrium values
#' * psiHat: fitted values of psi
#' * meshK: sampled values of K
#' @md
#' @examples
#' meshK_EL4P(lK, uK, EL4P_pops, PAR, tMax = 500, tol = 0.1)
meshK_EL4P <- function(lK, uK, EL4P_pops, PAR, tMax = 500, tol = 0.1){

  # sample K on mesh in log space; transform to normal space
  meshK = exp(seq(log(lK),log(uK),length.out=length(EL4P_pops)))
  psiHat = vector(mode="numeric",length=length(EL4P_pops))

  # fit EL4P: set values of psi so lambda = K at (p,G)
  for(ix in 1:length(EL4P_pops)){
    print(paste0("fitting psi for site ix: ",ix, ", target K: ",meshK[ix]))
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

#' Objective Function for psi Optimization
#'
#' Given a input psi, \code{x}, run a single aquatic population and output objective function; squared error of lambda around the given value of K \eqn{\left ( \lambda-K \right )^{2}}.
#' This will typically be called by \code{\link{meshK_EL4P}} or \code{optimize}; generally it is the objective function that will be passed to \code{optimize(...)}
#'
#' @param x value of psi from \code{optimize(...)}
#' @param ix index of this aquatic population
#' @param ixEL4P EL4P aquatic population object corresponding to \code{ix} (see \code{link\{EL4P})
#' @param ixKmesh the value of K for this aquatic population
#' @param ixAlpha the value of alpha for this aquatic population
#' @param PAR aquatic parameters passed from \code{\link{makePAR_EL4P}}
#' @param tMax 150; time to run \code{\link{runOne_GEL4P}} before calculating value of objective function
#' @return value of objective function at \code{x}
psiFit <- function(x, ix, ixEL4P, ixKmesh, PAR, tMax = 150){
  psi = abs(x)
  EL4P = runOne_GEL4P(ix = ix, psi = psi, pop = ixEL4P, PAR = PAR, tMax = tMax)
  return((EL4P$lambda - ixKmesh)^2) # return squared error of empirical lambda around value of meshK
}

#' Run a Single Aquatic Population for \code{psiFit}
#'
#' Run a single aquatic population's daily difference equation dynamics for \code{\link{psiFit}}
#'
#' @param ix index of this aquatic population
#' @param psi value of psi
#' @param pop EL4P aquatic population object (see \code{link\{EL4P})
#' @param PAR aquatic parameters passed from \code{\link{makePAR_EL4P}}
#' @param tMax 150; time to run the dynamic model
#' @return the aquatic population
runOne_GEL4P <- function(ix, psi, pop, PAR, tMax = 150){
  for(i in 1:30){
    pop = oneDay_GEL4P(ix = ix,psi = psi,pop = pop,PAR = PAR)
  }
  # lambdaH = vector(mode="numeric",length = tMax+1)
  # lambdaH[1] = pop$lambda # record values of lambda
  for(i in 1:tMax){
    pop = oneDay_GEL4P(ix = ix,psi = psi,pop = pop,PAR = PAR)
    PAR$M[ix] = ((exp(-1/PAR$lifespan))*PAR$M[ix]) + pop$lambda # simulate adult population dynamics
    # lambdaH[i+1] = pop$lambda
  }
  # plot(lambdaH, type = "l", main = paste0("lambdaH for site: ",ix), col = "red", lwd = 2)
  return(pop)
}

#' Run a Single Daily Time-step of Aquatic Population Dynamics
#'
#' Run a single aquatic population's daily difference equation dynamics for a single time step. Eggs are deposited assuming adult population at equilibrium
#' \eqn{M[ix]*aquaEq[ix]*\left ( \frac{G}{lifespan} \right )}, see \code{link{makePAR_EL4P}} for these parameters.
#'
#' @param ix index of this aquatic population
#' @param psi value of psi
#' @param pop EL4P aquatic population object (see \code{link\{EL4P})
#' @param PAR aquatic parameters passed from \code{\link{makePAR_EL4P}}
#' @return the aquatic population state variables after one day
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

#' Run a Single Aquatic Population to Equilibrium
#'
#' Run a single aquatic population to equilibrium (where variance in emergence, lambda is less than \col{tol}).
#' This function first runs the population through a burnin period (see \code{\link{burnin_GEL4P}}) then runs aquatic dynamics with simulated adult
#' dynamics from derived Ross-MacDonald parameters (see \code{\link{G2K_GEL4P}}). Then the dynamics are run while variance in lambda is above \code{tol}
#' (see \code{\link{checkDX_GEL4P}).
#'
#' @param ix index of this aquatic population
#' @param psi value of psi
#' @param pop EL4P aquatic population object (see \code{link\{EL4P})
#' @param PAR aquatic parameters passed from \code{\link{makePAR_EL4P}}
#' @param tMax 800; time to run the model to equilibrium
#' @param tol target minimum variance in lambda for aquatic populations equilibrium
#' @return the aquatic population
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

#' Run a Single Aquatic Population through Burnin
#'
#' Run a single aquatic population through a burnin period (until initial oscillations dampen).
#'
#' @param ix index of this aquatic population
#' @param psi value of psi
#' @param pop EL4P aquatic population object (see \code{link\{EL4P})
#' @param PAR aquatic parameters passed from \code{\link{makePAR_EL4P}}
#' @param tMax 800; time to run the model through burnin period.
#' @return the aquatic population
burnin_GEL4P <- function(ix, psi, pop, PAR, tMax = 800){
  for(i in 1:tMax){
    pop = oneDay_GEL4P(ix = ix, psi = psi, pop = pop, PAR = PAR)
  }
  return(pop)
}

#' Run a Single Aquatic Population with Simulated Adult Dynamics
#'
#' Run a single aquatic population's dynamics using simulated adult dynamics and egg laying from derived Ross-MacDonald parameters (see \code{\link{cohortBionomics}}).
#' Adult dynamics follow the simple equation \eqn{\left ( e^{-1/lifespan} \right *M[ix]) + \lambda} where lambda parameter is from the aquatic population being simulated.
#'
#' @param ix index of this aquatic population
#' @param psi value of psi
#' @param pop EL4P aquatic population object (see \code{link\{EL4P})
#' @param PAR aquatic parameters passed from \code{\link{makePAR_EL4P}}
#' @param tMax 800; time to run the model through burnin period.
#' @return the aquatic population
G2K_GEL4P <- function(ix, psi, pop, PAR, tMax = 800){
  PAR$M[ix] = pop$lambda + 1
  for(i in 1:tMax){
    pop = oneDay_GEL4P(ix = ix, psi = psi, pop = pop, PAR = PAR)
    PAR$M[ix] = ((exp(-1/PAR$lifespan))*PAR$M[ix]) + pop$lambda # simulate adult population dynamics
  }
  return(pop)
}

#' Run a Single Aquatic Population with Simulated Adult Dynamics and Output Daily Lambda
#'
#' Run a single aquatic population's dynamics using simulated adult dynamics and egg laying from derived Ross-MacDonald parameters (see \code{\link{cohortBionomics}}).
#' Adult dynamics follow the simple equation \eqn{( e^{-1/lifespan} *M[ix]) + \lambda} where lambda parameter is from the aquatic population being simulated.
#' This function will output the record of daily lambda emergence, used to check variance in \code{\link{run2Eq_GEL4P}}.
#'
#' @param ix index of this aquatic population
#' @param psi value of psi
#' @param pop EL4P aquatic population object (see \code{link\{EL4P})
#' @param PAR aquatic parameters passed from \code{\link{makePAR_EL4P}}
#' @param tMax 800; time to run the model through burnin period.
#' @return vector of lambda values over time \code{tMax}
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

#' Regress K on Psi and Return Coefficients
#'
#' Run a linear regression of K on iverse of psi and extract coefficients; optionally plot the regression. Psi and K should follow a linear relationship, and
#' the coefficients of the linear regression will give the parameter of the 1-dimensional response surface that relates the two parameters.
#'
#' @param meshK vector of sampled K values
#' @param psiHat vector of psi calculated from \link{\code{psiFit}} which provides reasonable initial starting values
#' @param plot plot the regression line?
#' @return coefficients of regression
psi2K_cf <- function(meshK, psiHat, plot = FALSE, ...){
  psiInv = 1/psiHat
  psi2K = lm(psiInv~meshK+0)
  cf = coef(psi2K)
  # plot unfitted values of psi
  if(plot){
    psi2K_plot(lmFit = psi2K,K = meshK,psi = psiHat, ...)
  }
  return(cf)
}

#' Plot Regression of K on Psi
#'
#' Write docs.
#'
#' @param lmFit NULL
#' @param K NULL
#' @param psi NULL
#' @param main NULL
#' @return nothing
psi2K_plot <- function(lmFit, K, psi, main = NULL){
  pCol = ggCol(n = 1,alpha = 0.8)
  plot(K,1/psi,type="p",pch=16,cex=1.15,col=pCol, ylab = expression(paste(1/psi," (density-dependent mortality)")), xlab = "K (carrying capacities)",main = main)
  legend(x = "topleft",legend = expression(paste("Regression of ",1/psi," on K")),bty = "n")
  grid()
  abline(lmFit)
  points(K, coef(lmFit)*K, col = "red", pch = 3)
  points(K, 1/K2psi(K = K,cf = coef(lmFit)),col = "purple")
}

#' Convert K and Regression Coefficient into Psi
#'
#' Write docs.
#'
#' @param K NULL
#' @param cf NULL
#' @return psi
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

    if(exists(x = "meshK")){
      plot(meshK,main="K mesh (sample in log-space)",type="p",pch=16,col=col[4],cex=2,xaxt="n",xlab="") # meshK
      grid()
      axis(1, at = seq(from=1,to=length(meshK),length.out = 10),labels = round(meshK[seq(from=1,to=length(meshK),length.out = 10)]))
    } else {
      plot(K,main="K (values from landscape)",type="p",pch=16,col=col[4],cex=2,xlab="Aquatic Habitats") # meshK
      grid()
    }

    par(mfrow=c(1,1))
  })

}
