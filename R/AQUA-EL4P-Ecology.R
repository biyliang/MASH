#################################################################
#
#   MASH/MBITES
#   Aquatic Ecology
#   EL4P
#   R Version
#   Sean Wu
#   March 22, 2017
#
#################################################################

# EL4P: object to store aquatic forms
EL4P <- function(eggs = 0, L1 = 0, L2 = 0, L3 = 0, L4 = 0, P = 0, lambda = 0){
  list(
    eggs   = eggs,   # eggs
    L1     = L1,     # larvae 1
    L2     = L2,     # larvae 2
    L3     = L3,     # larvae 3
    L4     = L4,     # larvae 4
    P      = P,      # pupae
    lambda = lambda  # emerging adults
  )
}

# EL4Pt: daily difference equations
# D: total larval density for entire pool
EL4Pt <- function(psi, alpha, p, pop){
  L1o = pop$L1; L2o=pop$L2; L3o=pop$L3; L4o=pop$L4
  D   = sum(L1o+L2o+L3o+L4o)
  s1  = exp(-alpha)
  s2  = exp(-(alpha+ psi*D))
  pop$lambda = s1*pop$P
  pop$P  = s2*p*L4o
  pop$L4 = s2*(p*L3o + (1-p)*L4o)
  pop$L3 = s2*(p*L2o + (1-p)*L3o)
  pop$L2 = s2*(p*L1o + (1-p)*L2o)
  pop$L1 = pop$eggs + s2*(1-p)*L1o
  pop$eggs = 0
  return(pop)
}

# oneSiteEL4Pt: run daily difference equations for a single aquatic habitat
oneSiteEL4Pt <- function(ix, tNow){

  addEggs2EL4P(ix, tNow) # move eggs from EggQ to EL4P (EL4P-Ecology)

  with(LANDSCAPE$aquaSites[[ix]],{ # daily difference equations in EL4P
    LANDSCAPE$aquaSites[[ix]]$EL4P <<- EL4Pt(psi = psi,alpha = alpha,p = p,pop = EL4P)
  })

  if(LANDSCAPE$aquaSites[[ix]]$EL4P$lambda > 0){ # if emerge, add to ImagoQ
    lambdaEmerge = rpois(n = 1,lambda = LANDSCAPE$aquaSites[[ix]]$EL4P$lambda)
    if(lambdaEmerge > 0){
      addAdults2Q(lambda = lambdaEmerge,tm = tNow,ix = ix,dam = 0,sire = 0) # aquaticEcology.R
    }
  }

}

# oneDay_EL4PL run EL4P aquatic ecology module for one day
oneDay_EL4P <- function(tNow){

  # log EggQ
  if(EggQ_TRACK){
    trackEggQ(con = .GlobalEnv$EggQCon)
  }

  for(ixA in 1:LANDSCAPE$nA){ # iterate over aquatic habitats
    oneSiteEL4Pt(ix = ixA, tNow = tNow)
  }

}


#################################################################
# Management of interaction between EggQ, EL4P, ImagoQ
#################################################################

# addEggs2EL4P: manage eggQ to EL4P transition for site ix
addEggs2EL4P <- function(ix, tNow){

  with(LANDSCAPE$aquaSites[[ix]],{

    ixQ = which(sapply(EggQ,function(x)({(x$tm < tNow) & (x$tm > 0) & (x$N > 0)}))) # filled EggQ

    for(ixEgg in ixQ){ # ixEgg iterates through ixQ
      addBatch2EL4P(ixQ = ixEgg,ix = ix)
      zeroBatch(ixQ = ixEgg,ix = ix)
    }

  })

}

# addBatch2EL4P
addBatch2EL4P <- function(ixQ, ix){
  LANDSCAPE$aquaSites[[ix]]$EL4P$eggs <<- LANDSCAPE$aquaSites[[ix]]$EL4P$eggs + LANDSCAPE$aquaSites[[ix]]$EggQ[[ixQ]]$N
}
