

#######################################################################

# Optimise investment for combined immune strategies

#######################################################################
# Optimisation
#######################################################################
# Trapezium rule
trap.rule <- function(x,dt){
  return((sum(x)-0.5*(x[1]+x[length(x)]))*dt)
}

# Function to optimise
# c.f has three components: (c,nu_L,nu_P)
f.optim.fn.combine <- function (c.f,parms.f){
  if (1-c.f[2]-c.f[3]<0){
    # Check nuT not negative
    return(-10^6)
  } else {
    # Calculate trajectory
    traj <- f.traj.combine(c("c"=c.f[1],"nuL"=c.f[2],"nuP"=c.f[3],parms.f),t.save)
    if (any(traj$C<=0)){
      # If host dies, artificially reduce mean condition by 10^6
      return(trap.rule(traj$C,t.save[2]-t.save[1])/max(t.save)-10^6)
    } else{
      return(trap.rule(traj$C,t.save[2]-t.save[1])/max(t.save))
    }
  }
}

# Optimisation function
f.optim.combine <- function(parms.f,c0=c(0.1,0.33,0.33)){
  opt <- optim(c0,f.optim.fn.combine,
               parms.f=parms.f,
               method="L-BFGS-B",
               lower=c(0,0,0),
               upper=c(Inf,1,1),
               control=list("fnscale"=-1))
  opt.traj <- f.traj.combine(c("c"=opt[["par"]][1],"nuL"=opt[["par"]][2],"nuP"=opt[["par"]][3],parms.f),t.save)
  return(data.frame("c"=opt[["par"]][1],"nuL"=opt[["par"]][2],"nuP"=opt[["par"]][3],as.list(parms.f),
                    "C.Mean"=opt[["value"]],
                    "M"=opt.traj$M[length(t.save)],
                    "L"=opt.traj$L[length(t.save)],
                    "R"=opt.traj$R[length(t.save)],
                    "I"=opt.traj$I[length(t.save)],
                    "convergence"=opt["convergence"]))
}

######################################################################
# Vary SR, SL, hL, hI
######################################################################

# Number of points
n.SR <- length(SRs)
n.SL <- length(SLs)
n.g <- length(gs)
n.hL <- length(hLs)
n.hI <- length(hIs)

# Optimise
c.opt <- data.frame()

for (j.SR in 1:n.SR){
  for (j.SL in 1:n.SL){
    for (j.g in 1:n.g){
      for (j.hL in 1:n.hL){
        for (j.hI in 1:n.hI){
          c.opt <- rbind(c.opt,
                         f.optim.combine(parms.f=c("SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                                   "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=hIs[j.hI],"hI1"=hIs[j.hI]),
                                         c0=c(0.1,0.33,0.33)
                                         )
                         )
        }
      }
    }
  }
}
