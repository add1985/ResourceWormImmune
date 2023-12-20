
# Optimise investment for each strategy

#######################################################################
# Optimisation
#######################################################################
# Trapezium rule
trap.rule <- function(x,dt){
  return((sum(x)-0.5*(x[1]+x[length(x)]))*dt)
}

# Function that needs optimising
# Optimise over c.f
f.optim.fn <- function (c.f,strat,parms.f){
  # Calculate trajectory
  traj <- f.traj(strat,c("c"=c.f,parms.f),t.save)
  # If host dies, return zero
  # Otherwise, calculate mean condition
  if (any(traj$C<=0)){
    # If host dies, artificially reduce mean condition by 10^6
    return(trap.rule(traj$C,t.save[2]-t.save[1])/max(t.save)-10^6)
  } else{
    return(trap.rule(traj$C,t.save[2]-t.save[1])/max(t.save))
  }
}

# Wrapper to perform optimisation
# parms.f contains parameters the need to be specified
f.optim <- function(strat,parms.f,c0=0.1){
  opt <- optim(c0,f.optim.fn,
               strat=strat,
               parms.f=parms.f,
               method="L-BFGS-B",
               lower=0,
               control=list("fnscale"=-1))
  opt.traj <- f.traj(strat,c("c"=opt[["par"]],parms.f),t.save)
  return(data.frame("strat"=as.character(strat),"c"=opt[["par"]],as.list(parms.f),
                    "C.Mean"=opt[["value"]],
                    "R"=opt.traj$R[length(t.save)],
                    "M"=opt.traj$M[length(t.save)],
                    "L"=opt.traj$L[length(t.save)],
                    "I"=opt.traj$I[length(t.save)],
                    "C"=opt.traj$C[length(t.save)],
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
    for(j.g in 1:n.g){
      for (j.hL in 1:n.hL){
        c.opt <- rbind(c.opt,
                     # Anorexia and tolerance don't have immunopathology
                     bind_rows(lapply(c(1,4),f.optim,
                                      parms.f=c("SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                                "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=0,"hI1"=0),
                                      c0=0.5)
                                 ),
                       # No strategy
                       # No need to optimise; calculate directly using f.optim.fn with anorexia, c=0
                       data.frame("strat"="0","c"=0,
                                  as.list(c("SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                            "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=0,"hI1"=0)),
                                  "C.Mean"=f.optim.fn(0,1,
                                                      c("SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                                        "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=0,"hI1"=0)
                                                      ),
                                  f.traj(1,
                                         c("c"=0,"SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                           "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=0,"hI1"=0),
                                         t.save
                                         )[length(t.save),c("R","M","L","I","C")],
                                  "convergence"=0
                                  ),
                       # Starvation
                       # No need to optimise; calculate directly using f.optim.fn with anorexia, c very large
                       data.frame("strat"="1.5","c"=0,
                                  as.list(c("SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                            "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=0,"hI1"=0)),
                                  "C.Mean"=f.optim.fn(10^7,1,
                                                      c("SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                                        "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=0,"hI1"=0)
                                  ),
                                  f.traj(1,
                                         c("c"=10^7,"SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                           "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=0,"hI1"=0),
                                         t.save
                                  )[length(t.save),c("R","M","L","I","C")],
                                  "convergence"=0
                       )
        )
        for (j.hI in 1:n.hI){
          c.opt <- rbind(c.opt,
                         # Immunopathology for prevention and clearance
                         bind_rows(lapply(c(2,3),f.optim,
                                          parms.f=c("SR"=SRs[j.SR],"SL"=SLs[j.SL],"g"=gs[j.g],
                                                    "hL"=hLs[j.hL],"hM"=hMs[j.hL],"hI0"=hIs[j.hI],"hI1"=hIs[j.hI]),
                                          c0=0.5)
                                   )
                         )
        }
      }
    }
  }
}
