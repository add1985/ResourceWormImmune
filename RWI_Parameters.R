
#######################################################################
# Parameters
#######################################################################

# Time points
tmax <- 90
t <- seq(0,tmax,0.05)
t.save <- seq(0,tmax,0.5)

# Parameters which take a range
cs <- seq(0.0,1,length.out=201) #
cs.i.scale <- 1
cs.i <- cs*cs.i.scale
SRs <- c(seq(0.1,5,0.05)) #seq(0.5,5,0.5) #
SLs <- c(0.5,2) #c(seq(0.,4,0.1)) # 0.5 #
hLs <- c(0.2,0.4,0.6)
hMs <- c(0.6,0.4,0.2)
gs <- c(0.1,0.5)
hIs <- c(0.1,0.25,0.5)

# Default values for other parameters
parms.default <- c(
  "SR" = 4,
  "r" = 1, # 1. 5
  "c" = 0.1,
  "p" = 0.01, # IC for immune response
  "q" = 0.1, # 0.1, 0.5, 1
  "s" = 1,
  "v" = 0.5, # 0.5, 2.5
  "l" = 0.1,
  "SL" = 0.5,
  "g" = 0.1, # 0.5
  "dL" = 0.1,
  "dM" = 0.02,
  "a" = 2, # 2, 4
  "b" = 1,
  "w" = 1, # 1, 2
  "hL" = 0.4,
  "hM" = 0.4,
  "hI0" = 0.25,
  "hI1" = 0.25,
  "k0" = 1, # 1, 5
  "k1" = 0.5, # 0.5, 2.5,
  "nuL" = 0.2,
  "nuP" = 0.2
)

# Calculate time to host death in absence of resources and parasites
# No immune response, initially "well-resourced"
f.check_w <- function(SR,w,a,b,r){
  f.C_starve <- function(t,y,parms){return(list(a*SR*exp(-r*t)/(1+b*y)-w))}
  traj.C_starve <- data.frame(
    ode(y = c(y=(a*SR/w-1)/b),
        times = seq(0,90,0.1),
        func = f.C_starve,
        method = "lsode")
    )
  return(traj.C_starve$time[min(which(traj.C_starve$y<=0))])
}
SR_WellResourced <- 5
print(f.check_w(SR_WellResourced,parms.default["w"],parms.default["a"],parms.default["b"],parms.default["r"]))

# Must be positive for initially alive host
print(parms.default["a"]*SRs-parms.default["w"])
SRs <- SRs[parms.default["a"]*SRs-parms.default["w"]>0]

#######################################################################
# Calculate trajectory
#######################################################################

# ODe system
# This includes strategy in parms
f.sys <- function(times,state,parms){
  # Variables
  R <- state["R"]
  L <- state["L"]
  M <- state["M"]
  I <- state["I"]
  C <- state["C"]
  
  # Strategy
  strat <- parms["strat"]
  # Parameters
  SR <- parms["SR"]
  r <- parms["r"]
  s <- parms["s"]
  v <- parms["v"]
  SL <- parms["SL"]
  g <- parms["g"]
  dL <- parms["dL"]
  dM <- parms["dM"]
  l <- parms["l"]
  p <- parms["p"]
  q <- parms["q"]
  a <- parms["a"]
  b <- parms["b"]
  w <- parms["w"]
  hL <- parms["hL"]
  hM <- parms["hM"]
  hI0 <- parms["hI0"]
  hI1 <- parms["hI1"]
  cc <- parms["c"]
  k <- parms["k0"]*cc/(1+parms["k1"]*cc)
  
  # ODE
  if (strat%in%c(2,3,4)){
    FIP <- p+q*(L+s*M)*I
    DR <- SR - r*R - cc*FIP*R/(1+v*FIP*R) 
    DL <- SL - (g+dL+(strat==2)*k*I)*L
    DM <- g*L-(dM+(strat==3)*k*I)*M
    DI <- FIP*R/(1+v*FIP*R) - l*I
    # Let condition be negative for the purposes of optimisation
    DC <- (a*r*R/(1+b*C)-w-(hL*L+hM*M)/(1+(strat==4)*k*I)-(hI0*k+hI1*k^2)*I)
  } else if (strat==1) {
    # No immune response if anorexia (strat=1)
    DR <- SR/(1+ cc*(L+s*M)) - r*R
    DL <- SL/(1+ cc*(L+s*M)) - (g+dL)*L
    DM <- g*L-dM*M
    DI <- 0
    # Let condition take negative values for the purposes of optimisation
    # Set negative values to host death post hoc
    DC <- (a*r*R/(1+b*C)-w-(hL*L+hM*M))
  }
  
  return(list(c(DR,DL,DM,DI,DC)))
}

# Wrapper to calculate trajectory and return specified points
# parms.f contains parameteres we wish to specify
f.traj <- function(strat,parms.f,t.spec){
  # Named vector of parameters
  # Any parameters not specified in parms.f are given by parms.default
  # Include strategy for passing to f.sys
  parms <- c(parms.f,
         parms.default[!names(parms.default)%in%names(parms.f)],
         "strat" = strat
         )
  # Check there is no immunopathology for tolerance or anorexia
  if (strat %in% c(1,4)){
    parms["hI0"] <- 0
    parms["hI1"] <- 0 
  }
    # Parasite-free initial conditions
  SR <- parms["SR"]
  r <- parms["r"]
  v <- parms["v"]
  p <- parms["p"]
  l <- parms["l"]
  a <- parms["a"]
  b <- parms["b"]
  w <- parms["w"]
  hI0 <- parms["hI0"]
  hI1 <- parms["hI1"]
  cc <- parms["c"]
  k <- parms["k0"]*cc/(1+parms["k1"]*cc)
  if (strat%in%c(2,3,4) & p==0){
    # Immune-mediated response
    R0 <- SR/r
    I0 <- p*r*R0/l
  }else{
    # No immune response for anorexia
    R0 <- SR/r
    I0 <- 0
  }
  C0 <- (a*r*R0/(w+(hI0*k+hI1*k^2)*I0)-1)/b
  state.init <- c(R0,0,0,I0,C0)
  names(state.init) <- c("R","L","M","I","C")
  # Calculate trajectory
  traj <- data.frame(ode(y = state.init,
                         times = t,
                         func = f.sys,
                         parms = parms,
                         method = "lsode"))
  
  # Return
  return(traj[traj$time %in% t.spec,])
}
# 
ggplot(pivot_longer(f.traj(1,c("SR"=5,"SL"=0.5,"c"=0.5),t.save),-time,values_to="Value",names_to="Variable"),
      aes(x=time,y=Value,color
          =Variable))+
   geom_line()

#######################################################################
# Calculate trajectory, combining strategies
# No anorexia in this version
#######################################################################

f.sys.combine <- function(times,state,parms){
  # Variables
  R <- state["R"]
  L <- state["L"]
  M <- state["M"]
  I <- state["I"]
  C <- state["C"]

  # Parameters
  SR <- parms["SR"]
  r <- parms["r"]
  s <- parms["s"]
  v <- parms["v"]
  SL <- parms["SL"]
  g <- parms["g"]
  dL <- parms["dL"]
  dM <- parms["dM"]
  l <- parms["l"]
  p <- parms["p"]
  q <- parms["q"]
  a <- parms["a"]
  b <- parms["b"]
  w <- parms["w"]
  hL <- parms["hL"]
  hM <- parms["hM"]
  hI0 <- parms["hI0"]
  hI1 <- parms["hI1"]
  cc <- parms["c"]
  k <- parms["k0"]*cc/(1+parms["k1"]*cc)
  nuP <- parms["nuP"]
  nuL <- parms["nuL"]
  nuT <- 1-nuP-nuL

  # ODE
  FIP <- p+q*(L+s*M)*I
  DR <- SR - r*R - cc*FIP*R/(1+v*FIP*R)
  DL <- SL - (g+dL+nuL*k*I)*L
  DM <- g*L-(dM+nuP*k*I)*M
  DI <- FIP*R/(1+v*FIP*R) - l*I
  # Let condition be negative for the purposes of optimisation
  DC <- (a*r*R/(1+b*C)-w-(hL*L+hM*M)/(1+nuT*k*I)-(hI0*(nuL+nuP)*k+hI1*(nuL^2+nuP^2)*k^2)*I) #(nuL+nuP)*(hI0*k+hI1*k^2)


  return(list(c(DR,DL,DM,DI,DC)))
}

# Wrapper to calculate trajectory and return specified points
# cc, SR, SL function arguments
# Other parameters in environment
# hI is strategy-dependent
f.traj.combine <- function(parms.f,t.spec){
  # Named vector of parameters
  # Any parameters not specified in parms.f are given by parms.default
  parms <- c(parms.f,
             parms.default[!names(parms.default)%in%names(parms.f)]
  )                                     
  # Parasite-free initial conditions
  SR <- parms["SR"]
  r <- parms["r"]
  v <- parms["v"]
  p <- parms["p"]
  l <- parms["l"]
  a <- parms["a"]
  b <- parms["b"]
  w <- parms["w"]
  hI0 <- parms["hI0"]
  hI1 <- parms["hI1"]
  cc <- parms["c"]
  k <- parms["k0"]*cc/(1+parms["k1"]*cc)
  nuP <- parms["nuP"]
  nuL <- parms["nuL"]
  nuT <- 1-nuP-nuL
  
  # Initial conditions
  if(p>0){
    # Immune-mediated response, with constitutive represented by p>0
    R0 <- (-r-p*(cc-v*SR)+sqrt((r+p*(cc-v*SR))^2+4*r*v*p*SR))/(2*r*v*p)
    I0 <- p*R0/(1+v*p*R0)/l
  }else if (p==0){
    # Immune-mediated response, with constitutive represented by I0>0, p=0
    R0 <- SR/r
    I0 <- 0.01*r*R0/l
  }
  C0 <- (a*r*R0/(w+(hI0*(nuL+nuP)*k+hI1*(nuL^2+nuP^2)*k^2)*I0)-1)/b
  #C0 <-  pmin(pmax(0,C0),1)
  state.init <- c(R0,0,0,I0,C0)
  names(state.init) <- c("R","L","M","I","C")
  # Calculate trajectory
  traj <- data.frame(ode(y = state.init,
                         times = t,
                         func = f.sys.combine,
                         parms = parms,
                         method = "lsode"))

  # Return
  return(traj[traj$time %in% t.spec,])
}
#
ggplot(pivot_longer(f.traj.combine(c("SR"=3,"nuP"=0.1,"SL"=2),t.save),-time,values_to="Value",names_to="Variable"),
       aes(x=time,y=Value,color
           =Variable))+
  geom_line()

