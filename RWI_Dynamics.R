# Time series; heat maps and single solutions

setwd("C:/Users/andydean/Google Drive/SoaySheep/WithinHostModel")
# setwd("G:/My Drive/SoaySheep/WithinHostModel")

library(tidyverse)
library(deSolve)
theme_set(theme_bw())

rm(list = ls()) 

#######################################################################
#######################################################################

# Dimensional parameters
source("RWI_Parameters.R")

# Number of points
n.c <- length(cs)
SLs <- 2
n.SL <- length(SLs)
SRs <- c(3,5)
n.SR <- length(SRs)

# Initialise solution dataframes
sol.i <- data.frame()
sol.ii <- data.frame()
sol.iii <- data.frame()
sol.iv <- data.frame()


# Function to remove dead hosts
f.DeadHost <- function(traj){
  t.dead <- which(traj$C<=0)[1] # Point of death
  if(is.na(t.dead)){
    return(traj) # Do nothing if host doesn't die
  } else{
    traj.na <- traj
    traj.na[t.dead:nrow(traj),c("R","L","M","I","C" )] <- NA # Else remove data after death
    return(traj.na)
  }
}

# Loop through parameters
for (j1 in 1:n.c){
  for (j2 in 1:n.SR){
    for (j3 in 1:n.SL){
      sol.i <- rbind(sol.i,
                     f.DeadHost(data.frame("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3],
                                f.traj(1,parms.f=c("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3]),t.save))))
      sol.ii <- rbind(sol.ii,
                      f.DeadHost(data.frame("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3],
                                 f.traj(2,parms.f=c("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3]),t.save))))
      sol.iii <- rbind(sol.iii,
                       f.DeadHost(data.frame("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3],
                                  f.traj(3,parms.f=c("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3]),t.save))))
      sol.iv <- rbind(sol.iv,
                      f.DeadHost(data.frame("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3],
                                 f.traj(4,parms.f=c("c"=cs[j1],"SR"=SRs[j2],"SL"=SLs[j3]),t.save))))
    }
  }
}

#######################################################################

# Normalise 
vars.max <- c("R"=max(max(sol.i$R,sol.ii$R,sol.iii$R,sol.iv$R,na.rm=T)),
              "L"=max(max(sol.i$L,sol.ii$L,sol.iii$L,sol.iv$L,na.rm=T)),
              "M"=max(max(sol.i$M,sol.ii$M,sol.iii$M,sol.iv$M,na.rm=T)),
              "I"=max(max(sol.i$I,sol.ii$I,sol.iii$I,sol.iv$I,na.rm=T)),
              "C"=max(max(sol.i$C,sol.ii$C,sol.iii$C,sol.iv$C,na.rm=T))
)
f.normalise <- function(sol,maxes){
  sol$R <- sol$R/maxes["R"]
  sol$L <- sol$L/max(maxes["L"],maxes["M"])
  sol$M <- sol$M/max(maxes["L"],maxes["M"])
  sol$I <- sol$I/maxes["I"]
  sol$C <- sol$C/maxes["C"]
  return(sol)
}
sol.i.norm <- f.normalise(sol.i,vars.max)
sol.ii.norm <- f.normalise(sol.ii,vars.max)
sol.iii.norm <- f.normalise(sol.iii,vars.max)
sol.iv.norm <- f.normalise(sol.iv,vars.max)

# Long form
f.long <- function(sol){
  sol <- pivot_longer(sol,-c(time,c,SR,SL),names_to="Variable",values_to="Value")
}
sol.i <- f.long(sol.i)
sol.ii <- f.long(sol.ii)
sol.iii <- f.long(sol.iii)
sol.iv <- f.long(sol.iv)
sol.i.norm <- f.long(sol.i.norm)
sol.ii.norm <- f.long(sol.ii.norm)
sol.iii.norm <- f.long(sol.iii.norm)
sol.iv.norm <- f.long(sol.iv.norm)

#######################################################################

# Optimum strategies

# Load optimisation dataframe
c.opt <- read.csv("opt_harm0p8.csv")
# Round SR to avoid missing values
c.opt$SR <- round(c.opt$SR,4)
# Filter out irrelevant parameters
c.opt <- c.opt[c.opt$SR %in% SRs &
                 c.opt$SL %in%SLs &
                 c.opt$g==parms.default["g"] &
                 c.opt$hL==parms.default["hL"] &
                 c.opt$hM==parms.default["hM"] &
                 ( (c.opt$strat%in%c(2,3) & c.opt$hI0==parms.default["hI0"] & c.opt$hI1==parms.default["hI1"]) |
                     c.opt$strat%in%c(1,4) )
               ,]


#######################################################################

# Merge four strategies
sol.norm <- rbind(
   data.frame(strat="1",sol.i.norm),
   data.frame(strat="2",sol.ii.norm),
   data.frame(strat="3",sol.iii.norm),
   data.frame(strat="4",sol.iv.norm)
 )

# Plot one variable only
sol.plot.t.all.X <- function(X,SR,SL){
  solplot <- ggplot(sol.norm[sol.norm$Variable%in%c(X)
                          & sol.norm$SR %in% SR & sol.norm$SL %in% SL
                          ,],aes(x=time,y=c,fill=Value))+
    geom_raster()+
    # Optimal investment; only if the host survices
    #geom_hline(data=c.opt[c.opt$SR%in%SR & c.opt$SL%in%SL & c.opt$C.Mean>0,c("strat","c","SR")],aes(yintercept=c),size=1)+
    scale_fill_distiller(palette="Spectral")+#,limits=c(0,1)
    facet_grid(SR~strat, # 
               #scales="free",
               labeller=labeller(strat=c("1"="Anorexia","2"="Prevention","3"="Clearance","4"="Tolerance"),
                                 Variable=c(C="Condition",R="Resource",L="Larvae",M="Parasite",I="Immune"),
                                 SR=c("3"="Lower resource avilability","5"="Higher resource availability")))+
    labs(x="Time",y="Investment")+
    scale_x_continuous(breaks=c(0,30,60,90),expand=c(0,0))+
    scale_y_continuous(breaks=seq(0,2.5,0.5),expand=c(0,0))+
    guides(fill="none")+
    theme(
      panel.grid = element_blank(),
      strip.background=element_blank(),
      panel.spacing=unit(1,"lines"),
      text=element_text(size=20))
  #multiplot(p234,layout=matrix(c(1,2,2,2), nrow=1, byrow=TRUE))
  return(solplot)
}
sol.plot.t.all.X("C",5,2)
sol.plot.t.all.X("C",SRs,0.5)
sol.plot.t.all.X("R",SRs,2)
sol.plot.t.all.X("R",SRs,0.5)
sol.plot.t.all.X("M",5,2)
sol.plot.t.all.X("M",SRs,0.5)
sol.plot.t.all.X("L",SRs,2)
sol.plot.t.all.X("L",SRs,0.5)
sol.plot.t.all.X("I",SRs,2)
sol.plot.t.all.X("I",SRs,0.5)
