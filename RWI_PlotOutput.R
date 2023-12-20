# Plot optimisation output

library(tidyverse)
library(deSolve)
theme_set(theme_bw())

rm(list = ls()) 

# Plotting functions
source("f_plot_opt.R")
source("f_plot_opt_nu.R")


##########################################################################

# Load optimisation results
c.opt <- read.csv("opt_harm0p8.csv")
c.opt.combine <- read.csv("opt_harm0p8_combine.csv")

# Remove anorexia if c=0
c.opt <- c.opt[!(c.opt$strat=="1" & c.opt$c<0.0001),]


# Data frame of nu values
c.opt.combine.nu <- c.opt.combine[,names(c.opt.combine)%in%c("C.Mean","convergence","g","SR","SL","hL","hI0","hI1","c","nuL","nuP","Value")]
c.opt.combine.nu$nuT <- 1-c.opt.combine.nu$nuL-c.opt.combine.nu$nuP
c.opt.combine.nu <- pivot_longer(c.opt.combine.nu,c("nuP","nuL","nuT"),names_to="Parameter",values_to="ParVal")
c.opt.combine.nu$ParValScaled <- c.opt.combine.nu$c*c.opt.combine.nu$ParVal

# Normal
# SL=2, hI=0.25
f.plot.opt(c.opt,C.Mean,2,c(0.25),"",c(1,5),c(0,5.1),"Optimised mean host condition",c.opt.combine)
f.plot.opt(c.opt,c,2,c(0.25),"",c(1,5),c(0,0.8),"Investment in response",c.opt.combine)
f.plot.opt(c.opt,M,2,c(0.25),"",c(1,5),c(0,70),"Final adult parasite burden",c.opt.combine)
f.plot.opt(c.opt,L,2,c(0.25),"",c(1,5),c(0,10.2),"Final larval parasite burden",c.opt.combine)
f.plot.opt(c.opt,I,2,c(0.25),"",c(1,5),c(0,22),"Final immune response level",c.opt.combine)
f.plot.opt(c.opt,R,2,c(0.25),"",c(1,5),c(0,10),"Final resource level",c.opt.combine)

# Plot combination parameters
f.plot.opt.nu(dat=c.opt.combine.nu,ParVal,2,0.25,title=NULL,XLim=c(1,5),YLim=c(0,1),YLab="Proportion of response allocated")
f.plot.opt.nu(dat=c.opt.combine.nu,ParValScaled,2,0.25,title=NULL,XLim=c(1,5),YLim=c(0,0.4),YLab="Proportion of response allocated")

# Normal; for talk (remove equal harm)
# SL=2, hI=0.25
f.plot.opt(c.opt[(c.opt$hL==0.6 & c.opt$g== 0.1) | (c.opt$hL==0.2 & c.opt$g== 0.5),],C.Mean,2,
           c(0.25),"",c(1,5),c(0,5.1),"Optimised mean host condition",
           c.opt.combine[(c.opt.combine$hL==0.6 & c.opt.combine$g== 0.1) | (c.opt.combine$hL==0.2 & c.opt.combine$g== 0.5),])
f.plot.opt.nu(dat=c.opt.combine.nu[c.opt.combine.nu$hL!=0.4,],ParVal,2,0.25,
              title=NULL,XLim=c(1,5),YLim=c(0,1),YLab="Proportion of response allocated")
f.plot.opt(c.opt[c.opt$hL==0.6 & c.opt$g==0.1,],c,2,c(0.25),"",c(1,5),c(0,0.8),
           "Investment in response",c.opt.combine[c.opt.combine$hL==0.6 & c.opt.combine$g==0.1,])
f.plot.opt(c.opt[c.opt$hL==0.6 & c.opt$g==0.1,],M,2,c(0.25),"",c(1,5),c(0,70),
           "Final adult parasite burden",c.opt.combine[c.opt.combine$hL==0.6 & c.opt.combine$g==0.1,])



# SL=0.5
f.plot.opt(c.opt,C.Mean,0.5,c(0.1,0.25,0.5),"",c(1,5),c(0,7),"Optimised mean host condition",c.opt.combine)
f.plot.opt(c.opt,c,0.5,c(0.1,0.25,0.5),"",c(1,5),c(0,0.8),"Investment in response",c.opt.combine)
f.plot.opt(c.opt,M,0.5,c(0.1,0.25,0.5),"",c(1,5),c(0,20),"Final adult parasite burden",c.opt.combine)

# One week
# Load optimisation results
c.opt <- read.csv("opt_harm0p8_t7.csv")

# Remove anorexia if c=0
c.opt <- c.opt[!(c.opt$strat=="1" & c.opt$c<0.0001),]

# Removed anorexia as it is really either none or starvation, but sometimes falsely converges at intermediate value
f.plot.opt(c.opt[c.opt$strat!=1,],C.Mean,2,c(0.25),"High infection pressure",c(1,5),c(0,9),"Optimised mean host condition")
f.plot.opt(c.opt[c.opt$strat!=1,],c,2,c(0.25),"High infection pressure",c(0,2),"Investment in response")

# One week, just two panels
# No strategy same as clearance in right-hand panel
#f.plot.opt(c.opt[c.opt$strat!=1 & c.opt$g==0.1 & c.opt$hL %in% c(0.2,0.6),],C.Mean,2,
#           c(0.25),"High infection pressure",c(1,5),c(0,9),"Optimised mean host condition")
f.plot.opt(c.opt[!(c.opt$strat==1 | (c.opt$strat==3 & c.opt$hL==0.6))  & c.opt$g==0.1 & c.opt$hL %in% c(0.2,0.6),],
           C.Mean,2,c(0.25),"",c(1,5),c(0,8),"Optimised mean host condition")
f.plot.opt(c.opt[!(c.opt$strat%in%c(0,1,1.5))  & c.opt$g==0.1 & c.opt$hL %in% c(0.2,0.6),],
           c,2,c(0.25),"",c(1,5),c(0,1.5),"Investment in response") # Haven't plotted starvation as it's infinite
f.plot.opt(c.opt[!(c.opt$strat==1 | (c.opt$strat==3 & c.opt$hL==0.6))  & c.opt$g==0.1 & c.opt$hL %in% c(0.2,0.6),],
           M,2,c(0.25),"",c(1,5),c(0,4),"Final adult parasite burden")

#################################################################################

# Parameter sensitivity

# Vary immune production
c.opt <- read.csv("opt_harm0p8.csv")
c.opt <- read.csv("opt_harm0p8_q0p5.csv")
c.opt <- read.csv("opt_harm0p8_q1.csv")
c.opt <- read.csv("opt_harm0p8_v2p5.csv")
c.opt <- read.csv("opt_harm0p8_v2p5_q0p5.csv")
c.opt <- read.csv("opt_harm0p8_v2p5_q1.csv")

f.plot.opt(c.opt,C.Mean,2,c(0.1,0.25,0.5),"",c(1,5),c(0,3.5),"Optimised mean host condition",NULL)

# Vary strength/investment relationship
c.opt <- rbind(
  data.frame(k0=1,k1=0.5,read.csv("opt_harm0p8.csv")),
  data.frame(k0=5,k1=0.5,read.csv("opt_harm0p8_k0_5.csv")),
  data.frame(k0=1,k1=2.5,read.csv("opt_harm0p8_k1_2p5.csv")),
  data.frame(k0=5,k1=2.5,read.csv("opt_harm0p8_k0_5_k1_2p5.csv"))
)

f.plot.opt.multi(c.opt,C.Mean,2,c(0.1,0.25,0.5),"",c(1,5),c(0,6),"Optimised mean host condition",NULL)

# Vary metabolism
c.opt <- read.csv("opt_harm0p8.csv")
c.opt <- read.csv("opt_harm0p8_a4w2.csv")
c.opt <- read.csv("opt_harm0p8_r5.csv")
c.opt <- read.csv("opt_harm0p8_a4w2_r5.csv")

f.plot.opt(c.opt[c.opt$strat!=1,],C.Mean,2,c(0.1,0.25,0.5),"",c(1,5),c(0,5),"Optimised mean host condition",NULL)



