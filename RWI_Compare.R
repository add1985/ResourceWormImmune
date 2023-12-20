
# Run optimisation code for single or combined strategies

library(tidyverse)
library(deSolve)
theme_set(theme_bw())

rm(list = ls()) 

#######################################################################
Optimise single strategies
#######################################################################

# Dimensional parameters
source("RWI_Parameters.R")

# Optimum investments
source("RWI_Optimise.R")
# write_csv(c.opt,"opt_harm0p8.csv")
# Check convergence; retry unconverged runs with different starting point
c.opt[which(c.opt$convergence!=0),]
c.opt.rerun <- bind_rows(
  apply(c.opt[which(c.opt$convergence!=0),],
        1,
        function(z) f.optim(strat=unlist(as.numeric(z["strat"])),
                            parms.f=setNames( # setNames allows conversion from dataframe into named vector
                              as.numeric(z[!names(z)%in%c("strat","c","C.Mean","R","M","L","I","C","convergence")]),
                              names(z)[!names(z)%in%c("strat","c","C.Mean","R","M","L","I","C","convergence")]),
                            c0=0.01) # 0.01, 10
        )
)
c.opt.rerun[which(c.opt.rerun$convergence!=0),]
c.opt[which(c.opt$convergence!=0),] <- c.opt.rerun
c.opt[which(c.opt$convergence!=0),]

# Save
#write.csv(c.opt,file="./c_opt_90_SL_hL_g0p5_Combine.csv",row.names=F)

# Or load previous run
#c.opt <- read.csv("opt_harm0p8.csv")
#c.opt$strat <- as.character(c.opt$strat)

# Add minimum survival point
f.MinSurvival <- function(strat,dat=c.opt){
  # Data for strategy
  f.dat <- dat[dat$strat==as.character(strat) & dat$Value>0,]
  # Min values of SR which host survives
  f.dat.MinSurvive <- data.frame()
  for (j.SL in unique(f.dat$SL)){
    for (j.hL in unique(f.dat$hL)){
      # Current value of SL, positive value
    f.dat.j <- f.dat[f.dat$SL==j.SL & f.dat$hL==j.hL & f.dat$Value>0,]
    f.dat.MinSurvive <- rbind(f.dat.MinSurvive,
                              f.dat.j[f.dat.j$Value==min(f.dat.j$Value),]
                              )
    }
    
  }
  # Corresponding points on x-axis
  f.dat.MinSurvive$Value <- 0
  f.dat.MinSurvive$c <- 0
  f.dat.MinSurvive$M <- 0
  return(f.dat.MinSurvive)
}
# Add to c.opt
c.opt <- rbind(bind_rows(lapply(1:4,f.MinSurvival)),c.opt)

# Remove massive investment values
#c.opt$c[c.opt$strat==1 & c.opt$c>10] <- NA

# Plot optimum investments for different environments and strategies
#c.opt <- read.csv("opt_harm1.csv") #
c.opt$strat <- as.character(c.opt$strat)
c.opt$hI0 <- as.character(c.opt$hI0)
c.opt[c.opt$C.Mean>=0 & c.opt$convergence==0 &c.opt$SL==0.5 
      & ((c.opt$hI0==0.25)|(c.opt$strat%in%c(1,4,0,1.5))) 
      & !c.opt$strat%in%c(1),] %>% #
  #pivot_longer(-c("strat","SR","SL","hL","hI0","hI1","convergence"),names_to="Variable",values_to="Value") %>%
  ggplot(aes(x=SR,y=c,color=strat))+#,color=strat))+,alpha=hI0
  geom_line(size=1)+
  #geom_line(aes(y=c),size=1,linetype="dashed")+
  facet_grid(g~hL,
               labeller=labeller(Variable=c("Value"="Mean host condition", 
                                         "M"="Final adult parasite load", 
                                         "c"="Optimal investment"),
                                 SL=c("0.5"="Low infection","2"="High infection"),
                                 hL=c("0.2"="Larvae less harmful","0.4"="Larvae equally harmful",
                                      "0.6"="Larvae more harmful"),
                                 g=c("0.5"="Rapid maturation","0.1"="Slower maturation")),
               scales="free_y")+
  labs(x="Resource availability",y="")+
  scale_x_continuous(limits=c(0,5), breaks=0:5,expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,2.5))+
  guides(fill="none")+ 
  theme(
     panel.grid = element_blank()
    )+
  scale_color_manual(values=c("#E74C3C","#F1948A","#F39C12", "#27AE60", "#2980B9","#17202A"), #
                      name="Strategy",
                      breaks=as.character(c(1,1.5,2:4,0)),
                      labels=c("Anorexia","Starvation" ,"Prevention", "Clearance", "Tolerance","None"))+
  ggtitle("g=0.5, v=2.5, q=0.5")

#######################################################################################



#######################################################################################
# Optimise combined strategies
#######################################################################################

source("RWI_Combine.R")
# write.csv(c.opt,file="./opt_harm0p8_Combine.csv",row.names=F)
# Check convergence; retry unconverged runs with different starting point
c.opt[which(c.opt$convergence!=0),]
c.opt.rerun <- bind_rows(
  apply(c.opt[which(c.opt$convergence!=0),],
        1,
        function(z) f.optim.combine(parms.f=setNames( # setNames allows conversion from dataframe into named vector
                              as.numeric(z[!names(z)%in%c("c","nuL","nuP","C.Mean","R","M","L","I","C","convergence")]),
                              names(z)[!names(z)%in%c("c","nuL","nuP","C.Mean","R","M","L","I","C","convergence")]),
                            c0=c(0.5,0.1,0.1)) # 0.01, 10
  )
)
c.opt.rerun[which(c.opt.rerun$convergence!=0),]
c.opt[which(c.opt$convergence!=0),] <- c.opt.rerun
c.opt[which(c.opt$convergence!=0),]


c.opt$hI0 <- as.character(c.opt$hI0)
c.opt$nuT <- 1-c.opt$nuL-c.opt$nuP
c.opt[c.opt$C.Mean>=0 & c.opt$convergence==0 &c.opt$SL==0.5,] %>% #
  #pivot_longer(-c("strat","SR","SL","hL","hI0","hI1","convergence"),names_to="Variable",values_to="Value") %>%
  ggplot(aes(x=SR,y=C.Mean,color=hI0))+#,color=strat))+,alpha=hI0
  geom_line(size=1)+
  #geom_line(aes(y=c),size=1,linetype="dashed")+
  facet_grid(g~hL,
             labeller=labeller(Variable=c("Value"="Mean host condition", 
                                          "M"="Final adult parasite load", 
                                          "c"="Optimal investment"),
                               SL=c("0.5"="Low infection","2"="High infection"),
                               hL=c("0.2"="Larvae less harmful","0.4"="Larvae equally harmful",
                                    "0.6"="Larvae more harmful"),
                               g=c("0.5"="Rapid maturation","0.1"="Slower maturation")),
             scales="free_y")+
  labs(x="Resource availability",y="")+
  scale_x_continuous(limits=c(0,5), breaks=0:5,expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,7))+
  guides(fill="none")+ 
  theme(
    panel.grid = element_blank()
  )+
  ggtitle("Combined")

c.opt[c.opt$C.Mean>=0 & c.opt$convergence==0 &c.opt$SL==0.5 ,] %>% #& c.opt$hI0==0.25
  pivot_longer(-names(c.opt)[!names(c.opt)%in%c("c","nuL","nuP","nuT")],names_to="Variable",values_to="Value") %>%
  #pivot_longer(-names(c.opt)[!names(c.opt)%in%c("nuL.c","nuP.c","nuT.c")],names_to="Variable",values_to="Value") %>%
  ggplot(aes(x=SR,y=Value,color=Variable,alpha=hI0))+#,color=strat))+
  geom_line(size=1)+
  #geom_line(aes(y=c),size=1,linetype="dashed")+
  facet_grid(g~hL,
             labeller=labeller(Variable=c("Value"="Mean host condition", 
                                          "M"="Final adult parasite load", 
                                          "c"="Optimal investment"),
                               SL=c("0.5"="Low infection","2"="High infection"),
                               hL=c("0.2"="Larvae less harmful","0.4"="Larvae equally harmful",
                                    "0.6"="Larvae more harmful"),
                               g=c("0.5"="Rapid maturation","0.1"="Slower maturation")),
             scales="free_y")+
  labs(x="Resource availability",y="")+
  scale_x_continuous(limits=c(0,5), breaks=0:5,expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.2))+
  guides(fill="none")+ 
  theme(
    panel.grid = element_blank()
  )+
  ggtitle("Combined")


# Manually change those that converged in the wrong place
c.opt.wrong <- c.opt[c.opt$hL==0.6 & c.opt$g==0.5 & c.opt$hI0==0.25 & c.opt$SL==0.5 & c.opt$C.Mean>0,][1,] #
#c.opt.wrong <- c.opt[c.opt$g==0.1  & c.opt$SL==0.5 & c.opt$C.Mean>0 & c.opt$SR<2,]
c.opt.rerun <- bind_rows(
  apply(c.opt.wrong,
        1,
        function(z) f.optim.combine(parms.f=setNames( # setNames allows conversion from dataframe into named vector
          as.numeric(z[!names(z)%in%c("c","nuL","nuP","C.Mean","R","M","L","I","C","convergence")]),
          names(z)[!names(z)%in%c("c","nuL","nuP","C.Mean","R","M","L","I","C","convergence")]),
          c0=c(0.15,0.09,0.2)) # 0.01, 10#
  )
)
c.opt.rerun$nuT <- 1-c.opt.rerun$nuL-c.opt.rerun$nuP
c.opt[449,]<-c.opt.rerun[,names(c.opt)]
