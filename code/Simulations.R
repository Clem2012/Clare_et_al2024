######C.Garcia
#####06/07/2018
####Follow 'MCZ_Recovery'
###Competitive compensation

#This script performs the extinction and recolonisation model for MCZ work;
#The extinction compensation model from Solan et al and Thomsen et al is used here and augmented with a recolonisation;
#Ranking sensitivity to extinction & Recolonisation was calculated from traits by D. Clare;
#Community is full from the start (or finish);
#6 models are defined here: random extinction, sensitivity and recovery each with and without compensation
#***********************************Script start****************************************
rm(list=ls())

#library
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

## New recoverability/sensitivity scores - Updated by D. Clare 01/05/2020
fladen<-read.csv('fladen2.csv')

head(fladen, 10)
start<-fladen

#*******************************************************************************************************************
#Recovery model 4 scenarios: 
# recovery ranking : Random; Compensation : No
# recovery ranking : Recoverability; Compensation: No
# recovery ranking : Random; Compensation : Yes
# recovery ranking : Recoverability; Compensation: Yes

#*******************************************************************************************************************
# For convenience store the number of species & total biomass as a seperate variable
nsp <- nrow(start)
BioToT<-sum(start$Bavg)

# Set the number of simulations
nsims <- 1000


## Choose whether biomass compensation is active or not (TRUE / FALSE)**********
COMP <- T
#*******************************************************************************

# Add extra columns to record the abundances, biomass and extinction/compensation probabilities
# during each simulation
start$AiSim <- NA
start$BiSim <- NA
start$CiSim <- NA
start$PiSim <- NA
start$RecSim<- NA

######
# Set the recovery probabilities
######
##Recovery, the highest values of RecProb is the most likely to recolonise
start$RecProb <- start$Recoverability / sum(start$Recoverability)
#start$RecProb <- 1/nsp
#*******************************************************************************

# Set up output dataframe to record the results. The columns are:
# Simulation    - the simulation number;
# Nsp           - the number of species remaining in the community;
# RecSp         - the species that recovers;
# AbnRec        - abundance by which the species has recovered;
# BioRec        - biomass by which the species has recovered;
# Measure       - the ecosystem function that is being calculated,
#                      e.g. OrgC or Chla or MD etc;
# Value         - the value of that ecosystem function;
# Note that the ExtSp & CompSp columns contain the species that will go extinct/compensate
# at the NEXT STEP, i.e. the first row records the function of the full
# community and details who is going to go extinct/compensate next.

output <- expand.grid(Simulation = 1:nsims, Nsp = 1:nsp,
                      RecSp=NA, AbnRec=NA, BioRec=NA, 
                      BPc=NA, aRPD=NA)
# Sort the output dataframe for clarity
output <- output[order(output$Simulation),]

# With everything now set up we can run the simulations

for (sim_count in 1:nsims){
  #browser()
  cat("sim_count: ", sim_count, "\n")
  
  # Reset abundances and probabilities for the current simulation
  start$AiSim <- 0
  start$BiSim <- 0
  start$RecSim <- start$RecProb
  start$CiSim <- 0
  start$PiSim<- 0
  # Count up from 0 to full community,
  
  for (sp_count in 1:nsp)
  {
    # Define randomly which species recover first based on probability
    Recovery <- which(cumsum(start$RecSim)>=runif(1))[1]
    
    # Call compensation function before we gain the recovered species
    if(COMP==TRUE)
    {
      # Define by how much the taxa is recovering by with compensation
      start$CiSim[Recovery] <- start[Recovery, 'Bavg']
      start$PiSim<-start$CiSim/sum(start$CiSim)
      start$BiSim<-BioToT * start$PiSim
      start$AiSim<-start$BiSim/start$Bind
    } else
    {
      start[Recovery, 'BiSim']<-start[Recovery, 'Bavg']
      start[Recovery, 'AiSim']<-start[Recovery, 'Aavg']
    }
    
    # Calculate total BPc of community
    BPc <- sum(start$AiSim * (sqrt(start$Bind) * start$Mi * start$Ri), na.rm=T)
    
    # Calculate functions for community
    aRPD <- exp((-0.4666*log(BPc))+4.2151)
    
    
    # Store these results in the output dataframe
    output[output$Simulation == sim_count & output$Nsp==sp_count, "BPc"] <- BPc
    output[output$Simulation == sim_count & output$Nsp==sp_count, "RecSp"] <- as.character(start$Taxa[Recovery])
    output[output$Simulation == sim_count & output$Nsp==sp_count, "AbnRec"] <- start$AiSim[Recovery]
    output[output$Simulation == sim_count & output$Nsp==sp_count, "BioRec"] <- start$BiSim[Recovery]
    output[output$Simulation == sim_count & output$Nsp==sp_count, "aRPD"] <- aRPD
    
    # Recovery happened! Probability of recovery species to 0
    start[Recovery,c('RecSim')] <- 0
    
    # Recovery probability reset for the rest of the species
    start$RecSim<- start$RecSim / sum(start$RecSim)
  }
}

output1<-output

# No Compensation
#RecoveryRANDOM_NOCOMP<-rbind(output1, output)
RecoveryRECOVERABILITY_NOCOMP<-rbind(output1, output)
#save.image('./output/RecoveryRANDOM_NOCOMP.RData')
save.image('./output/RecoveryRECOVERABILITY_NOCOMP.RData')

# Compensation
#RecoveryRANDOM_COMP<-rbind(output1, output)
RecoveryRECOVERABILITY_COMP<-rbind(output1, output)
#save.image('./output/RecoveryRANDOM_COMP.RData')
save.image('./output/RecoveryRECOVERABILITY_COMP.RData')

#****************************************************
rm(fladen, output, output1, start)
save.image('./output/RecoveryALL.RData')


#*******************************************************************************************************************
#Extinction model 1 Sensitivity ranking (order of extinction) and biomass proportion (competition) 
#*******************************************************************************************************************
start<-fladen
# For convenience store the number of species as a seperate variable
nsp <- nrow(start)
BioToT<-sum(start$Bavg)

COMP<-T
# Set the number of simulations
nsims <- 200

# Add extra columns to record the abundances, biomass and extinction/compensation probabilities
# during each simulation
start$AiSim <- NA
start$BiSim <- NA
start$CiSim <- NA
start$EPSim <- NA

######
# Set the extinction probabilities to :
# Exinction Model - Sensitivity
######
##The highest values of ExtProb is the most likely to go extinct
#start$ExtProb <- start$Sensitivity / sum(start$Sensitivity) # Extinction model - Sensitivity
start$ExtProb <- 1/nsp


# Set up output dataframe to record the results. The columns are:
# Simulation    - the simulation number;
# Nsp           - the number of species remaining in the community;
# ExtSp         - the species that will go extinct at the next step;
# CompSp        - the species that will compensate at the next step;           
# Measure       - the ecosystem function that is being calculated,
#                      e.g. OrgC or Chla or MD etc;
# Value         - the value of that ecosystem function;
# Note that the ExtSp & CompSp columns contain the species that will go extinct/compensate
# at the NEXT STEP, i.e. the first row records the function of the full
# community and details who is going to go extinct/compensate next.

output <- expand.grid(Simulation = 1:nsims, Nsp = nsp:1,
                      ExtSp=NA, BioExt=NA, 
                      BPc=NA, aRPD=NA)
# Sort the output dataframe for clarity
output <- output[order(output$Simulation),]

# With everything now set up we can run the simulations
for (sim_count in 1:nsims){
  #browser()
  cat("sim_count: ", sim_count, "\n")
  # Reset abundances and probabilities for the current simulation
  start$AiSim <- start$Aavg
  start$BiSim <- start$Bavg
  start$EPSim <- start$ExtProb
  start$PiSim <- start$Bavg/sum(start$Bavg)
  start$CiSim <- 0
  
  
  # Count down from full number of species to 1,
  # knocking out species according to extinction probability 
  
  for (sp_count in nsp:1)
  {
    # Calculate total BPc & BQPc potential of community
    BPc <- sum(start$AiSim * (sqrt(start$Bind) * start$Mi * start$Ri), na.rm=T)
    
    # Calculate functions for community
    aRPD <- exp((-0.4666*log(BPc))+4.2151)
    
    
    # Store these results in the output dataframe
    output[output$Simulation == sim_count & output$Nsp==sp_count, "BPc"] <- BPc
    output[output$Simulation == sim_count & output$Nsp==sp_count, "aRPD"] <- aRPD
    
    # Randomly pick a species to go extinct based on probability
    Extinct <- which(cumsum(start$EPSim)>=runif(1))[1]
    
    # How much BIOMASS will be lost with the doomed species
    BiomassLost <- start[Extinct,"BiSim"]
    # Species going extinct cannot compensate
    start[Extinct,c("AiSim", "BiSim", 'EPSim',"PiSim")] <- 0
    
    # Call compensation function before we lose the extinct species
    if(COMP==TRUE)
    {
      start$PiSim<- start$PiSim / sum(start$PiSim) # Normalise
      
      #Compensate
      start$CiSim<- BiomassLost * start$PiSim 
      #start$BiSim<- start$BiSim + start$CiSim
      start$BiSim<- start$BiSim + (BiomassLost * start$PiSim)
      start$AiSim<- start$BiSim/start$Bind
    } 
    
    # Record ID of who has gone extinct and compensated
    output[output$Simulation == sim_count & output$Nsp==sp_count,"ExtSp"] <- as.character(start$Taxa[Extinct])
    
    # Record amount of compensation biomass
    output[output$Simulation == sim_count & output$Nsp==sp_count,"BioExt"] <-BiomassLost
    
    
    # Extinction probability reset for newly compensating species
    start$EPSim<- start$EPSim / sum(start$EPSim)
  }
}

output1<-output

# No Compensation
#ExtinctionRANDOM_NOCOMP<-rbind(output1, output)
ExtinctionSensitivity_NOCOMP<-rbind(output1, output)
rm(list=setdiff(ls(), "ExtinctionRANDOM_NOCOMP"))

#save.image('./output/ExtinctionRANDOM_NOCOMP.RData')
save.image('./output/ExtinctionSensitivity_NOCOMP.RData')

# Compensation
ExtinctionRANDOM_COMP<-rbind(output1, output)
#ExtinctionSensitivity_COMP<-rbind(output1, output)
save.image('./output/ExtinctionRANDOM_COMP.RData')
#save.image('./output/ExtinctionSensitivity_COMP.RData')

#****************************************************

save.image('./output/ExtinctionALL.RData')


#******************************************

#******************************************
ggplot(ExtinctionRANDOM_NOCOMP ,aes(x=Nsp,y=log(BPc)))+
  geom_point(colour="grey",alpha=0.1)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  theme_bw()+
  theme(legend.position="none")


ggsave('./plot/sensitivity_competitiveness.png', 
       width = 20, height = 20, units = "cm")
#*******************************************************************************************************************
#Extinction model 2 random ranking (order of extinction) and biomass proportion (competitivity) 
#*******************************************************************************************************************
start<-fladen
# For convenience store the number of species as a seperate variable
nsp <- nrow(start)
BioToT<-sum(start$Bavg)
COMP<-TRUE
# Set the number of simulations
nsims <- 500

# Add extra columns to record the abundances, biomass and extinction/compensation probabilities
# during each simulation
start$AiSim <- NA
start$BiSim <- NA
start$CiSim <- NA
start$EPSim <- NA

######
# Set the extinction probabilities to :
# Exinction Model - Sensitivity
######
##The highest values of ExtProb is the most likely to go extinct
start$ExtProb <- 1/nsp # Extinction model - Random


# Set up output dataframe to record the results. The columns are:
# Simulation    - the simulation number;
# Nsp           - the number of species remaining in the community;
# ExtSp         - the species that will go extinct at the next step;
# CompSp        - the species that will compensate at the next step;           
# Measure       - the ecosystem function that is being calculated,
#                      e.g. OrgC or Chla or MD etc;
# Value         - the value of that ecosystem function;
# Note that the ExtSp & CompSp columns contain the species that will go extinct/compensate
# at the NEXT STEP, i.e. the first row records the function of the full
# community and details who is going to go extinct/compensate next.

output <- expand.grid(Simulation = 1:nsims, Nsp = nsp:1,
                      ExtSp=NA, BioExt=NA, 
                      BPc=NA, aRPD=NA)
# Sort the output dataframe for clarity
output <- output[order(output$Simulation),]

# With everything now set up we can run the simulations
for (sim_count in 1:nsims){
  #browser()
  cat("sim_count: ", sim_count, "\n")
  # Reset abundances and probabilities for the current simulation
  start$AiSim <- start$Aavg
  start$BiSim <- start$Bavg
  start$EPSim <- start$ExtProb
  start$PiSim <- start$Bavg/sum(start$Bavg)
  start$CiSim <- 0
  
  
  # Count down from full number of species to 1,
  # knocking out species according to extinction probability 
  
  for (sp_count in nsp:2)
  {
    # Calculate total BPc & BQPc potential of community
    BPc <- sum(start$AiSim * (sqrt(start$Bind) * start$Mi * start$Ri), na.rm=T)
    
    # Calculate functions for community
    aRPD <- exp((-0.4666*log(BPc))+4.2151)
    
    
    # Store these results in the output dataframe
    output[output$Simulation == sim_count & output$Nsp==sp_count, "BPc"] <- BPc
    output[output$Simulation == sim_count & output$Nsp==sp_count, "aRPD"] <- aRPD
    
    # Randomly pick a species to go extinct based on probability
    Extinct <- which(cumsum(start$EPSim)>=runif(1))[1]
    
    # Call compensation function before we lose the extinct species
    if(COMP==TRUE)
    {
      # How much BIOMASS will be lost with the doomed species
      BiomassLost <- start[Extinct,"BiSim"]
      
      # Picking the compensatory species
      # Species going extinct cannot compensate
      start[Extinct,"PiSim"] <- 0
      start$PiSim<- start$PiSim / sum(start$PiSim) # Normalise
      
      #Compensate
      start$CiSim<- start$CiSim + (BiomassLost * start$PiSim) 
      start$BiSim<- start$BiSim + (BiomassLost * start$PiSim)
      start$AiSim<- start$AiSim + start$BiSim/start$Bind
    } 
    
    # Record ID of who has gone extinct and compensated
    output[output$Simulation == sim_count & output$Nsp==sp_count,"ExtSp"] <- as.character(start$Taxa[Extinct])
    
    # Record amount of compensation biomass
    output[output$Simulation == sim_count & output$Nsp==sp_count,"BioExt"] <-BiomassLost
    
    # Extinction happens! Set abundance and probability of doomed species to 0
    start[Extinct,c("AiSim", "BiSim", 'EPSim')] <- 0
    
    # Extinction probability reset for newly compensating species
    start$EPSim<- start$EPSim / sum(start$EPSim)
  }
}
random_competitiveness<-output

save.image('./output/outputRevrecovery.RData')
#******************************************
ggplot(random_competitiveness ,aes(x=Nsp,y=log(BPc)))+
  geom_point(colour="grey",alpha=0.1)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  theme_bw()+
  theme(legend.position="none")

ggsave('./plot/random_competitiveness.png', 
       width = 20, height = 20, units = "cm")

#*******************************************************************************************************************
#Recovery model: random ranking (order of appearance) and competitivity (biomass proportion of the full community) 
#*******************************************************************************************************************
start<-fladen
# For convenience store the number of species & total biomass as a seperate variable
nsp <- nrow(start)
BioToT<-sum(start$Bavg)

# Set the number of simulations
nsims <- 500
# Choose whether biomass compensation is active or not (TRUE / FALSE)
COMP <- FALSE

# Add extra columns to record the abundances, biomass and extinction/compensation probabilities
# during each simulation
start$AiSim <- NA
start$BiSim <- NA
start$CiSim <- NA
start$PiSim <- NA
start$RecSim<- NA

######
# Set the recovery probabilities
######
##Recovery, the highest values of RecProb is the most likely to recolonise
#start$RecProb <- 1/nsp # Recovery model - Random
#start$RecProb <- start$Recovery / sum(start$Recovery) # Recovery model - Recovery ranking
start$RecProb <- (1/start$Sensitivity) / sum(1/start$Sensitivity) # Recovery model  - Reverse Sensitivity



# Set up output dataframe to record the results. The columns are:
# Simulation    - the simulation number;
# Nsp           - the number of species remaining in the community;
# RecSp         - the species that recovers;
# AbnRec        - abundance by which the species has recovered;
# BioRec        - biomass by which the species has recovered;
# Measure       - the ecosystem function that is being calculated,
#                      e.g. OrgC or Chla or MD etc;
# Value         - the value of that ecosystem function;
# Note that the ExtSp & CompSp columns contain the species that will go extinct/compensate
# at the NEXT STEP, i.e. the first row records the function of the full
# community and details who is going to go extinct/compensate next.

output <- expand.grid(Simulation = 1:nsims, Nsp = 1:nsp,
                      RecSp=NA, AbnRec=NA, BioRec=NA, 
                      BPc=NA, aRPD=NA)
# Sort the output dataframe for clarity
output <- output[order(output$Simulation),]

# With everything now set up we can run the simulations

for (sim_count in 1:nsims){
  #browser()
  cat("sim_count: ", sim_count, "\n")
  
  # Reset abundances and probabilities for the current simulation
  start$AiSim <- 0
  start$BiSim <- 0
  start$RecSim <- start$RecProb
  start$CiSim <- 0
  start$PiSim<- 0
  # Count up from 0 to full community,
  
  for (sp_count in 1:nsp)
  {
    # Define randomly which species recover first based on probability
    Recovery <- which(cumsum(start$RecSim)>=runif(1))[1]
    
    # Call compensation function before we lose the extinct species
    if(COMP==TRUE)
    {
      # Define by how much the taxa is recovering by with compensation
      start$CiSim[Recovery] <- start[Recovery, 'Bavg']
      start$PiSim<-start$CiSim/sum(start$CiSim)
      start$BiSim<-BioToT * start$PiSim
      start$AiSim<-start$BiSim/start$Bind
    } else
    {
      start[Recovery, 'BiSim']<-start[Recovery, 'Bavg']
      start[Recovery, 'AiSim']<-start[Recovery, 'Aavg']
    }
    
    # Calculate total BPc of community
    BPc <- sum(start$AiSim * (sqrt(start$Bind) * start$Mi * start$Ri), na.rm=T)
    
    # Calculate functions for community
    aRPD <- exp((-0.4666*log(BPc))+4.2151)
    
    
    # Store these results in the output dataframe
    output[output$Simulation == sim_count & output$Nsp==sp_count, "BPc"] <- BPc
    output[output$Simulation == sim_count & output$Nsp==sp_count, "RecSp"] <- as.character(start$Taxa[Recovery])
    output[output$Simulation == sim_count & output$Nsp==sp_count, "AbnRec"] <- start$AiSim[Recovery]
    output[output$Simulation == sim_count & output$Nsp==sp_count, "BioRec"] <- start$BiSim[Recovery]
    output[output$Simulation == sim_count & output$Nsp==sp_count, "aRPD"] <- aRPD
    
    # Recovery happened! Probability of recovery species to 0
    start[Recovery,c('RecSim')] <- 0
    
    # Recovery probability reset for the rest of the species
    start$RecSim<- start$RecSim / sum(start$RecSim)
  }
}
Recovery_RevSensitivity<-output
save.image('./output/outputRecovery.RData')

#******************************************
ggplot(Recovery_RevSensitivity ,aes(x=Nsp,y=log(BPc)))+
  geom_point(colour="grey",alpha=0.1)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  theme_bw()+
  theme(legend.position="none")

ggsave('./plot/Recovery_RevSensitivityNOCOMP.png', 
       width = 20, height = 20, units = "cm")

