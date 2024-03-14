######C.Garcia
#####10/08/2018
####Follow 'Simulations.R' script
###Plotting

#This script performs the plots of all the models calculated previously;
#8 modelled scenarios are plotted here: extinction (random & sensitivity, with and without compensation) & recovery (random and recovery ranking with and without compensation)

#***********************************Script start****************************************
rm(list=ls())
library(tidyverse)
library(cowplot)
library(mgcv)
library(hexbin)
library(RColorBrewer)



setwd('C:/Users/cg05/OneDrive - CEFAS/Science/Project Cefas - Other/MCZ/Extinction (Dave Clare 2017 - 20XX)/')
load('./output/RecoveryRANDOM_COMP.RData') # All 4 recovery models
load('./output/RecoveryRANDOM_NOCOMP.RData') # All 4 recovery models
load('./output/RecoveryRECOVERABILITY_COMP.RData') # All 4 recovery models
load('./output/RecoveryRECOVERABILITY_NOCOMP.RData') # All 4 recovery models

load('./output/ExtinctionRANDOM_COMP.RData') # All 4 extinction models
load('./output/ExtinctionRANDOM_NOCOMP.RData') # All 4 extinction models
load('./output/ExtinctionSensitivity_COMP.RData') # All 4 extinction models
load('./output/ExtinctionSensitivity_NOCOMP.RData') # All 4 extinction models

rm(fladen, output, output1, start)

## Colours
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(1000)

#***************************************************************************************
## Using geom_tile and density
#***************************************************************************************
RecRdmNC<-ggplot(RecoveryRANDOM_NOCOMP, aes(Nsp, log10(BPc))) +
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x, k=4), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Recovery Random & No Compensation')

## Recovery in recoverability order & No compensation
RecRecNC<-ggplot(RecoveryRECOVERABILITY_NOCOMP ,aes(x=Nsp,y=log10(BPc)))+
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x, k=4), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Recovery Recoverability & No Compensation')

## Random recovery & competitive compensation
RecRdmC<-ggplot(RecoveryRANDOM_COMP ,aes(x=Nsp,y=log10(BPc)))+
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x, k=4), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #geom_smooth(method=gam, alpha=0.5, se=TRUE, colour='red', linetype="twodash", size=1)+
  ggtitle('Recovery Random & Competitive Compensation')

## Recovery in 'recovery' order & competitive compensation
RecRecC<-ggplot(RecoveryRECOVERABILITY_COMP ,aes(x=Nsp,y=log10(BPc)))+
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x, k=4), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Recovery Recoverability & Competitive Compensation')


# Extinction plots
## Random extinction & No compensation
ExtRdmNC<-ggplot(ExtinctionRANDOM_NOCOMP ,aes(x=Nsp,y=log10(BPc)))+
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Extinction Random & No Compensation')


## Extinction by sensitivity & No compensation
ExtSenNC<-ggplot(ExtinctionSensitivity_NOCOMP ,aes(x=Nsp,y=log10(BPc)))+
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Extinction Sensitivity & No Compensation')


## Random extinction & competitive compensation
ExtRdmC<-ggplot(ExtinctionRANDOM_COMP ,aes(x=Nsp,y=log10(BPc)))+
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Extinction Random & Competitive Compensation')


## Extinction by sensitivity & competitive compensation
ExtSenC<-ggplot(ExtinctionSensitivity_COMP ,aes(x=Nsp,y=log10(BPc)))+
  geom_jitter(colour="grey",alpha=0.5)+
  stat_density2d(aes(fill=..level.., alpha=..level..),
                 size=3, bins=20, geom='polygon') +
  scale_fill_gradientn(colours=r)+
  #scale_y_continuous(limits = c(-3, 2.5))+
  stat_smooth(method = "gam", formula = y ~ s(x), colour="black", size = 1)+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Extinction Sensitivity & Competitive Compensation')


#Group Plots
plot_grid(RecRdmNC, ExtRdmNC, RecRecNC, ExtSenNC,  ncol=2)

#ggsave('./plot/No Compensation_Log10.png', 
#       width = 20, height = 20, units = "cm")


plot_grid(RecRdmC, ExtRdmC, RecRecC, ExtSenC,  ncol=2)


#ggsave('./plot/Compensation_Log10.png', 
#       width = 20, height = 20, units = "cm")
