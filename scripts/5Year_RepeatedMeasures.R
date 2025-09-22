library(tidyverse)
library(lmerTest)
library(emmeans)
library(readxl)
library(sjPlot)
library(jtools)
library(ggeffects)
library(ggplot2)
library(multcomp)
library(MuMIn)
library(gtsummary)
library(viridis)
library(plyr)
library(pbkrtest)
library(glmmTMB)
library(DHARMa)
library(FD)
library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(ggordiplots)
library(data.table)

#read in cover for each year, add year column
cov.20<-read.csv(file.choose(), header=T, stringsAsFactors = F)
cov.20$Year<-2020

cov.21<-read.csv(file.choose(), header=T, stringsAsFactors = F)
cov.21$Year<-2021

cov.22<-read.csv(file.choose(), header=T, stringsAsFactors = F)
cov.22$Year<-2022

cov.23<-read.csv(file.choose(), header=T, stringsAsFactors = F)
cov.23$Year<-2023

cov.24<-read.csv(file.choose(), header=T, stringsAsFactors = F)
cov.24$Year<-2024


#combine for total cover across years and treatments
cov.CLIFF<-rbind.fill(cov.20, cov.21, cov.22, cov.23, cov.24)
cov.CLIFF$Year<-as.factor(cov.CLIFF$Year)
cov.CLIFF$plot<-as.factor(cov.CLIFF$plot)
cov.CLIFF$severity<-as.factor(cov.CLIFF$severity)
cov.CLIFF[is.na(cov.CLIFF)] <- 0 #Make NAs 0


#Trait analysis
# cov.traits<-read.csv(file.choose())#cover of trait species
# cov.stand<-wisconsin(cov.traits)#wisconsin standardized cover
traits<-read.csv(file.choose(), header=T, row.names=1)#traits
# traits<-log(traits[,1:3])#log transform

fd.out<-dbFD(traits, cov.traits)#run dbFD to get CWM
cwm.traits<-fd.out$CWM#isolate CWM, inspect
traits.dist<-vegdist(cwm.traits, method="euclidean")

cov.CLIFF<-cbind(cov.CLIFF, cwm.traits)#bind to mega dataframe


#This chunk of code controls how permutations are carried out and is essential
#for the repeated mesures analysis
CTRL.t <- how(within = Within(type = "free"), #restrict permutations for repeat measures
              plots = Plots(type = "none"),
              blocks=cov.CLIFF$plot,
              nperm = 999,
              observed = TRUE)

adonis.out<-adonis2(cwm.traits~plot+severity*Year, #treats time as split plot factor, plot as sample unit per Bakker 2024
                    data=cov.CLIFF, 
                    method="euclidean", 
                    permutations=CTRL.t,
                    by="margin")
adonis.out #interaction is significant, create new factor for pairwise adonis

cov.CLIFF$sev.year<-paste(cov.CLIFF$severity,cov.CLIFF$Year)

pairwise.out<-pairwise.adonis2(cwm.traits~sev.year, #treats time as split plot factor, plot as sample unit per Bakker 2024
                               data=cov.CLIFF, 
                               method="euclidean", 
                               by="margin")
pairwise.out #multiple significant pairwise comparisons with bonferonni correction

#Plot MDS of traits
#extract points and centroids for each severity/year combo
traits.mds<-metaMDS(traits.dist, weakties=F)
traits.scores = data.frame(MDS1 = traits.mds$points[,1], MDS2 = traits.mds$points[,2],severity=cov.CLIFF$severity, Year=cov.CLIFF$Year)

trait.centroids <- traits.scores %>%
  group_by(severity,Year) %>%
  dplyr::summarize(MDS1=mean(MDS1), MDS2=mean(MDS2))

#traits in ENVfit
trt.env<-envfit(traits.mds, cwm.traits, perm=999)
trt.env<-as.data.frame(trt.env$vectors$arrows*sqrt(trt.env$vectors$r))
trt.env$traits<-rownames(trt.env)


#Graph it!  Creates ordinations for each year with all severities and overlays
#envfit
traits.scores %>%
  ggplot(aes(x=MDS1, y=MDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.1) +
  geom_point(shape=1, size=2) +
  geom_point(data=trait.centroids, aes(MDS1, MDS2), size=1, shape=19) +
  facet_wrap(~Year,ncol=1)+
  scale_linetype_manual(values=c("solid", "dashed", "longdash"))+
  scale_color_manual(values=c("firebrick", "darkorange3", "darkolivegreen4" ))+
  scale_fill_manual(values=c("firebrick", "darkorange3", "darkolivegreen4" ))+
  geom_segment(data=trt.env,aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),colour="darkgrey", inherit.aes=FALSE) + 
  geom_text(data=trt.env,aes(x=NMDS1,y=NMDS2,label=traits),size=3, inherit.aes=FALSE)+
  theme_classic()


trt.env

#Graphs trat means by year and severity for each trait
####Graph trait means
cwms<-cov.CLIFF[, c(88,85:87)] 
cwms[,2:4]<-exp(cwms[,2:4])
cwms$sev.year<-as.factor(cwms$sev.year)
detach(package:plyr)

group.means<-cwms %>%
  group_by(sev.year) %>%
  summarize(mean.ldmc = mean(ldmc), SD.ldmc = sd(ldmc),
     mean.sla = mean(sla), SD.sla = sd(sla),
     mean.height=mean(height), SD.height=sd(height))
group.means<-separate_wider_delim(group.means, cols = sev.year, delim = " ", names = c("severity", "year"))

#########LDMC

ldmc.bar<-ggplot(group.means) +
  geom_bar( aes(x=severity, y=mean.ldmc),stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=severity, ymin=mean.ldmc-SD.ldmc, ymax=mean.ldmc+SD.ldmc), width=0.4, colour="black", alpha=0.9)+
  facet_wrap(~year, ncol=1)+
  labs(x = "Severity", y = bquote (LDMC~(gg^-1)))+
  theme_classic()
ldmc.bar

ldmc.aov<-aov(ldmc~severity*Year, data=cov.CLIFF)
emmeans(ldmc.aov, pairwise~severity|Year, adj=("mvt"))

####SLA


sla.bar<-ggplot(group.means) +
  geom_bar( aes(x=severity, y=mean.sla),stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=severity, ymin=mean.sla-SD.sla, ymax=mean.sla+SD.sla), width=0.4, colour="black", alpha=0.9)+
  facet_wrap(~year, ncol=1)+
  labs(x = "Severity", y = bquote (log(SLA)))+
  theme_classic()
sla.bar

sla.aov<-aov(sla~severity*Year, data=cov.CLIFF)
emmeans(sla.aov, pairwise~severity|Year, adj=("mvt"))

###########Height

height.bar<-ggplot(group.means) +
  geom_bar( aes(x=severity, y=mean.height),stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=severity, ymin=mean.height-SD.height, ymax=mean.height+SD.height), width=0.4, colour="black", alpha=0.9)+
  facet_wrap(~year, ncol=1)+
  labs(x = "Severity", y = bquote (log(Height)))+
  theme_classic()
height.bar

height.aov<-aov(height~severity*Year, data=cov.CLIFF)
emmeans(height.aov, pairwise~severity|Year, adj=("mvt"))
