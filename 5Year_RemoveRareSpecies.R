library(tidyverse)
library(lmerTest)
library(emmeans)
library(readxl)
library(sjPlot)
library(jtools)
library(ggeffects)
library(ggplot2)
install.packages("multcomp")
library(multcomp)
library(MuMIn)
install.packages('gtsummary')
library(gtsummary)
library(viridis)
library(plyr)
library(pbkrtest)
install.packages("glmmTMB")
library(glmmTMB)
install.packages("DHARMa")
library(DHARMa)
library(FD)
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(ggordiplots)
library(data.table)

# Code that removes rare species from all years
# Binds new species and leaf traits to old trait table from 2024 data

cov.20<-read.csv("data/CoverSep2020.csv")
cov.20$Year<-2020

cov.21<-read.csv("data/CoverSep2021.csv")
cov.21$Year<-2021

cov.22<-read.csv("data/CoverSep2022.csv")
cov.22$Year<-2022

cov.23<-read.csv("data/CoverSep2023.csv")
cov.23$Year<-2023

cov.24<-read.csv("data/CoverSep2024.csv")
cov.24$Year<-2024

#Combine for total cover across years and treatments
library(plyr)
cov.CLIFF<-rbind.fill(cov.20, cov.21, cov.22, cov.23, cov.24)
detach("package:plyr", unload = TRUE)
cov.CLIFF$Year<-as.factor(cov.CLIFF$Year)
cov.CLIFF$plot<-as.factor(cov.CLIFF$plot)
cov.CLIFF$severity<-as.factor(cov.CLIFF$severity)
cov.CLIFF[is.na(cov.CLIFF)] <- 0 #Make NAs 0
cov.CLIFF <- cov.CLIFF[,colSums(cov.CLIFF != 0) > 0] #Remove species with zero cover
cov.CLIFF <- cov.CLIFF[,order(colnames(cov.CLIFF))] #Put columns in order to help reading species

#Removing rare species (species occurring on <5% of plots, rounded up)
yr <- c("2020", "2021", "2022", "2023", "2024") #Set vector of years
sp <- vector() #Create empty object
for(i in yr){

x <- cov.CLIFF %>% filter(Year==i) #Create trimmed master community matrix
  
plot_0.05 <- ceiling(nrow(x)*0.05) #Calculate 5% plot threshold

freq <- x %>%
  select(-plot, -severity, -Year) %>% #Remove plot, severity, and year columns for frequency calcuation
  mutate(across(everything(), ~replace(., .>0, 1))) %>% #Turn abundance into presence-absence
  colSums() #Add up presences to calculate frequency

sp[i] <- rbind(list(names(freq[freq > plot_0.05]))) #Filter species that meet frequency threshold
}

#Species list for each year
sp$'2020'
sp$'2021'
sp$'2022'
sp$'2023'
sp$'2024'

#Species list including all years
spp <- unique(unlist(sp))

#Find new species not already included from 2024
new.sp <- setdiff(spp,sp$'2024')

#Pull traits for new species from master trait matrix (separate from trait table below)
new.traits <- read_csv("data/TraitMatrix.csv") %>%
  filter(SPP %in% new.sp)
colnames(new.traits)[1] <- "spp" #Make sure first column is the same as the trait table

#Bind new species onto 2024 trait table
traits <- read_csv("data/TraitTable.csv")
library(plyr)
traits <- rbind.fill(traits,new.traits)
detach("package:plyr", unload = TRUE)


write.csv(traits,"data/5Year_TraitTable.csv")
