# about -------------------------------------------------------------------
# author: maddie wallace
# date: 2-19-25

# run indicator species analysis across burn severity gradient
# fern edited and added a relative cover per severity calculation 4/7/2025
# maddie added year 1 - 5 analysis on 10/2/2025

# packages ----------------------------------------------------------------
#install.packages("indicspecies")
#install.packages("tidyverse")
library(indicspecies)
library(tidyverse)
library(vegan)

# analysis ----------------------------------------------------------------
# get data for y1-5
cov.24 <- read.csv("data/CoverSep2024.csv")
cov.23 <- read.csv("data/CoverSep2023.csv")
cov.22 <- read.csv("data/CoverSep2022.csv")
cov.21 <- read.csv("data/CoverSep2021.csv")
cov.20 <- read.csv("data/CoverSep2020.csv")


# reorganize data for pakige
abund.24 <- cov.24[,3:ncol(cov.24)]
abund.24[is.na(abund.24)] <- 0
severity.24 <- cov.24$severity

abund.23 <- cov.23[,3:ncol(cov.23)]
abund.23[is.na(abund.23)] <- 0
severity.23 <- cov.23$severity

abund.22 <- cov.22[,3:ncol(cov.22)]
abund.22[is.na(abund.22)] <- 0
severity.22 <- cov.22$severity

abund.21 <- cov.21[,3:ncol(cov.21)]
abund.21[is.na(abund.21)] <- 0
severity.21 <- cov.21$severity

abund.20 <- cov.20[,3:ncol(cov.20)]
abund.20[is.na(abund.20)] <- 0
severity.20 <- cov.20$severity

# analyze 
indicspp.24 <- multipatt(abund.24, severity.24, 
                      func = "r.g", 
                      restcomb = c(1,2,3),
                      control = how(nperm=9999))
summary(indicspp.24, alpha = 0.1)


indicspp.23 <- multipatt(abund.23, severity.23, 
                         func = "r.g", 
                         restcomb = c(1,2,3),
                         control = how(nperm=9999))
summary(indicspp.23, alpha = 0.1)


indicspp.22 <- multipatt(abund.22, severity.22, 
                         func = "r.g", 
                         restcomb = c(1,2,3),
                         control = how(nperm=9999))
summary(indicspp.22, alpha = 0.1)


indicspp.21 <- multipatt(abund.21, severity.21, 
                         func = "r.g", 
                         restcomb = c(1,2,3),
                         control = how(nperm=9999))
summary(indicspp.21, alpha = 0.1)


indicspp.20 <- multipatt(abund.20, severity.20, 
                         func = "r.g", 
                         restcomb = c(1,2,3),
                         control = how(nperm=9999))
summary(indicspp.20, alpha = 0.1)

# get relative cover calculated for each year!
#2024
ind_spp.24 <- c("VETH", "LOWR", "CEFE", "LIDA", "FEAR",
             "PSMA", "ELEL", "MUVI", "PIPR")
ind_sev.24 <- c("H", "H", "H", "H", "H", "H", "H", "L", "U")
ind.24 <- data.frame(ind_spp.24, ind_sev.24)
names(ind.24) <- c("spp", "severity")

ind_rel_cov.24 <- cov.24[,ind_spp.24] %>% 
  decostand(method = "total")
ind_rel_cov.24$severity <- cov.24$severity

ind_cov_means.24 <- ind_rel_cov.24 %>% 
  pivot_longer(cols = "VETH":"PIPR", names_to = "spp") %>% 
  group_by(severity, spp) %>% 
  summarise(cov = mean(value), cov_sd = sd(value)) %>% 
  mutate(cov_cv = cov_sd/cov)

ind_cov.24 <- right_join(ind_cov_means.24, ind.24)

#2023
ind_spp.23 <- c("SCSC", "MUST", "LIMU", "BADI", "COCA",
                "VETH", "LASE", "SATR", "LIDA", "CHAL", "LOWR")
ind_sev.23 <- c("U", "L", "L", "L", "H", "H", "H", "H", "H", "H", "H")
ind.23 <- data.frame(ind_spp.23, ind_sev.23)
names(ind.23) <- c("spp", "severity")


ind_rel_cov.23 <- cov.23[, ind_spp.23] %>%
  replace(is.na(.), 0) %>%
  decostand(method = "total") %>%
  as.data.frame()

ind_rel_cov.23$severity <- cov.23$severity

ind_cov_means.23 <- ind_rel_cov.23 %>% 
  pivot_longer(cols = all_of(ind_spp.23), names_to = "spp") %>% 
  group_by(severity, spp) %>% 
  summarise(cov = mean(value), cov_sd = sd(value)) %>% 
  mutate(cov_cv = cov_sd/cov)

ind_cov.23 <- right_join(ind_cov_means.23, ind.23)


#2022
ind_spp.22 <- c("PIPR", "HOWR", "MUST", "BADI", "ERDI",
                "MUMO", "SATR", "VETH", "COCA", "CRGR", "CEFE")
ind_sev.22 <- c("U", "U", "L", "L", "L", "L", "H", "H", "H", "H", "H")
ind.22 <- data.frame(ind_spp.22, ind_sev.22)
names(ind.22) <- c("spp", "severity")


ind_rel_cov.22 <- cov.22[, ind_spp.22] %>%
  replace(is.na(.), 0) %>%
  decostand(method = "total") %>%
  as.data.frame()

ind_rel_cov.22$severity <- cov.22$severity

ind_cov_means.22 <- ind_rel_cov.22 %>% 
  pivot_longer(cols = all_of(ind_spp.22), names_to = "spp") %>% 
  group_by(severity, spp) %>% 
  summarise(cov = mean(value), cov_sd = sd(value)) %>% 
  mutate(cov_cv = cov_sd/cov)

ind_cov.22 <- right_join(ind_cov_means.22, ind.22)


#2021
ind_spp.21 <- c("PIPR", "MUST", "MUMO", "HEMU", "SATR",
                "VETH", "LASE", "CEFE", "LOWR", "CHAL")
ind_sev.21 <- c("U", "L", "L", "L", "H", "H", "H", "H", "H", "H")
ind.21 <- data.frame(ind_spp.21, ind_sev.21)
names(ind.21) <- c("spp", "severity")


ind_rel_cov.21 <- cov.21[, ind_spp.21] %>%
  replace(is.na(.), 0) %>%
  decostand(method = "total") %>%
  as.data.frame()

ind_rel_cov.21$severity <- cov.21$severity

ind_cov_means.21 <- ind_rel_cov.21 %>% 
  pivot_longer(cols = all_of(ind_spp.21), names_to = "spp") %>% 
  group_by(severity, spp) %>% 
  summarise(cov = mean(value), cov_sd = sd(value)) %>% 
  mutate(cov_cv = cov_sd/cov)

ind_cov.21 <- right_join(ind_cov_means.21, ind.21)

#2020
ind_spp.20 <- c("PIPR", "MUMO", "CEFE", "ELEL", "LOWR")
ind_sev.20 <- c("U", "L", "H", "H", "H")
ind.20 <- data.frame(ind_spp.20, ind_sev.20)
names(ind.20) <- c("spp", "severity")


ind_rel_cov.20 <- cov.20[, ind_spp.20] %>%
  replace(is.na(.), 0) %>%
  decostand(method = "total") %>%
  as.data.frame()

ind_rel_cov.20$severity <- cov.20$severity

ind_cov_means.20 <- ind_rel_cov.20 %>% 
  pivot_longer(cols = all_of(ind_spp.20), names_to = "spp") %>% 
  group_by(severity, spp) %>% 
  summarise(cov = mean(value), cov_sd = sd(value)) %>% 
  mutate(cov_cv = cov_sd/cov)

ind_cov.20 <- right_join(ind_cov_means.20, ind.20)


--------------------------------------------------------------------------------
# visualize
ggplot(ind_cov_means, aes(x = severity, y = cov*100))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymax = 100*(cov + cov_sd), ymin = 100*cov))+
  facet_wrap(~spp, scales = "free_x")


#visualize: proportion of plots in each seveerity with a presence of each indicator?
long_matrix <- cov.24 |> 
  pivot_longer(cols = 3:ncol(cov.24), 
               names_to = "species", 
               values_to = "presence") %>% 
  mutate(pres = case_when(presence == 0 ~ "n",
                          presence > 0 ~ "y"))

ind_data <- long_matrix |> 
  filter(species %in% ind_spp)

# number of plots:
# U 20
# L 18
# H 19

ind_sev_n <- ind_data %>% 
  group_by(species, severity) %>% 
  filter(pres == "y") %>% 
  summarise(n_pres = n())

ind_spp <- c("VETH", "LOWR", "CEFE", "LIDA", "FEAR",
             "PSMA", "ELEL", "MUVI", "PIPR")

ggplot(ind_data, aes(x = severity, y = presence, fill = severity)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  facet_wrap(~ species) +
  labs(x = "severity", y = "mean proportion present")




