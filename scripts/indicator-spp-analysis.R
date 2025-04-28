# about -------------------------------------------------------------------
# author: maddie wallace
# date: 2-19-25
# run indicator species analysis across burn severity gradient

# fern edited and added a relative cover per severity calculation 4/7/2025

# packages ----------------------------------------------------------------
#install.packages("indicspecies")
#install.packages("tidyverse")
library(indicspecies)
library(tidyverse)
library(vegan)

# analysis ----------------------------------------------------------------
matrix <- read.csv("data/CommunityMatrix.csv")

# reorganize data for pakige
abund <- matrix[,3:ncol(matrix)]
severity <- matrix$Severity

# analyze
indicspp <- multipatt(abund, severity, 
                      func = "r.g", 
                      restcomb = c(1,2,3),
                      control = how(nperm=9999))
summary(indicspp, alpha = 0.1)

# get relative cover calculated

ind_spp <- c("VETH", "LOWR", "CEFE", "LIDA", "FEAR",
             "PSMA", "ELEL", "MUVI", "PIPR")
ind_sev <- c("H", "H", "H", "H", "H", "H", "H", "L", "U")
ind <- data.frame(ind_spp, ind_sev)
names(ind) <- c("spp", "severity")

ind_rel_cov <- matrix[,ind_spp] %>% 
  decostand(method = "total")
ind_rel_cov$severity <- matrix$Severity

ind_cov_means <- ind_rel_cov %>% 
  pivot_longer(cols = "VETH":"PIPR", names_to = "spp") %>% 
  group_by(severity, spp) %>% 
  summarise(cov = mean(value), cov_sd = sd(value)) %>% 
  mutate(cov_cv = cov_sd/cov)


ind_cov <- right_join(ind_cov_means, ind)

ggplot(ind_cov_means, aes(x = severity, y = cov*100))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymax = 100*(cov + cov_sd), ymin = 100*cov))+
  facet_wrap(~spp, scales = "free_x")



#visualize: proportion of plots in each seveerity with a presence of each indicator?
long_matrix <- matrix |> 
  pivot_longer(cols = 3:ncol(matrix), 
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
  group_by(species, Severity) %>% 
  filter(pres == "y") %>% 
  summarise(n_pres = n())




ind_spp <- c("VETH", "LOWR", "CEFE", "LIDA", "FEAR",
             "PSMA", "ELEL", "MUVI", "PIPR")



ggplot(ind_data, aes(x = Severity, y = presence, fill = Severity)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  facet_wrap(~ species) +
  labs(x = "severity", y = "mean proportion present")




