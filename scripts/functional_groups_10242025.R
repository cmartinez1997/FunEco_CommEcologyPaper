library(tidyverse)
library(vegan)
library(FD)
library(wesanderson)
library(cowplot)
library(pairwiseAdonis)


#read in cover for each year, add year column
cov.20<-read_csv("data/CoverSep2020.csv")
cov.20$Year<-2020

cov.21<-read_csv("data/CoverSep2021.csv")
cov.21$Year<-2021

cov.22<-read_csv("data/CoverSep2022.csv")
cov.22$Year<-2022

cov.23<-read_csv("data/CoverSep2023.csv")
cov.23$Year<-2023

cov.24<-read_csv("data/CoverSep2024.csv")
cov.24$Year<-2024


#combine for total cover across years and treatments
library(plyr) 
cov.CLIFF<-rbind.fill(cov.20, cov.21, cov.22, cov.23, cov.24)
detach("package:plyr", unload = TRUE) ##### included this because i (ian) have had issues with plyr functions conflicting with functions in other packages 
cov.CLIFF$Year<-as.factor(cov.CLIFF$Year)
cov.CLIFF$plot<-as.factor(cov.CLIFF$plot)
cov.CLIFF$severity<-as.factor(cov.CLIFF$severity)
cov.CLIFF[is.na(cov.CLIFF)] <- 0 #Make NAs 0

year <- cov.CLIFF$Year
plot <- cov.CLIFF$plot
severity <- cov.CLIFF$severity

#read in spp that are not rare (that we have traits for)

traits <- read_csv("data/5Year_TraitTable.csv") %>% 
  mutate(spp = case_when(spp == "MUVI" ~ "MUST", .default = spp)) %>% 
  arrange(spp) %>%
  column_to_rownames(var="spp") %>% ##### ian was here #####
  select(sla, height, seedmass, resprouting) %>% ##### ian was here #####
  mutate(resprouting=resprouting+1) ##### ian was here, again. making resprouting 1-2 #####

spp_pool <- rownames(traits) ##### ian was here #####
cov.CLIFF <- cov.CLIFF[,(which(names(cov.CLIFF) %in% spp_pool))]
cov.CLIFF <- cov.CLIFF[,order(names(cov.CLIFF))]
cov.CLIFF$Year <- year ; cov.CLIFF$plot <- plot 
cov.CLIFF$severity <- severity

#steps:
# for each year, i want to aggregate values by functional type (either groups or nativity)
# then for each year, i will decostand() the data
# then i can merge all the data to put into ggplot and facet and such

year <- c(2020, 2021, 2022, 2023, 2024)
cov.stand <- NULL
for(i in year){
  cov.CLIFF_deco <- cov.CLIFF %>%
    filter(Year==i)
  deco <- cov.CLIFF_deco %>%
    select(-c(plot, Year, severity)) %>%
    decostand(method = "total")
  cov.stand<-rbind(cov.stand, deco)
}


# Stats: PermANOVA

traits <- traits %>% 
  mutate(spp = rownames(.)) %>% 
  mutate(gram = case_when(spp %in% c("CARO", "ELEL", "FEAR",
                                     "MUMO", "MUST", "PIPR", 
                                     "SCSC") ~ 2, .default = 1),
         shrub = case_when(spp == "CEFE" ~ 2, .default = 1),
         tree = case_when(spp == "QUGA" ~ 2, .default = 1),
         forb = case_when(!(spp %in% c("CARO", "ELEL", "FEAR",
                            "MUMO", "MUVI", "PIPR", 
                            "SCSC", "CEFE", "QUGA")) ~ 2, .default = 1),
         exotic = case_when(spp %in% c("VETH", "LIDA", "SATR") ~ 2,
                              .default = 1))

traits_exo <- traits %>% 
  select(exotic)

traits_type <- traits %>% 
  select(gram, shrub, tree, forb)

# Figures:
cov.stand$year <- cov.CLIFF$Year
cov.stand$severity <- cov.CLIFF$severity
cov.stand$plot <- cov.CLIFF$plot

group_cover <- cov.stand %>% 
  pivot_longer(cols = ARLU:VETH, names_to = "spp") %>% 
  mutate(group = case_when(spp %in% c("CARO", "ELEL", "FEAR",
                                      "MUMO", "MUVI", "PIPR", 
                                      "SCSC") ~ "Graminoid",
                           spp == "CEFE" ~ "Shrub",
                           spp == "QUGA" ~ "Tree",
                           .default = "Forb"),
         nativity = case_when(spp %in% c("VETH", "LIDA", "SATR") ~ "Exotic",
                              .default = "Native"))

# Functional groups:
group_cover2a <- group_cover %>% 
  group_by(severity, plot, group, year) %>% 
  summarise(fun_cov = sum(value)) %>% 
  ungroup()

fun_cover <- group_cover2a %>% 
  pivot_wider(names_from = group, values_from = fun_cov)

rel_fun_cover <- fun_cover %>% 
  select(-severity, -plot, -year) %>% 
  decostand(method = "total")

rel_fun_cover$severity <- fun_cover$severity
rel_fun_cover$year <- fun_cover$year
rel_fun_cover$plot <- fun_cover$plot

CTRL.t <- how(within = Within(type = "free"), #restrict permutations for repeat measures
              plots = Plots(type = "none"),
              blocks=fun_cover$plot,
              nperm = 999,
              observed = TRUE)

fun <- select(rel_fun_cover, Forb, Tree, Shrub, Graminoid)

adonis.out<-adonis2(fun~plot+severity*year, #treats time as split plot factor, plot as sample unit per Bakker 2024
                    data = rel_fun_cover, # CAN SOMEONE CONFIRM THAT THIS IS CORRECT?
                    method = "bray", # changed from euclidean to bray
                    permutations = CTRL.t,
                    by = "margin")

adonis.out #interaction is significant, create new factor for pairwise adonis :)

rel_fun_cover$sev.year<-paste(rel_fun_cover$severity,rel_fun_cover$year)

pairwise.out<-pairwise.adonis2(fun~sev.year, #treats time as split plot factor, plot as sample unit per Bakker 2024
                               data=rel_fun_cover,
                               method="bray",
                               by="margin")
pairwise.out #multiple significant pairwise comparisons with bonferonni (?) correction

write.csv(pairwise.out, "funPERMresults.csv", row.names = F)

cover_df <- rel_fun_cover %>% 
  pivot_longer(cols = Forb:Tree, names_to = "group") %>% 
  group_by(severity, group, year) %>% 
  summarise(cov = mean(value), cov_sd = sd(value))

cover_df$severity <- factor(cover_df$severity, c("U", "L", "H"))
severity_colors <- c("U" = "#5bbcd695", "L" = "#f9840295", "H" = "#fb040495")

ggplot(cover_df, aes(x = year, y = cov*100))+
  theme_light(base_size = 18)+
  geom_bar(stat = "identity", alpha = 1,
           color = "black", aes(group = interaction(year, severity), 
                                                   fill = severity),
           position = position_dodge())+
  geom_errorbar(aes(ymax = (cov*100 + cov_sd*100), ymin = cov*100, group = severity),
                position = position_dodge())+
  labs(x = "Year", y = "Relative cover (%)",
       fill = "Severity")+
  scale_fill_manual(values = severity_colors)+
  facet_wrap(~group, nrow = 4)+
  theme(strip.background = element_rect(color = "black", fill = "white"))+
  theme(strip.text = element_text(colour = 'black'))
# ggsave("outputs/5yr_fun_cover.png", last_plot(),
#        width = 5, height = 12, units = "in", dpi = 600)

# functional group stats:
# 
# TukeyHSD(aov(Forb ~ severity, data = rel_fun_cover)) # forbs
# TukeyHSD(aov(Shrub ~ severity, data = rel_fun_cover)) # shrubs
# TukeyHSD(aov(Graminoid ~ severity, data = rel_fun_cover)) # grammies
# TukeyHSD(aov(Tree ~ severity, data = rel_fun_cover)) # trees


# Nativity:
group_cover2b <- group_cover %>% 
  group_by(severity, plot, year, nativity) %>% 
  summarise(type_cov = sum(value)) %>% 
  ungroup()

nat_cover <- group_cover2b %>% 
  pivot_wider(names_from = nativity, values_from = type_cov)

rel_nat_cover <- nat_cover %>% 
  select(-severity, -plot, -year) %>% 
  decostand(method = "total")

rel_nat_cover$severity <- nat_cover$severity
rel_nat_cover$year <- nat_cover$year
rel_nat_cover$plot <- nat_cover$plot

type_df <- rel_nat_cover %>% 
  pivot_longer(cols = Exotic:Native, names_to = "group") %>% 
  group_by(severity, group, year) %>% 
  summarise(cov = mean(value), cov_sd = sd(value))
type_df$severity <- factor(type_df$severity, c("U", "L", "H"))

CTRL.t <- how(within = Within(type = "free"), #restrict permutations for repeat measures
              plots = Plots(type = "none"),
              blocks=nat_cover$plot,
              nperm = 999,
              observed = TRUE)

nat <- select(rel_nat_cover, Exotic, Native)

adonis.out<-adonis2(nat~plot+severity*year, #treats time as split plot factor, plot as sample unit per Bakker 2024
                    data = rel_nat_cover, # CAN SOMEONE CONFIRM THAT THIS IS CORRECT?
                    method = "bray", # changed from euclidean to bray
                    permutations = CTRL.t,
                    by = "margin")

adonis.out #interaction is significant, create new factor for pairwise adonis :)

rel_nat_cover$sev.year<-paste(rel_nat_cover$severity,rel_nat_cover$year)

pairwise.out<-pairwise.adonis2(nat~sev.year, #treats time as split plot factor, plot as sample unit per Bakker 2024
                               data=rel_nat_cover,
                               method="bray",
                               by="margin")
pairwise.out

ggplot(type_df, aes(x = year, y = cov*100))+
  theme_light(base_size = 18)+
  geom_bar(stat = "identity", alpha = 1,
           color = "black", aes(group = interaction(year, severity), 
                                fill = severity),
           position = position_dodge())+
  geom_errorbar(aes(ymax = (cov*100 + cov_sd*100), ymin = cov*100, group = severity),
                position = position_dodge())+
  labs(x = "Year", y = "Relative cover (%)",
       fill = "Severity")+
  scale_fill_manual(values = severity_colors)+
  facet_wrap(~group, nrow = 4)+
  theme(strip.background = element_rect(color = "black", fill = "white"))+
  theme(strip.text = element_text(colour = 'black'))
# ggsave("outputs/5yr_exo_cover.png", last_plot(),
#        width = 5, height = 7, units = "in", dpi = 600)


# get these pairwise comparisons done
# extract p values
# do bonferroni correction
# then add letters

TukeyHSD(aov(Native ~ year:severity, data = rel_nat_cover)) # native
TukeyHSD(aov(Exotic ~ severity, data = rel_nat_cover)) # exotic

TukeyHSD(aov(Graminoid ~ severity*year, data = rel_fun_cover)) # native
TukeyHSD(aov(Exotic ~ severity, data = rel_nat_cover)) # exotic

# 
# 
# 
# # exotics permanova
# # but i gotta figure out how to make this go within years and what not
# perm_exo <- adonis2(vegdist(cwm_exo, method="euclidean") ~ severity, permutations=9999)
# pair_exo <- pairwise.adonis(vegdist(cwm_exo, method="euclidean"), severity, perm=9999)
# perm_exo
# pair_exo
# 
# # Test of beta dispersion and post-hoc pairwise test of beta dispersion
# anova(betadisper(vegdist(cwm_exo, method="euclidean"), data$Severity, type="centroid"))
# TukeyHSD(betadisper(vegdist(cwm_exo, method="euclidean"), data$Severity, type="centroid"))
# 
# # functional group permanova
# perm_type <- adonis2(vegdist(cwm_type, method="euclidean") ~ data$Severity, permutations=9999)
# pair_type <- pairwise.adonis(vegdist(cwm_type, method="euclidean"), data$Severity, perm=9999)
# perm_type
# pair_type
# 
# # Test of beta dispersion and post-hoc pairwise test of beta dispersion
# anova(betadisper(vegdist(cwm_type, method="euclidean"), data$Severity, type="centroid"))
# TukeyHSD(betadisper(vegdist(cwm_type, method="euclidean"), data$Severity, type="centroid"))
# 
# # cow <- plot_grid(a,b, ncol=1, align = "v", axis="1")
# # cow
# # ggsave("outputs/cover_plot.png", last_plot(),
# #        width = 8.5, height = 6, units = "in", dpi = 300)
