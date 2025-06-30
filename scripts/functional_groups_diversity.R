library(tidyverse)
library(vegan)
library(FD)
library(wesanderson)
library(cowplot)
library(pairwiseAdonis)


data <- read_csv("data/CommunityMatrix.csv")
comm <- data %>%
  column_to_rownames(var = "Plot") %>% 
  select(ARLU:VETH) %>%
  wisconsin()

# Stats: PermANOVA

traits <- read_csv("data/TraitTable.csv")%>%
  filter(spp %in% colnames(data)) %>%
  mutate(gram = case_when(spp %in% c("CARO", "ELEL", "FEAR",
                                     "MUMO", "MUVI", "PIPR", 
                                     "SCSC") ~ 2, .default = 1),
         shrub = case_when(spp == "CEFE" ~ 2, .default = 1),
         tree = case_when(spp == "QUGA" ~ 2, .default = 1),
         forb = case_when(!(spp %in% c("CARO", "ELEL", "FEAR",
                            "MUMO", "MUVI", "PIPR", 
                            "SCSC", "CEFE", "QUGA")) ~ 2, .default = 1),
         exotic = case_when(spp %in% c("VETH", "LIDA", "SATR") ~ 2,
                              .default = 1))

traits_exo <- traits %>% 
  select(spp, exotic) %>% 
  column_to_rownames(var="spp")

traits_type <- traits %>% 
  select(spp, gram, shrub, tree, forb) %>% 
  column_to_rownames(var="spp")

cwm_exo <- dbFD(traits_exo, comm)$CWM
cwm_type <- dbFD(traits_type, comm)$CWM


perm_exo <- adonis2(vegdist(cwm_exo, method="euclidean") ~ data$Severity, permutations=9999)
pair_exo <- pairwise.adonis(vegdist(cwm_exo, method="euclidean"), data$Severity, perm=9999)
perm_exo
pair_exo

# Test of beta dispersion and post-hoc pairwise test of beta dispersion
anova(betadisper(vegdist(cwm_exo, method="euclidean"), data$Severity, type="centroid"))
TukeyHSD(betadisper(vegdist(cwm_exo, method="euclidean"), data$Severity, type="centroid"))


perm_type <- adonis2(vegdist(cwm_type, method="euclidean") ~ data$Severity, permutations=9999)
pair_type <- pairwise.adonis(vegdist(cwm_type, method="euclidean"), data$Severity, perm=9999)
perm_type
pair_type

# Test of beta dispersion and post-hoc pairwise test of beta dispersion
anova(betadisper(vegdist(cwm_type, method="euclidean"), data$Severity, type="centroid"))
TukeyHSD(betadisper(vegdist(cwm_type, method="euclidean"), data$Severity, type="centroid"))

# Figures:
group_cover <- data %>% 
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
  group_by(Severity, Plot, group) %>% 
  summarise(fun_cov = sum(value)) %>% 
  ungroup()

fun_cover <- group_cover2a %>% 
  pivot_wider(names_from = group, values_from = fun_cov)

rel_fun_cover <- fun_cover %>% 
  select(-Severity, -Plot) %>% 
  decostand(method = "total")

rel_fun_cover$severity <- fun_cover$Severity

cover_df <- rel_fun_cover %>% 
  pivot_longer(cols = Forb:Tree, names_to = "group") %>% 
  group_by(severity, group) %>% 
  summarise(cov = mean(value), cov_sd = sd(value))

cover_df$severity <- factor(cover_df$severity, c("U", "L", "H"))
  
ggplot(cover_df, aes(x = severity, y = cov*100, fill = group))+
  theme_light(base_size = 18)+
  geom_bar(stat = "identity", color = "black")+
  geom_errorbar(aes(ymax = (cov*100 + cov_sd*100), ymin = cov*100),
                width = 0.1)+
  scale_y_continuous(breaks = c(25, 50, 75, 100),
                   labels = c("25", "50", "75", "100"),
                   limits = c(0, 110))+
  labs(x = "Burn severity", y = "Relative Cover (%)",
       fill = "Functional type")+
  scale_fill_manual(values = wes_palette("Cavalcanti1", 4))+
  facet_wrap(~group)+
  theme(strip.background = element_rect(color = "black", fill = "white"))+
  theme(strip.text = element_text(colour = 'black'))+
  theme(legend.position = "none")

# ggsave("outputs/functional_group_cover_facet.png", last_plot(),
#        width = 8, height = 8, units = "in")

# functional group stats:

TukeyHSD(aov(Forb ~ severity, data = rel_fun_cover)) # forbs
TukeyHSD(aov(Shrub ~ severity, data = rel_fun_cover)) # shrubs
TukeyHSD(aov(Graminoid ~ severity, data = rel_fun_cover)) # grammies
TukeyHSD(aov(Tree ~ severity, data = rel_fun_cover)) # trees


# Nativity:
group_cover2b <- group_cover %>% 
  group_by(Severity, Plot, nativity) %>% 
  summarise(type_cov = sum(value)) %>% 
  ungroup()

nat_cover <- group_cover2b %>% 
  pivot_wider(names_from = nativity, values_from = type_cov)

rel_nat_cover <- nat_cover %>% 
  select(-Severity, -Plot) %>% 
  decostand(method = "total")

rel_nat_cover$severity <- nat_cover$Severity

type_df <- rel_nat_cover %>% 
  pivot_longer(cols = Exotic:Native, names_to = "group") %>% 
  group_by(severity, group) %>% 
  summarise(cov = mean(value), cov_sd = sd(value))
type_df$severity <- factor(type_df$severity, c("U", "L", "H"))

ggplot(type_df, aes(x = severity, y = cov*100, fill = group))+
  theme_light(base_size = 18)+
  geom_bar(stat = "identity", color = "black")+
  geom_errorbar(aes(ymax = (cov*100 + cov_sd*100), ymin = cov*100), width = 0.1)+
  scale_y_continuous(breaks = c(25, 50, 75, 100),
                     labels = c("25", "50", "75", "100"),
                     limits = c(0, 110))+
  labs(x = "Burn severity", y = "Relative Cover (%)",
       fill = "Nativity")+
  scale_fill_manual(values = c("#ae94a3", "#4d795f"))+
  facet_wrap(~group)+
  theme(strip.background = element_rect(color = "black", fill = "white"))+
  theme(strip.text = element_text(colour = 'black'))+
  theme(legend.position = "none")
# ggsave("outputs/nativity_cover_facet.png", last_plot(),
#        width = 8, height = 4, units = "in")

TukeyHSD(aov(Native ~ severity, data = rel_nat_cover)) # native
TukeyHSD(aov(Exotic ~ severity, data = rel_nat_cover)) # exotic

# cow <- plot_grid(a,b, ncol=1, align = "v", axis="1")
# cow
# ggsave("outputs/cover_plot.png", last_plot(),
#        width = 8.5, height = 6, units = "in", dpi = 300)
