library(tidyverse)
library(vegan)
library(FD)
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(ggfortify)
library(ggordiplots)
library(sjPlot)
library(pagedown)

################################################################################
# Ok so a little bit of a code garbage heap...                                 #
# NMDS plots with trait vectors at the bottom                                  #  
################################################################################

# unlogging!
data <- read_csv("data/CommunityMatrix.csv")
traits <- read_csv("data/TraitTable.csv")%>%
  #mutate(ldmc=log(ldmc)) %>%
  #mutate(height=log(height)) %>%
  #mutate(sla=log(sla)) %>%
  #mutate(seedmass=log(seedmass)) %>%
  filter(spp %in% colnames(data)) %>%
  select(spp, sla, height, resprouting, seedmass) %>%
  column_to_rownames(var="spp")

# binary 0,1 doesnt work for some reason in the dbFD function so recode the traits
traits$resprouting <- ifelse(traits$resprouting == 1, 2, 1)

comm <- data %>%
  column_to_rownames(var = "Plot") %>% 
  select(ARLU:VETH) %>%
  wisconsin()

cwm <- dbFD(traits, comm)$CWM
str(data)

################################################################################
# TAXONOMIC PERMANOVA ##########################################################
################################################################################

# NMDS stress
tax_nmds <- metaMDS(comm, distance="bray", k=2, trymax=100)
tax_nmds$stress
stressplot(tax_nmds)

# PERMANOVA and post-hoc pairwise PERMANOVA
perm_tax <- adonis2(vegdist(comm, method = "bray") ~ data$Severity, permutations = 9999)
pair_tax <- pairwise.adonis(vegdist(comm, method="bray"), data$Severity, perm=9999)


# Test of beta dispersion and post-hoc pairwise test of beta dispersion
anova(betadisper(vegdist(comm, method="bray"), data$Severity, type="centroid"))
TukeyHSD(betadisper(vegdist(comm, method="bray"), data$Severity, type="centroid"))

# Let's get to plotting
tax_scores <- as_tibble(scores(tax_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "Unburned",
    sites %in% c(21:40) ~ "Low",
    sites %in% c(41:60) ~ "High"
  ))
  
tax_centroids <- tax_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  group_by(severity) %>%
  dplyr::summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

tax_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point(shape=1, size=2) +
  geom_point(data=tax_centroids, aes(NMDS1, NMDS2), size=3, shape=19) +
  theme_classic()

################################################################################
# FUNCTIONAL PERMANOVA ##########################################################
################################################################################

# Test if SLA and LDMC are correlated
cor.test(traits$sla, traits$ldmc, method="spearman")
ggplot(traits, aes(sla, ldmc)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme_classic()

# NMDS stress
fun_nmds <- metaMDS(cwm, distance="euclidean", k=2, trymax=100)
fun_nmds$stress
stressplot(fun_nmds)

# PERMANOVA and post-hoc pairwise PERMANOVA
perm_fun <- adonis2(vegdist(cwm, method="euclidean") ~ data$Severity, permutations=9999)
pair_fun <- pairwise.adonis(vegdist(cwm, method="euclidean"), data$Severity, perm=9999)

# Test of beta dispersion and post-hoc pairwise test of beta dispersion
anova(betadisper(vegdist(cwm, method="euclidean"), data$Severity, type="centroid"))
TukeyHSD(betadisper(vegdist(cwm, method="euclidean"), data$Severity, type="centroid"))

# Let's get to plotting
fun_scores <- as_tibble(scores(fun_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "Unburned",
    sites %in% c(21:40) ~ "Low",
    sites %in% c(41:60) ~ "High"
  ))

fun_centroids <- fun_scores %>%
  group_by(severity) %>%
  summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

fun_scores %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point(shape=1, size=2) +
  geom_point(data=fun_centroids, aes(NMDS1, NMDS2), size=3, shape=19) +
  theme_classic()

################################################################################
# ENVFIT PLOTS #################################################################
################################################################################

# envfit alpha set to 0.9 so non-significant traits are shown

tax_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  gg_envfit(tax_nmds, env=cwm, groups=data$Severity, alpha=0.9)

fun_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  gg_envfit(fun_nmds, env=cwm, groups=data$Severity, alpha=0.9)



################################################################################
# ENVFIT + ELLIPSES PLOTS ######################################################
################################################################################

# functional nmds with ellipses AND arrows so help me god
png("outputs/fun_nmds.png", width = 7.5, height = 5, units = "in", res = 300)

fun_scores$severity <- factor(fun_scores$severity, levels = c("Unburned", "Low", "High"))
envfit_fun <- envfit(fun_nmds, cwm, perm=999)

print(fun_scores, n = 57)

# plot range to try and stretch points out
x_range <- range(c(fun_scores$NMDS1, fun_centroids$NMDS1, envfit_fun$vectors$arrows[,1]))
y_range <- range(c(fun_scores$NMDS2, fun_centroids$NMDS2, envfit_fun$vectors$arrows[,2]))

buffer_x <- diff(x_range) * 0.2
buffer_y <- diff(y_range) * 0.2

par(pin = c(5.5, 3))

# plot
plot(fun_nmds, display="sites", type="n")
     #xlim = c(min(x_range) - buffer_x, max(x_range) + buffer_x), 
     #ylim = c(min(y_range) - buffer_y, max(y_range) + buffer_y))

# colors
severity_colors_points <- c("Unburned" = "#5bbcd695", "Low" = "#f9840295", "High" = "#fb040495")
severity_colors_centroids <- c("Unburned" = "#5bbcd6", "Low" = "#f98402", "High" = "#fb0404")
severity_symbols <- c("Unburned" = 15, "Low" = 16, "High" = 17)
severity_lty <- c("Unburned" = 1, "Low" = 2, "High" = 3)

# points and centroids
points(fun_scores$NMDS1, fun_scores$NMDS2, 
       col = severity_colors_points[fun_scores$severity],
       bg = severity_colors_points[fun_scores$severity],
       pch = severity_symbols[fun_scores$severity],
       cex=1)

points(fun_centroids$NMDS1, fun_centroids$NMDS2, 
       col = severity_colors_centroids[fun_centroids$severity], 
       pch = 19, cex=1.5)

# add outline to centroids
points(fun_centroids$NMDS1, fun_centroids$NMDS2, 
       pch = 21,                        
       bg = severity_colors_centroids[fun_centroids$severity],  
       col = "black",                   
       lwd = 0.6,                       
       cex = 1.5)    

# ellipses
fun_ellipses <- ordiellipse(fun_nmds, 
            groups = fun_scores$severity, 
            kind = "sd", 
            conf = 0.95, 
            col = severity_colors_centroids[unique(fun_scores$severity)],
            lty = severity_lty[unique(fun_scores$severity)],
            lwd = 2,
            label = FALSE)
# legend
#legend("topleft", legend=names(severity_colors_centroids), 
       #col=severity_colors_centroids, pch=severity_symbols)

# add envfit arrows
plot(envfit_fun, p.max = 0.05, col = "black",
     xlim = c(min(x_range) - buffer_x, max(x_range) + buffer_x), 
     ylim = c(min(y_range) - buffer_y, max(y_range) + buffer_y))

# save
dev.off()


###############################################################################
# taxonomic nmds with ellipses AND arrows so help me god
png("outputs/tax_nmds.png", width = 7.5, height = 5, units = "in", res = 300)

tax_scores$severity <- factor(tax_scores$severity, levels = c("Unburned", "Low", "High"))
envfit_tax <- envfit(tax_nmds, cwm, perm=999)

# plot range to try and stretch points out
x_range <- range(c(tax_scores$NMDS1, tax_centroids$NMDS1, envfit_tax$vectors$arrows[,1]))
y_range <- range(c(tax_scores$NMDS2, tax_centroids$NMDS2, envfit_tax$vectors$arrows[,2]))

buffer_x <- diff(x_range) * 0.3
buffer_y <- diff(y_range) * 0.3

par(pin = c(5.5, 3))

# plot
plot(tax_nmds, display="sites", type="n",
     xlim = c(min(x_range) - buffer_x, max(x_range) + buffer_x), 
     ylim = c(min(y_range) - buffer_y, max(y_range) + buffer_y))

# colors
severity_colors_points <- c("Unburned" = "#5bbcd695", "Low" = "#f9840295", "High" = "#fb040495")
severity_colors_centroids <- c("Unburned" = "#5bbcd6", "Low" = "#f98402", "High" = "#fb0404")
severity_symbols <- c("Unburned" = 15, "Low" = 16, "High" = 17)
severity_lty <- c("Unburned" = 1, "Low" = 2, "High" = 3)

# points and centroids
points(tax_scores$NMDS1, tax_scores$NMDS2, 
       col = severity_colors_points[tax_scores$severity],
       bg = severity_colors_points[tax_scores$severity],
       pch = severity_symbols[tax_scores$severity],
       cex=1)

points(tax_centroids$NMDS1, tax_centroids$NMDS2, 
       col = severity_colors_centroids[tax_centroids$severity], 
       pch = 19, cex=1.5)

# add outline to centroids
points(tax_centroids$NMDS1, tax_centroids$NMDS2, 
       pch = 21,                        
       bg = severity_colors_centroids[tax_centroids$severity],  
       col = "black",                   
       lwd = 0.6,                       
       cex = 1.5)  

# ellipses
tax_ellipses <- ordiellipse(tax_nmds, 
                            groups = tax_scores$severity, 
                            kind = "sd", 
                            conf = 0.95, 
                            col = severity_colors_centroids[unique(tax_scores$severity)],
                            lty = severity_lty[unique(tax_scores$severity)],
                            lwd = 2,
                            label = FALSE)
# legend
legend("topleft", legend=names(severity_colors_centroids), 
       col=severity_colors_centroids, pch=severity_symbols)

# add envfit arrows
plot(envfit_tax, p.max = 1, col = "black",
     xlim = c(min(x_range) - buffer_x, max(x_range) + buffer_x), 
     ylim = c(min(y_range) - buffer_y, max(y_range) + buffer_y))

# save
dev.off()

