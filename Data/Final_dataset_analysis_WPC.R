library(ggplot2)
library(lme4)
library(lmerTest)
library(nlme)
library(dplyr)
library(gridExtra)
library(cowplot)

setwd("D:/Ecophys PhD/Whole Plant Conductance (Ch.1)/Data")
trans <- read.csv('Compiled_hourly_transpiration_WP_weather_data.csv')

# Look at diurnal patterns in a subset of data
IDsam <- sample(unique(trans$ID), 30)
trans.sub <- subset(trans, Hour > 7 & Hour < 17 & trans$ID %in% IDsam)
ggplot(trans.sub, aes(x = WP, y = deltaM.Shoot, color = Hour)) + geom_point() +
  facet_wrap(~ ID)

ggplot(trans.sub, aes(x = VPD, y = deltaM, color = Hour)) + geom_point() +
  facet_wrap(~ ID)

# There is no obvious hour effect within this hourly range
# Obtain daily mean daily values of deltaM and the environmental values
# This is a first pass to test for confounding effects of radiation and VPD
trans.ag <- aggregate(cbind(VPD, Rad_weighted, WP, deltaM.Leaf, deltaM.Shoot) ~
                        ID + yday + Batch + Species, trans.sub, mean)
# Get rid of single negative transpiration value to allow log-transformation
trans.ag <- subset(trans.ag, deltaM.Leaf > 0)
trans.ag$yday <- as.factor(trans.ag$yday)
trans.ag$Batch <- as.factor(trans.ag$Batch)
# Models do not converge well with Batch random effect - leave out
# This is simply to test the assumption that WP drives the patterns we see
mod.RadVPD <- lmer(log(deltaM.Leaf) ~ Rad_weighted + VPD +
                  (1 | Species) + (1 | yday) + (1 | ID),
                  trans.ag, REML = FALSE)
mod.VPD <- lmer(log(deltaM.Leaf) ~ VPD +
                  (1 | Species) + (1 | yday) + (1 | ID),
                   trans.ag, REML = FALSE)
mod.Null <- lmer(log(deltaM.Leaf) ~ 1 + 
                  (1 | Species) + (1 | yday) + (1 | ID),
                  trans.ag, REML = FALSE)
mod.WP <- lmer(log(deltaM.Leaf) ~ log(-WP) + 
                  (1 | Species) + (1 | yday) + (1 | ID),
                  trans.ag, REML = FALSE)
AIC(mod.RadVPD, mod.VPD, mod.WP, mod.Null)
# Shows that only WP appears to influence transpiration rates in the chamber

# Fig. S1 ----------------------------------------------------------------------
# 4 panels showing daily transpiration traces (4 grasses, 4 trees)?
trace <- trans[,c(1,2,4,5,7,12,16,19)]
trace$time <- ((trace$Day*24 + trace$Hour))

trace.tree <- subset(trace, Group == 'Tree')
trace.tree <- subset(trace.tree, ID == c("Acaexu2","Domrot3","Diccin1","Berdis1"))
ggplot(trace.tree) +
  geom_line(aes(x = time, y = deltaM.Leaf)) +
  labs(x = "Hours in Chamber",
       y = expression(paste('E (g g'^'-1',' h'^'-1',')'))) +
  facet_wrap(~ ID) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"))
        
trace.grass <- subset(trace, Group == 'Grass')
trace.grass <- subset(trace.grass, ID == c("Trigra2","Erasup1","Arimer2","Perpat2"))
ggplot(trace.grass) +
  geom_line(aes(x = time, y = deltaM.Leaf)) +
  labs(x = "Hours in Chamber",
       y = expression(paste('E (g g'^'-1',' h'^'-1',')'))) +
  facet_wrap(~ ID) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"))

# Fig. S2
# Radiation fig. Would be good to show that radiation in the chamber was not a factor
# influencing transpiration
figS2 <- ggplot(trans.ag) + 
  geom_point(aes(x = Rad_weighted, y = log(deltaM.Leaf)), size = 4) +
  labs(x = expression(paste('Height-weighted PAR (',mu,'mol m'^'-2','s'^'-1',')')),
       y =  expression(paste('log E (g g'^'-1',' h'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.position = "none")
plotpath <- file.path("X:/Projects/Whole Plant Conductance 2019/Paper/Figs/FigS2.pdf")
pdf(file = plotpath)
figS2
dev.off()

# Calculate mean and max daily E, as well as the within-day change in E
trans.day <- subset(trans, Hour > 7 & Hour < 17)
IDs <- unique(trans.day$ID)
N <- length(IDs)
for (i in 1:N){
  sub1 <- subset(trans.day, ID == IDs[i])
  days <- unique(sub1$Day)
  K <- length(days)
  for (j in 1:K){
    sub2 <- subset(sub1, Day == days[j])
    # Exclude partial days (6 hours of data minimum)
    if (nrow(sub2) > 5){
      df <- data.frame(IDs[i])
      df$Day <- days[j]
      df$WP <- mean(sub2$WP)
      df$Emean <- mean(sub2$deltaM.Leaf, na.rm = TRUE)
      df$Emax <- max(sub2$deltaM.Leaf, na.rm = TRUE)
      reg <- lm(deltaM.Leaf ~ Hour, sub2)
      df$Slope <- coef(reg)[2]
      if (i == 1 & j == 1) trans.all <- df
      else trans.all <- rbind(trans.all, df)
    }
  }
}
names(trans.all)[1] <- 'ID'
# Create a dataframe binned by WP
trans.all <- trans.all %>% mutate(WPbin = case_when(WP > -0.5 ~ '< 0.5',
                WP > -1.0  & WP <= -0.5 ~ '0.5-1.0',
                WP > -1.5  & WP <= -1.0 ~ '1.0-1.5',
                WP <= -1.5 ~ '> 1.5'))
trans.bin1 <- aggregate(cbind(Emean, Emax, Slope) ~ ID + WPbin, trans.all, mean)

# Calculate daily change in E from evening to morning for each individual
# First, get rid of stuff we don't need
trans2 <- trans[,c(1,2,4,5,12,16,18,19)]
trans8am <- subset(trans2, Hour == 8)
trans4pm <- subset(trans2, Hour == 16)
trans4pm$Day <- trans4pm$Day + 1
delta <- merge(trans8am, trans4pm, by = c("ID", "Species", "Day", "Batch"))
delta$Ediff <- delta$deltaM.Leaf.x - delta$deltaM.Leaf.y
delta <- delta[,-c(6,8:11)]
names(delta)[6] <- 'WP'
delta <- delta %>% mutate(WPbin = case_when(WP > -0.5 ~ '< 0.5',
                WP > -1.0  & WP <= -0.5 ~ '0.5-1.0',
                WP > -1.5  & WP <= -1.0 ~ '1.0-1.5',
                WP <= -1.5 ~ '> 1.5'))
trans.bin2 <- aggregate(Ediff ~ ID + WPbin, delta, mean)

# Combine the two datasets
trans.final <- merge(trans.bin1, trans.bin2, by = c('ID', 'WPbin'), all = TRUE)
# Reincorporate metadata
meta <- delta[,c(1,2,4,5)]
meta <- unique(meta)
trans.final <- merge(meta, trans.final)
trans.final$WPbin <- factor(trans.final$WPbin,
    levels = c("< 0.5", "0.5-1.0", "1.0-1.5", "> 1.5" ))

# Fig 1: a single transpiration trace with the four daily metrics? -------------------
trace.fig1 <- subset(trace, ID == "Acaexu2")
trace.fig1$time <- (trace.fig1$Day*24 + trace.fig1$Hour)
fig1 <- ggplot(trace.fig1) +
  geom_line(aes(x = time, y = deltaM.Leaf)) +
  labs(x = "Hours in Chamber",
       y = expression(paste('E (g g'^'-1',' h'^'-1',')'))) +
  ylim(0, 3.25) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(2)),
        axis.title.x = element_text(vjust = 0, size = rel(2)),
        axis.text = element_text(size = rel(1.5)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm")) +
  annotate("text", x = 55, y = 2.78, label = "Emax", parse = FALSE, size = 4)+
  geom_segment(aes(x = 35, y = 2.78, xend = 50, yend = 2.78)) +
  annotate("text", x = 56.5, y = 2.66, label = "Emean", parse = FALSE, size = 4)+
  geom_segment(aes(x = 33.5, y = 2.66, xend = 50, yend = 2.66)) +
  geom_segment(aes(x = 60, y = 2.0, xend = 69, yend = 1.25)) +
  geom_segment(aes(x = 44, y = 2.35, xend = 57, yend = 1.9)) +
  annotate("text", x = 55.5, y = 2.25, label = "Ebd", parse = FALSE, size = 4) +
  annotate("text", x = 70, y = 1.8, label = "Ewd", parse = FALSE, size = 4)



# Plot daily mean transpiration by functional type and bin
fig2a <- ggplot(trans.final) +
  geom_boxplot(aes(x = WPbin, y = Emean, fill = Group)) + 
  labs(x = Psi[s]~'(-MPa)', 
       y =  expression(paste('E'['mean'],' (g g'^'-1',' h'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

fig2b <- ggplot(trans.final) +
  geom_boxplot(aes(x = WPbin, y = Emax, fill = Group)) + 
  labs(x = Psi[s]~'(-MPa)', 
       y =  expression(paste('E'['max'],' (g g'^'-1',' h'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

fig2c <- ggplot(trans.final) +
  geom_boxplot(aes(x = WPbin, y = Slope, fill = Group)) + 
  labs(x = Psi[s]~'(-MPa)', 
       y =  expression(paste('E'['wd'],' (g g'^'-1',' h'^'-2',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

fig2d <- ggplot(trans.final) +
  geom_boxplot(aes(x = WPbin, y = Ediff, fill = Group)) + 
  labs(x = Psi[s]~'(-MPa)', 
       y =  expression(paste('E'['bd'],' (g g'^'-1',' h'^'-1',')'))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1.2)),
        axis.title.x = element_text(vjust = 0, size = rel(1.2)),
        axis.text = element_text(size = rel(1)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)))

prow <- plot_grid(
  fig2a + theme(legend.position="none"),
  fig2b + theme(legend.position="none"),
  fig2c + theme(legend.position="none"),
  fig2d + theme(legend.position="none"),
  align = 'vh',
  labels = c('(a)','(b)','(c)','(d)'),
  hjust = 0,
  nrow = 2 )
legend <- get_legend(
  # Create some space to the left of the legend
  fig2a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

fig2 <- plot_grid(prow, legend, rel_widths = c(3, .4), scale = 1)
plotpath <- file.path("X:/Projects/Whole Plant Conductance 2019/Paper/Figs/Fig2.pdf")
pdf(file = plotpath)
fig2
dev.off()

# Stats
# Compare Emean values of WP > -0.5 Mpa
mod.Emean1 <- lmer(Emean ~ Group + (1 | Batch) + (1 | Species),
                   subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
mod.Emean1.null <- lmer(Emean ~ 1 + (1 | Batch) + (1 | Species),
                        subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Emean1, mod.Emean1.null)

# Compare Emean values of WP between -0.5 and -1.0  Mpa
mod.Emean2 <- lmer(Emean ~ Group + (1 | Batch) + (1 | Species),
                   subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
mod.Emean2.null <- lmer(Emean ~ 1 + (1 | Batch) + (1 | Species),
                        subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
anova(mod.Emean2, mod.Emean2.null)

# Compare Emean values of WP between -1.0 and -1.5  Mpa
mod.Emean3 <- lmer(Emean ~ Group + (1 | Batch) + (1 | Species),
                   subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
mod.Emean3.null <- lmer(Emean ~ 1 + (1 | Batch) + (1 | Species),
                        subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
anova(mod.Emean3, mod.Emean3.null)

# Compare Emean values of WP between < -1.5 Mpa
mod.Emean4 <- lmer(Emean ~ Group + (1 | Batch) + (1 | Species),
                   subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
mod.Emean4.null <- lmer(Emean ~ 1 + (1 | Batch) + (1 | Species),
                        subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
anova(mod.Emean4, mod.Emean4.null)

#-------------------------------------------------------------------------------

# Compare Emax values of WP > -0.5 Mpa
mod.Emax1 <- lmer(Emax ~ Group + (1 | Batch) + (1 | Species),
                   subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
mod.Emax1.null <- lmer(Emax ~ 1 + (1 | Batch) + (1 | Species),
                        subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Emax1, mod.Emax1.null)

# Compare Emax values of WP between -0.5 and -1.0  Mpa
mod.Emax2 <- lmer(Emax ~ Group + (1 | Batch) + (1 | Species),
                  subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
mod.Emax2.null <- lmer(Emax ~ 1 + (1 | Batch) + (1 | Species),
                       subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Emax2, mod.Emax2.null)

# Compare Emax values of WP between -1.0 and -1.5  Mpa
mod.Emax3 <- lmer(Emax ~ Group + (1 | Batch) + (1 | Species),
                  subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
mod.Emax3.null <- lmer(Emax ~ 1 + (1 | Batch) + (1 | Species),
                       subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Emax3, mod.Emax3.null)

# Compare Emax values of WP < -1.5 Mpa
mod.Emax4 <- lmer(Emax ~ Group + (1 | Batch) + (1 | Species),
                  subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
mod.Emax4.null <- lmer(Emax ~ 1 + (1 | Batch) + (1 | Species),
                       subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Emax4, mod.Emax4.null)

#-------------------------------------------------------------------------------

# Compare Eslope values of WP > -0.5 Mpa
mod.Eslope1 <- lmer(Slope ~ Group + (1 | Batch) + (1 | Species),
                  subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
mod.Eslope1.null <- lmer(Slope ~ 1 + (1 | Batch) + (1 | Species),
                       subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Eslope1, mod.Eslope1.null)

# Compare Eslope values of WP between -0.5 and -1.0  Mpa
mod.Eslope2 <- lmer(Slope ~ Group + (1 | Batch) + (1 | Species),
                    subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
mod.Eslope2.null <- lmer(Slope ~ 1 + (1 | Batch) + (1 | Species),
                         subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Eslope2, mod.Eslope2.null)

# Compare Eslope values of WP between -1.0 and -1.5  Mpa
mod.Eslope3 <- lmer(Slope ~ Group + (1 | Batch) + (1 | Species),
                    subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
mod.Eslope3.null <- lmer(Slope ~ 1 + (1 | Batch) + (1 | Species),
                         subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Eslope3, mod.Eslope3.null)

# Compare Eslope values of WP < -1.5 Mpa
mod.Eslope4 <- lmer(Slope ~ Group + (1 | Batch) + (1 | Species),
                    subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
mod.Eslope4.null <- lmer(Slope ~ 1 + (1 | Batch) + (1 | Species),
                         subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Eslope4, mod.Eslope4.null)

#-------------------------------------------------------------------------------

# Compare Ediff values of WP > -0.5 Mpa
mod.Ediff1 <- lmer(Ediff ~ Group + (1 | Batch) + (1 | Species),
                    subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
mod.Ediff1.null <- lmer(Ediff ~ 1 + (1 | Batch) + (1 | Species),
                         subset(trans.final, WPbin == '< 0.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Ediff1, mod.Ediff1.null)

# Compare Ediff values of WP between -0.5 and -1.0  Mpa
mod.Ediff2 <- lmer(Ediff ~ Group + (1 | Batch) + (1 | Species),
                    subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
mod.Ediff2.null <- lmer(Ediff ~ 1 + (1 | Batch) + (1 | Species),
                         subset(trans.final, WPbin == '0.5-1.0'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Ediff2, mod.Ediff2.null)

# Compare Ediff values of WP between -1.0 and -1.5  Mpa
mod.Ediff3 <- lmer(Ediff ~ Group + (1 | Batch) + (1 | Species),
                    subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
mod.Ediff3.null <- lmer(Ediff ~ 1 + (1 | Batch) + (1 | Species),
                         subset(trans.final, WPbin == '1.0-1.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Ediff3, mod.Ediff3.null)

# Compare Ediff values of WP < -1.5 Mpa
mod.Ediff4 <- lmer(Ediff ~ Group + (1 | Batch) + (1 | Species),
                    subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
mod.Ediff4.null <- lmer(Ediff ~ 1 + (1 | Batch) + (1 | Species),
                         subset(trans.final, WPbin == '> 1.5'), REML = FALSE)
# Likelihood ratio test (Group model vs. null model)
anova(mod.Ediff4, mod.Ediff4.null)

#-------------------------------------------------------------------------------

# Analysis of WRF turgor loss point data
osm <- read.csv('WRF TLP data/WRF_tlp.csv')
osm$Group2 <- ifelse(osm$Group == 'Grass', 'Grass', 'Tree')
# Aggregate by species for plotting purposes
osm.ag <- aggregate(psi_s_MPa ~ Code + Group2, osm, mean)
tlp.mod <- lme(psi_s_MPa ~ Group2, osm, random = ~ 1 | Code)
anova(tlp.mod)

# Analysis of WUE
area <- read.csv('WUE data/Cross_sectional_area.csv')
# Combine with biomass data
mass <- read.csv("WUE data/Plant_mass_data.csv")
# Photo mass equals leaf mass for trees, leaf + stem mass for grasses
# Compare two regressions for grasses: area vs. leaf and area vs. leaf + stem (shoot)
mass$Shoot_mass <- mass$Stem_mass + mass$Leaf_mass
final.area <- subset(area, Period == 'Final')
mass <- merge(mass, final.area, by = 'UniqueID')
mass$logleaf <- log(mass$Leaf_mass)
mass$logshoot <- log(mass$Shoot_mass)
mass$logarea <- log(mass$Area)
grass <- subset(mass, FT == 'Grass')
grass.mod1 <- lm(logleaf ~ logarea, grass)
summary(grass.mod1)
grass.mod2 <- lm(logshoot ~ logarea, grass)
summary(grass.mod2)
# The shoot regression has a higher R-squared
tree <- subset(mass, FT == 'Tree' & Species != "Comapi")
tree.mod1 <- lm(logleaf ~ logarea, tree)
summary(tree.mod1)
# Apply regressions to initial data to impute leaf and shoot biomass
init.area <- subset(area, Period == 'Initial')
init.area$logarea <- log(init.area$Area)
FT <- mass[,1:3]
init.area <- merge(init.area, FT, by = 'UniqueID')
init.grass <- subset(init.area, FT == 'Grass')
init.tree <- subset(init.area, FT == 'Tree')
init.grass$Shoot_mass <- exp(predict(grass.mod2, newdata = init.grass))
init.tree$Leaf_mass <- exp(predict(tree.mod1, newdata = init.tree))
# Use final leaf:stem ratio to infer initial tree shoot mass
tree$leaf_stem <- tree$Leaf_mass / tree$Stem_mass
treels <- tree[,c(1,16)]
init.tree <- merge(init.tree, treels)
init.tree$Stem_mass <- init.tree$Leaf_mass / init.tree$leaf_stem
init.tree$Shoot_mass <- init.tree$Leaf_mass + init.tree$Stem_mass
tree <- tree[,c(1:3,10)]
init.tree <- init.tree[,c(1,10)]
tree.growth <- merge(tree, init.tree, by = 'UniqueID')
names(tree.growth)[4:5] <- c('Shoot_final','Shoot_init')
# Now calculate grass growth
# Need to add on shoot mass removed from 3 grasses before photos
grass$Shoot_mass <- ifelse(!is.na(grass$Excluded_mass),
                           grass$Shoot_mass + grass$Excluded_mass, grass$Shoot_mass)
grass <- grass[,c(1:3,10)]
init.grass <- init.grass[,c(1,7)]
grass.growth <- merge(grass, init.grass, by = 'UniqueID')
names(grass.growth)[4:5] <- c('Shoot_final','Shoot_init')
growth <- rbind(tree.growth,grass.growth)
growth$deltaShoot <- growth$Shoot_final - growth$Shoot_init

# Now calculate water use
water <- read.csv('WUE data/Water_addition.csv')
water$Transpiration_mass <- water$Initial_mass - water$Adj_final_mass + water$Water_added
water <- water[,c(1,7)]
wue <- merge(growth, water, by = 'UniqueID')
wue$WUE_shoot_unadj <- wue$deltaShoot / wue$Transpiration_mass
# Calculate whole-plant mass by using root:shoot ratios
# Use to estimate error introduced by ignoring plant mass change
# use 'bucket' data to obtain ratios for species used here
rmr <- read.csv('WUE data/RMR.csv')
rmr.ag <- aggregate(RMR ~ Species, rmr, mean)
# Use Combretum collinum as a stand-in for Terminalia sericea
rmr.ag <- rmr.ag %>% mutate(Species = ifelse(Species == 'Comcol', 'Terser', Species))
wue <- merge(wue, rmr.ag)
wue$Root_init <- wue$Shoot_init * wue$RMR / (1 - wue$RMR)
wue$Root_final <- wue$Shoot_final * wue$RMR / (1 - wue$RMR)
wue$Mass_init <- wue$Shoot_init + wue$Root_init
wue$Mass_final <- wue$Shoot_final + wue$Root_final
wue$deltaMass <- wue$Mass_final - wue$Mass_init
# Calculate error introduced by growth
# Total transpiration is reduced by net growth
# Add the seedling growth mass to the net mass loss (note: it is dry mass)
wue$Transpiration_adj <- wue$Transpiration_mass + wue$deltaMass
wue$WUE_shoot_adj <- wue$deltaShoot / wue$Transpiration_adj
# Adjusting for plant mass gain makes a negligible difference to calculations
# Change 'FT' to 'Group' for consustency with other datasets
names(wue)[3] <- 'Group'
wue.ag <- aggregate(WUE_shoot_adj ~ Species + Group, wue, mean)
wue.mod.shoot <- lme(WUE_shoot_adj ~ Group, wue, random = ~ 1 | Species)
anova(wue.mod.shoot)

# RGR analysis
con <- read.csv('RGR data/Dry_wet_mass_ratios.csv')
con$Ratio <- con$Tmass_dry / con$Tmass_wet
con.ag <- aggregate(Ratio ~ Functional_Type, con, mean)
con.ag.sd <- aggregate(Ratio ~ Functional_Type, con, sd)
rgr <- read.csv('RGR data/RGR_data.csv')
# Calculate growth time (Harvest - Transplantation)
rgr <- merge(rgr, con.ag)
Dtrans <- strptime(rgr$Trans_date, format = '%m/%d/%Y')
Dharv <- strptime(rgr$Harv_date, format = '%m/%d/%Y')
rgr$Minit <- rgr$Mass_init_wet * rgr$Ratio
rgr$Mfinal <- rgr$AG_drymass + rgr$BG_drymass
rgr$deltaT <- Dharv$yday - Dtrans$yday
rgr$RGR <- (log(rgr$Mfinal) - log(rgr$Minit)) / rgr$deltaT
names(rgr)[1] <- 'Group'
rgr.ag <- aggregate(RGR ~ Group + Species, rgr, mean)
rgr.mod <- lme(RGR ~ Group, rgr, random = ~ 1 | Species)
anova(rgr.mod)

# Plot tlp, WUE and RGR together
fig3a <- ggplot(data = osm.ag, aes(x = Group2, y = psi_s_MPa)) +
  geom_boxplot(outlier.size = 3) +
  labs(x = '', 
       y =  expression(paste(pi[0],' (MPa)'))) +
  ylim(-5,0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8)))
fig3a <- fig3a + annotate("text", x = 1.5, y = 0, label = "P < 0.005", parse = FALSE, size = 4)

fig3b <- ggplot(data = wue.ag, aes(x = Group, y = WUE_shoot_adj)) + 
  geom_boxplot(outlier.size = 3) +
  labs(x = '', 
       y =  expression(paste('WUE (g g'^'-1',')'))) +
  ylim(0.002,0.008) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8)))
fig3b <- fig3b + annotate("text", x = 1.5, y = 0.008, label = "P < 0.005", parse = FALSE, size = 4)

fig3c <- ggplot(data = rgr.ag, aes(x = Group, y = RGR)) + 
  geom_boxplot(outlier.size = 3) +
  labs(x = '', 
       y =  expression(paste('RGR (log g d'^'-1',')'))) +
  ylim(0,0.25) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(0.8)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8)))
fig3c <- fig3c + annotate("text", x = 1.5, y = 0.25, label = "P < 0.0001", parse = FALSE, size = 4)

fig3 <- plot_grid(fig3a, fig3b, fig3c,
  align = 'v',
  labels = c('(a)','(b)','(c)'),
  hjust = 0,
  nrow = 1
)
plotpath <- file.path("X:/Projects/Whole Plant Conductance 2019/Paper/Figs/Fig3.pdf")
pdf(file = plotpath, height = 3)
fig3
dev.off()