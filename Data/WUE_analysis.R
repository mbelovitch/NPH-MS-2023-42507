library(dplyr)

setwd("X:/Projects/WUE 2020/Final data analysis")
area <- read.csv('Cross_sectional_area.csv')
# Combine with biomass data
mass <- read.csv("Plant_mass_data.csv")
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
water <- read.csv('Water_addition.csv')
water$Transpiration_mass <- water$Initial_mass - water$Adj_final_mass + water$Water_added
water <- water[,c(1,7)]
wue <- merge(growth, water, by = 'UniqueID')
wue$WUE_shoot_unadj <- wue$deltaShoot / wue$Transpiration_mass
# Calculate whole-plant mass by using root:shoot ratios
# use 'bucket' data to obatin ratios for species used here
rmr <- read.csv('RMR.csv')
rmr.ag <- aggregate(RMR ~ Species, rmr, mean)
# Use Combretum collinum as a stand-in for Terminalia sericea
rmr.ag <- rmr.ag %>% mutate(Species = ifelse(Species == 'Comcol', 'Terser', Species))
wue <- merge(wue, rmr.ag)
wue$Root_init <- wue$Shoot_init * wue$RMR / (1 - wue$RMR)
wue$Root_final <- wue$Shoot_final * wue$RMR / (1 - wue$RMR)
wue$Mass_init <- wue$Shoot_init + wue$Root_init
wue$Mass_final <- wue$Shoot_final + wue$Root_final
wue$deltaMass <- wue$Mass_final - wue$Mass_init
wue$WUE_total_unadj <- wue$deltaMass / wue$Transpiration_mass
# Calculate error introduced by growth
# Total transpiration is reduced by net growth
# Add the seedling growth mass to the net mass loss (note: it is dry mass)
wue$Transpiration_adj <- wue$Transpiration_mass + wue$deltaMass
wue$WUE_shoot_adj <- wue$deltaShoot / wue$Transpiration_adj
wue$WUE_total_adj <- wue$deltaMass / wue$Transpiration_adj
# Adjusting for plant mass gain makes a negligible difference to calculations
ggplot(wue) + geom_boxplot(aes(x = FT, y = wue$WUE_shoot_unadj))
ggplot(wue) + geom_boxplot(aes(x = Species, y = wue$WUE_shoot_unadj, fill = FT))