# Sterling Herron
# Chapter 2
# Last updated: Dec 2020
# 2017 data (St Louis annual vs. perennial comparative analysis)
# cleaned dataframes to use: comb.seed.clean, comb.acc.clean, comb.germ.clean, acc.germ.clean, comb.plant.clean.new
# Version info & packages ####
version$version.string # R version
RStudio.Version()$version # R studio version

library(readxl)
library(plyr) # have to load plyr before dplyr (subsequent packages require dplyr), since if plyr loads after dplyr it messes up some dplyr functions
library(tidyverse)
library(car)
library(xlsx)
library(Hmisc)
library(corrplot)
library(ggpubr)
library(XLConnect)
library(lme4)
library(lmerTest)
library(emmeans)
library(nlme)
library(afex)
library(viridis)
library(igraph)
library(factoextra)
library(missMDA)
library(RColorBrewer)
library(devtools)
library(ggbiplot) # have to have devtools loaded and install via install_github("vqv/ggbiplot")
library(psych)

# article on grid plots
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/#use-common-legend-for-combined-ggplots

# File upload ####

# ACCESSION-LEVEL BASE FOR ADDING DATA TO
comb.acc <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/combined_accession_individual_data2017.xlsx", 
                   sheet = "accession_base", na = "NA")
View(comb.acc)

# INDIVIDUAL SEED DATA
comb.seed <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/complete_seed_data2017.xlsx", 
                       sheet = "new_data", na = "NA")
View(comb.seed)

# GERMINATION ORIGINAL TOTALS DATA
orig.germ <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/germination_rep_info2017.xlsx", 
                        sheet = "totals_simplified", na = "NA")
View(orig.germ)

# GERMINATION 2ND PHASE DATA
phase.2.germ <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/germination_rep_info2017.xlsx", 
                        sheet = "second_phase", na = "NA")
View(phase.2.germ)

# GERMINATION REPLICATE DATA
comb.germ <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/combined_accession_individual_data2017.xlsx", 
                        sheet = "replicate_germ", na = "NA")
View(comb.germ)

# INDIVIDUAL PLANT DATA
comb.plant <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/combined_accession_individual_data2017.xlsx", 
                      sheet = "individual_veg", na = "NA")
View(comb.plant)

# INDIVIDUAL LEAF DATA
comb.leaf <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/combined_accession_individual_data2017.xlsx", 
                         sheet = "leaf_data", na = "NA")
View(comb.leaf)

# Germination initial summaries & merging ####

totals.germ <- orig.germ %>%
  group_by(genus, lifespan, species.short, accession) %>%
  summarise(
    total.original.sum = sum(total.original.seeds, na.rm=T),
    germ.start = mean(germ.start, na.rm=T),
    main.total.lost.sum = sum(main.total.lost, na.rm=T),
    dead.scar.sum = sum(dead.scar, na.rm=T),
    germ.scar.sum = sum(germ.scar, na.rm=T),
    germ.dead.sum = sum(germ.dead, na.rm=T),
    total.nonviable.sum = sum(total.nonviable.simple, na.rm=T),
    seeds.imbibed.1.sum = sum(seeds.imbibed.1, na.rm=T),
    seeds.imbibed.2.sum = sum(seeds.imbibed.2, na.rm=T),
    seeds.unimbibed.sum = sum(seeds.unimbibed, na.rm=T),
    germ.1st.sum = sum(germ.1st, na.rm=T),
    time.1st.mn = mean(time.1st, na.rm=T),
    germ.2nd.sum = sum(germ.2nd, na.rm=T),
    time.2nd.mn = mean(time.2nd, na.rm=T),
    germ.both.sum = sum(germ.both, na.rm=T),
    time.both.mn = mean(time.both, na.rm=T),
    germ.total.sum = sum(germ.total.count, na.rm=T),
    germ.total.mn = mean(germ.total.count)
  ) %>% ungroup ()
View(totals.germ)
write.xlsx(totals.germ, file = "C:/Users/SterlingH/Desktop/R_exports/germ_sum2017.xlsx")

# NEW GERMINATION merge to get a new summary after each new update
# All accessions should be present in this dataframe, even if they only have one type of data.
germ.sum <- Reduce(function(x,y) merge(x,y, all=T, 
                                       by=c("genus","species.short","accession","lifespan")), 
                   list(totals.germ, phase.2.germ))
write.xlsx(germ.sum, file = "C:/Users/SterlingH/Desktop/R_exports/germ_check.xlsx") # check all needed data is present
# there should be 99 rows for 99 accessions (even though some will not be used)

# add summaries for germination totals
germ.sum.new <- germ.sum %>%
  replace_na(list(total.original.sum=0,main.total.lost=0,dead.scar.sum=0,germ.scar.sum=0,germ.dead.sum=0,total.nonviable.sum=0,
                  seeds.imbibed.1.sum=0,seeds.imbibed.2.sum=0,seeds.unimbibed.sum=0,germ.1st.sum=0,germ.2nd.sum=0,
                  germ.both.sum=0,p2.unimbibed=0,p2.germinated=0,p2.nonviable=0,p2.scar.dead=0, 
                  p2.scar.germ=0,p2.ungerminated=0,p2.total.lost=0,p2.missing=0,germ.between=0,germ.viability=0,p2.imb.other=0,
                  tz.viable=0,tz.nonviable=0,tz.noembryo=0,p2.germ.sum=0,p2.nonviable.sum=0,total.imb.left=0)) %>% # have to replace NAs in each column with 0 or else the mutate function outputs NAs
  rowwise() %>% # keeps calculations by row instead of all rows
  mutate(FINAL.total.sum = total.original.sum-main.total.lost.sum-dead.scar.sum-p2.scar.dead-p2.total.lost-p2.missing,
         FINAL.germ.sum = sum(germ.both.sum,p2.germ.sum),
         FINAL.viable.sum = sum(FINAL.germ.sum,tz.viable),
         FINAL.germ.prop = FINAL.germ.sum/FINAL.total.sum,
         FINAL.viable.prop = FINAL.viable.sum/FINAL.total.sum,
         check.1 = sum(dead.scar.sum, germ.scar.sum, total.nonviable.sum, seeds.imbibed.1.sum, seeds.imbibed.2.sum, seeds.unimbibed.sum, germ.both.sum, main.total.lost.sum), # do the numbers from the first round of germination add up to the original total?
         # not tabulating germ.dead.sum since redundant with germ.both.sum
         check.2 = sum(p2.germinated, p2.nonviable, p2.scar.dead, p2.ungerminated), # does the number counted in the second phase add up to the original unimbibed in this phase?
         imb.ungerm.original = sum(seeds.imbibed.1.sum, seeds.imbibed.2.sum, p2.ungerminated)) # to use to see if the number of germinants can be found by subtracting imb.ungerm.viability 
write.xlsx(germ.sum.new, file = "C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/germ_summary.xlsx")

# DATA FILTERING for individual trait datasets ####

# TO REMOVE all accessions (except L. sativus) with serious problems:
bad.acc <- c("PI 255365", "PI 260405", "PI 557540", "PI 638794", "W6 18752", "PI 494138", "PI 661844", "W6 20158", "PI 628322", "PI 602377", "PI 346067", "PI 422499", "W6 12942", "W6 13157", "W6 4872", "W6 26292", "PI 442562")
# Some of these accessions were already manually removed from some data files, but
# excluding again here because not all files had these excluded.

# TO REMOVE L. sativus
bad.species <- c("L. sativus")

# TO REMOVE problematic seeds from ImageJ file:
bad.seeds <- c("VC", "VU", "O") # O = misoriented; VC = very cracked; VU = very unusual (as yet nonexistent, but may be used later)

# FOR SEED mass data
comb.acc.clean1 <- comb.acc[!comb.acc$accession %in% bad.acc,]
comb.acc.clean <- comb.acc.clean1[!comb.acc.clean1$species.short %in% bad.species,]

# FOR SEED size raw data
comb.seed.clean1 <- comb.seed[!comb.seed$accession %in% bad.acc,]
comb.seed.clean2 <- comb.seed.clean1[!comb.seed.clean1$species.short %in% bad.species,]
comb.seed.clean <- comb.seed.clean2[!comb.seed.clean2$condition %in% bad.seeds,]
write.xlsx(comb.seed.clean, file = "C:/Users/SterlingH/Desktop/R_exports/seed_check.xlsx")

# FOR GERMINATION raw replicate data for rate (will be in accession data for figures, etc in the future)
# the germination-specific accessions and reps that needed to be excluded
# were excluded manually in the comb_acc_ind_all spreadsheet
comb.germ.clean1 <- comb.germ[!comb.germ$accession %in% bad.acc,]
comb.germ.clean <- comb.germ.clean1[!comb.germ.clean1$species.short %in% bad.species,]
write.xlsx(comb.germ.clean, file = "C:/Users/SterlingH/Desktop/R_exports/germt50_check.xlsx")

# FOR GERMINATION accession data (for proportions)
bad.germ <- c("W6 2747", "PI 494749") # remove accessions with unmet germination needs
acc.germ.clean1 <- germ.sum.new[!germ.sum.new$accession %in% bad.acc,]
acc.germ.clean2 <- acc.germ.clean1[!acc.germ.clean1$accession %in% bad.germ,]
acc.germ.clean <- acc.germ.clean2[!acc.germ.clean2$species.short %in% bad.species,]

# FOR GROWTH raw data
comb.plant.clean1 <- comb.plant[!comb.plant$accession %in% bad.acc,]
comb.plant.clean2 <- comb.plant.clean1[!comb.plant.clean1$species.short %in% bad.species,]
# Growth special filters
# Exclude row if main stem was damaged by day 21
dam <- c("A","B","C")
comb.plant.clean.new <- comb.plant.clean2[!comb.plant.clean2$cut.1 %in% dam,]
# Make day 35 height & nodes, & height and node growth NA if cut.2=A (but keep day 21 data, since cut.1 bad rows were already excluded)
comb.plant.clean.new$height.35.new[comb.plant.clean.new$cut.2=="A"] <- NA
comb.plant.clean.new$leaves.35.new[comb.plant.clean.new$cut.2=="A"] <- NA
comb.plant.clean.new$height.perday.21.35[comb.plant.clean.new$cut.2=="A"] <- NA
comb.plant.clean.new$leaves.perday.21.35[comb.plant.clean.new$cut.2=="A"] <- NA
comb.plant.clean.new$height.perday.21.35.ln[comb.plant.clean.new$cut.2=="A"] <- NA
comb.plant.clean.new$leaves.perday.21.35.ln[comb.plant.clean.new$cut.2=="A"] <- NA

# Exclude PI 358829 outliers for just these traits
comb.plant.clean.new$height.21.new[comb.plant.clean.new$accession=="PI 358829"] <- NA 
comb.plant.clean.new$height.perday.21.35[comb.plant.clean.new$accession=="PI 358829"] <- NA 
comb.plant.clean.new$leaves.perday.21.35[comb.plant.clean.new$accession=="PI 358829"] <- NA 
comb.plant.clean.new$height.perday.21.35.ln[comb.plant.clean.new$accession=="PI 358829"] <- NA 
comb.plant.clean.new$leaves.perday.21.35.ln[comb.plant.clean.new$accession=="PI 358829"] <- NA 

# Exclude PI 250758 for DAP-21 data and growth rate due to having only one plant
comb.plant.clean.new$height.21.new[comb.plant.clean.new$accession=="PI 250758"] <- NA 
comb.plant.clean.new$leaves.21.new[comb.plant.clean.new$accession=="PI 250758"] <- NA 
comb.plant.clean.new$height.perday.21.35[comb.plant.clean.new$accession=="PI 250758"] <- NA 
comb.plant.clean.new$leaves.perday.21.35[comb.plant.clean.new$accession=="PI 250758"] <- NA 
comb.plant.clean.new$height.perday.21.35.ln[comb.plant.clean.new$accession=="PI 250758"] <- NA 
comb.plant.clean.new$leaves.perday.21.35.ln[comb.plant.clean.new$accession=="PI 250758"] <- NA 

write.xlsx(comb.plant.clean.new, file = "C:/Users/SterlingH/Desktop/R_exports/plant_check.xlsx")
# Exclude problematic accessions
# bad.growth.acc <- c("W6 17187") # KEEPING THIS ACCESSION IN FOR NOW
# comb.plant.clean.new <- comb.plant.clean.new[!comb.plant.clean.new$accession %in% bad.growth.acc,]

# SUMMARIES for analysis & merging ####

# For accession-level seed size data:
acc.seed.sum <- comb.seed.clean %>%
  group_by(genus,lifespan,species.short,accession) %>%
  summarise(
    seed.n = length(seed.length[!is.na(seed.length)]),
    seed.width.mn = mean(seed.width, na.rm=T),
    seed.width.sd = sd(seed.width, na.rm=T),
    seed.length.mn = mean(seed.length, na.rm=T),
    seed.length.sd = sd(seed.length, na.rm=T),
    seed.perim.mn = mean(seed.perim, na.rm=T),
    seed.perim.sd = sd(seed.perim, na.rm=T),
    seed.area.mn = mean(seed.area, na.rm=T),
    seed.area.sd = sd(seed.area, na.rm=T),
    seed.circ.mn = mean(seed.circ, na.rm=T),
    seed.circ.sd = sd(seed.circ, na.rm=T),
    seed.round.mn = mean(seed.roundness, na.rm=T),
    seed.round.sd = sd(seed.roundness, na.rm=T)
  ) %>% ungroup ()
View(acc.seed.sum)

write.xlsx(acc.seed.sum, file = "C:/Users/SterlingH/Desktop/R_exports/seednum.xlsx")

# For accession-level germination rate data:
acc.germrep.sum <- comb.germ.clean %>%
  group_by(genus,lifespan,species.short,accession) %>%
  summarise(
    ppd50.mn = mean(ppd50, na.rm=T),
    ppd50.sd = sd(ppd50, na.rm=T)
  ) %>% ungroup ()
View(acc.germrep.sum)

# For accession-level growth data:
acc.growth.sum <- comb.plant.clean.new %>%
  group_by(genus,lifespan,species.short,accession) %>%
  summarise(
    height.perday.mn = mean(height.perday.21.35, na.rm=T),
    height.perday.sd = sd(height.perday.21.35, na.rm=T),
    height.perday.n = length(height.perday.21.35[!is.na(height.perday.21.35)]),
    nodes.perday.mn = mean(nodes.perday.21.35, na.rm=T),
    nodes.perday.sd = sd(nodes.perday.21.35, na.rm=T),
    nodes.perday.n = length(nodes.perday.21.35[!is.na(nodes.perday.21.35)]),
    leaves.perday.mn = mean(leaves.perday.21.35, na.rm=T),
    leaves.perday.sd = sd(leaves.perday.21.35, na.rm=T),
    leaves.perday.n = length(leaves.perday.21.35[!is.na(leaves.perday.21.35)]),
    height.21.mn = mean(height.21.new, na.rm=T),
    height.21.sd = sd(height.21.new, na.rm=T),
    height.21.n = length(height.21.new[!is.na(height.21.new)]),
    nodes.21.mn = mean(nodes.21.new, na.rm=T),
    nodes.21.sd = sd(nodes.21.new, na.rm=T),
    nodes.21.n = length(nodes.21.new[!is.na(nodes.21.new)]),
    leaves.21.mn = mean(leaves.21.new, na.rm=T),
    leaves.21.sd = sd(leaves.21.new, na.rm=T),
    leaves.21.n = length(leaves.21.new[!is.na(leaves.21.new)]),
    height.35.mn = mean(height.35.new, na.rm=T),
    height.35.sd = sd(height.35.new, na.rm=T),
    height.35.n = length(height.35.new[!is.na(height.35.new)]),
    nodes.35.mn = mean(nodes.35.new, na.rm=T),
    nodes.35.sd = sd(nodes.35.new, na.rm=T),
    nodes.35.n = length(nodes.35.new[!is.na(nodes.35.new)]),
    leaves.35.mn = mean(leaves.35.new, na.rm=T),
    leaves.35.sd = sd(leaves.35.new, na.rm=T),
    leaves.35.n = length(leaves.35.new[!is.na(leaves.35.new)]),
    height.rgr.mn = mean(height.perday.21.35.ln, na.rm=T),
    height.rgr.sd = sd(height.perday.21.35.ln, na.rm=T),
    height.rgr.n = length(height.perday.21.35.ln[!is.na(height.perday.21.35.ln)]),
    leaves.rgr.mn = mean(leaves.perday.21.35.ln, na.rm=T),
    leaves.rgr.sd = sd(leaves.perday.21.35.ln, na.rm=T),
    leaves.rgr.n = length(leaves.perday.21.35.ln[!is.na(leaves.perday.21.35.ln)]),
  ) %>% ungroup ()
View(acc.growth.sum)
write.xlsx(acc.growth.sum, file = "C:/Users/SterlingH/Desktop/R_exports/growth_averaged.xlsx")

# For replicate-level growth data:
rep.growth.sum <- comb.plant.clean.new %>%
  group_by(genus,lifespan,species.short,accession,replicate) %>%
  summarise(
    height.perday.mn = mean(height.perday.21.35, na.rm=T),
    height.perday.sd = sd(height.perday.21.35, na.rm=T),
    height.perday.n = length(height.perday.21.35[!is.na(height.perday.21.35)]),
    nodes.perday.mn = mean(nodes.perday.21.35, na.rm=T),
    nodes.perday.sd = sd(nodes.perday.21.35, na.rm=T),
    nodes.perday.n = length(nodes.perday.21.35[!is.na(nodes.perday.21.35)]),
    leaves.perday.mn = mean(leaves.perday.21.35, na.rm=T),
    leaves.perday.sd = sd(leaves.perday.21.35, na.rm=T),
    leaves.perday.n = length(leaves.perday.21.35[!is.na(leaves.perday.21.35)]),
    height.21.mn = mean(height.21.new, na.rm=T),
    height.21.sd = sd(height.21.new, na.rm=T),
    height.21.n = length(height.21.new[!is.na(height.21.new)]),
    nodes.21.mn = mean(nodes.21.new, na.rm=T),
    nodes.21.sd = sd(nodes.21.new, na.rm=T),
    nodes.21.n = length(nodes.21.new[!is.na(nodes.21.new)]),
    leaves.21.mn = mean(leaves.21.new, na.rm=T),
    leaves.21.sd = sd(leaves.21.new, na.rm=T),
    leaves.21.n = length(leaves.21.new[!is.na(leaves.21.new)]),
    height.35.mn = mean(height.35.new, na.rm=T),
    height.35.sd = sd(height.35.new, na.rm=T),
    height.35.n = length(height.35.new[!is.na(height.35.new)]),
    nodes.35.mn = mean(nodes.35.new, na.rm=T),
    nodes.35.sd = sd(nodes.35.new, na.rm=T),
    nodes.35.n = length(nodes.35.new[!is.na(nodes.35.new)]),
    leaves.35.mn = mean(leaves.35.new, na.rm=T),
    leaves.35.sd = sd(leaves.35.new, na.rm=T),
    leaves.35.n = length(leaves.35.new[!is.na(leaves.35.new)]),
    tray = mean(tray) # used to get the tray number ID for that replicate (since only one tray # per replicate)
  ) %>% ungroup ()
View(rep.growth.sum)
write.xlsx(rep.growth.sum, file = "C:/Users/SterlingH/Desktop/R_exports/growth_rep.xlsx")
# remove any with no data for any of the vegetative traits
rep.nozero <- rep.growth.sum %>%
  filter(!is.na(height.21.mn) | !is.na(height.35.mn) | !is.na(leaves.21.mn) | !is.na(leaves.35.mn)) # removes only the reps with all of this data missing
View(rep.nozero)

# Summarize sampling for vegetative growth
growth.sampling <- rep.nozero %>%
  group_by(genus,lifespan,accession) %>%
  summarise(
    rep.n = n_distinct(replicate, na.rm=T),
    tray.n = n_distinct(tray, na.rm=T)
  ) %>% ungroup ()
View(growth.sampling)
write.xlsx(growth.sampling, file = "C:/Users/SterlingH/Desktop/R_exports/growth_sampling.xlsx")

# For accession-level leaf data:
acc.leaf.sum <- comb.leaf.clean %>%
  group_by(genus,lifespan,species.short,accession) %>%
  summarise(
    leaf.fresh.mn = mean(leaf.fresh, na.rm=T),
    leaf.fresh.sd = sd(leaf.fresh, na.rm=T),
    leaf.dry.mn = mean(leaf.dry, na.rm=T),
    leaf.dry.sd = sd(leaf.dry, na.rm=T),
    leaf.area.mn = mean(leaf.area, na.rm=T),
    leaf.area.sd = sd(leaf.area, na.rm=T),
    ldrm.mn = mean(ldrm, na.rm=T),
    ldrm.sd = sd(ldrm, na.rm=T),
    sla.mn = mean(sla, na.rm=T),
    sla.sd = sd(sla, na.rm=T),
  ) %>% ungroup ()
View(acc.leaf.sum)

# MERGES ####
# Full accession merge (for correlations)

# MERGE OF ACCESSION dataframes
# ADD each new piece to the dataset at the END of the list(), so it doesn't mess up the column numbers
# for the correlation matrices!

# FIRST Merge germ rate rep data (summarized by accession) with existing acc.germ.clean dataframe
merge.germ.new <- Reduce(function(x,y) merge(x,y, all=T, 
                        by=c("genus","lifespan","species.short","accession")), 
                         list(acc.germ.clean, acc.germrep.sum))

# USE THIS FOR FINAL MERGE!
merge.acc.all <- Reduce(function(x,y) merge(x,y, all=T, # all = T ensures that dataframes with no matching accessions are still added to the new merged dataframe
                        by=c("genus","lifespan","species.short","accession")), # have to do whole hierarchy or it creates new columns
                        list(comb.acc.clean, acc.seed.sum, merge.germ.new, acc.growth.sum))# put all data frames to merge here
View(merge.acc.all)
write.xlsx(merge.acc.all, file = "C:/Users/SterlingH/Desktop/R_exports/merg_acc_all.xlsx")

#######################################################

# FULL DATASET NETWORK & CORRELATIONS ####
# Have to start by making a correlation matrix
d <- merge.acc.all %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
d <- scale(d, scale=TRUE, center=TRUE)
x <- cor(d, use="complete.obs") # allows NA's in matrix, & removes them
res1 <- corr.test(d, adjust="holm", use="complete.obs") # Pearson is the default
# correlation coefficients will be the same as x if "complete.obs" is used for corr.test (same as "complete")
# Keep as use="complete.obs" or the p value will be different! Default is pairwise
colnames(res1$p) <- colnames(x)
rownames(res1$p) <- rownames(x)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(x) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(x) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(x, order="original", method="color", col=col,
         p.mat = res1$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Full dataset", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(x) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")

# NETWORK!!!
screen_mat <- res1$p
screen_mat[screen_mat > 0.05] <- 0
x[screen_mat == 0] <- 0
g_raw <- graph_from_adjacency_matrix(x, mode='upper', weighted=T, diag=F) # mode=upper: only the upper right triangle is used for the edge weights
g <- graph_from_adjacency_matrix(abs(x)^2, mode='upper', weighted=T, diag=F) # use r2 for the correlation values, since iGraph can't handle negative corr values.
E(g)$color[E(g_raw)$weight > 0] <- 'blue3'
E(g)$color[E(g_raw)$weight < 0] <- 'brown3'
plot(g, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white') # auto, line thickness increased 5x
plot(g, layout=layout.circle, edge.width=abs(E(g)$weight)*5, vertex.color='white') # circle
plot(g, layout=layout.sphere, edge.width=abs(E(g)$weight)*5, vertex.color='white') # sphere
plot(g, layout=layout.grid, edge.width=abs(E(g)$weight)*5, vertex.color='white') # grid

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg <- degree(g, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col <- colrange(max(deg)+1)
col <- col[deg+1]

# get coordinates of tkplot to replot with the norm plot function
tk <- tkplot(gv, layout=layout_with_fr, edge.width=abs(E(g)$weight)*7, vertex.color=col.vicia, 
                  vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.color="black", main="Full dataset")

# first move the nodes where you want them and then get the coordinates
# Keep the tk plot open!
layout.legume <- tkplot.getcoords(tk)

# matrix coords for future reference
legume.mat <- cbind(c(130,187,41,221,71,144,68,312,397,272,287,304,272,416,399,388,418), 
      c(199,3,127,124,2,370,314,384,356,282,183,90,0,281,184,88,0))
legume.mat==layout.legume # check that it's the same

# Set up grid of plots
par(mar=c(0,0,1,0), mfrow=c(2,3))

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(g, layout=legume.mat, edge.width=abs(E(g)$weight)*7, vertex.color=col, 
     vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("A  Full dataset", adj=0.1, line=-0.5)

# trying changing color based on betweenness
deg <- betweenness(g, directed=T, weights=NA)
colrange <- colorRampPalette(c("yellow", "red"))
col <- colrange(max(deg)+1)
col <- col[deg+1]

plot(g, layout=layout_with_fr, edge.width=abs(E(g)$weight)*5, vertex.color=col, 
     vertex.size=20, vertex.label.cex=0.45, vertex.label.font=2, vertex.label.color="black", main="Full dataset")

tkplot(g, layout=layout_with_fr, edge.width=abs(E(g)$weight)*5, vertex.color=col, 
       vertex.size=18, vertex.label.cex=0.45, vertex.label.font=2, vertex.label.color="black", main="Full dataset")

# ANNUAL NETWORK & CORRELATIONS ####
ann.tot <- merge.acc.all %>% filter(lifespan=="annual")
da <- ann.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
da <- scale(da, scale=TRUE, center=TRUE)
xa <- cor(da, use="complete.obs") # allows NA's in matrix, & removes them
resa <- corr.test(da, adjust="holm", use="complete.obs") # Pearson is the default
# correlation coefficients will be the same as x if "complete.obs" is used for corr.test.
colnames(resa$p) <- colnames(xa)
rownames(resa$p) <- rownames(xa)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xa) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xa) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xa, order="original", method="color", col=col,
         p.mat = resa$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Annual", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xa) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_mata <- resa$p
screen_mata[screen_mata > 0.05] <- 0
xa[screen_mata == 0] <- 0
g_rawa <- graph_from_adjacency_matrix(xa, mode='upper', weighted=T, diag=F)
ga <- graph_from_adjacency_matrix(abs(xa)^2, mode='upper', weighted=T, diag=F)
E(ga)$color[E(g_rawa)$weight > 0] <- 'blue3'
E(ga)$color[E(g_rawa)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.ann <- degree(ga, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.ann <- colrange(max(deg.ann)+1)
col.ann <- col.ann[deg.ann+1]

# Keep the tk plot open!
# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(ga, layout=legume.mat, edge.width=abs(E(ga)$weight)*7, vertex.color=col.ann, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("B  Annual", adj=0.1, line=-0.5)

# PERENNIAL NETWORK & CORRELATIONS ####
per.tot <- merge.acc.all %>% filter(lifespan=="perennial")
dper <- per.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dper <- scale(dper, scale=TRUE, center=TRUE)
xper <- cor(dper, use="complete.obs") # allows NA's in matrix, & removes them
resper <- corr.test(dper, adjust="holm", use="complete.obs") # Pearson is the default
# correlation coefficients will be the same as x if "complete.obs" is used for corr.test.
# Keep as use="complete.obs" or the p value will be different! Default is pairwise
colnames(resper$p) <- colnames(xper)
rownames(resper$p) <- rownames(xper)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xper) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xper) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xper, order="original", method="color", col=col,
         p.mat = resper$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Perennial", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xper) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matper <- resper$p
screen_matper[screen_matper > 0.05] <- 0
xper[screen_matper == 0] <- 0
g_rawper <- graph_from_adjacency_matrix(xper, mode='upper', weighted=T, diag=F)
gper <- graph_from_adjacency_matrix(abs(xper)^2, mode='upper', weighted=T, diag=F)
E(gper)$color[E(g_rawper)$weight > 0] <- 'blue3'
E(gper)$color[E(g_rawper)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.per <- degree(gper, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.per <- colrange(max(deg.per)+1)
col.per <- col.per[deg.per+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gper, layout=legume.mat, edge.width=abs(E(gper)$weight)*7, vertex.color=col.per, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("C  Perennial", adj=0.1, line=-0.5)

# LATHYRUS NETWORK & CORRELATIONS ####
lath.tot <- merge.acc.all %>% filter(genus=="Lathyrus")
dl <- lath.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dl <- scale(dl, scale=TRUE, center=TRUE)
xl <- cor(dl, use="complete.obs") # allows NA's in matrix, & removes them
resl <- corr.test(dl, adjust="holm", use="complete.obs") # Pearson is the default
# correlation coefficients will be the same as x if "complete.obs" is used for corr.test.
# Keep as use="complete.obs" or the p value will be different! Default is pairwise
colnames(resl$p) <- colnames(xl)
rownames(resl$p) <- rownames(xl)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xl) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xl) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xl, order="original", method="color", col=col,
         p.mat = resl$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Lathyrus", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xl) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matl <- resl$p
screen_matl[screen_matl > 0.05] <- 0
xl[screen_matl == 0] <- 0
g_rawl <- graph_from_adjacency_matrix(xl, mode='upper', weighted=T, diag=F)
gl <- graph_from_adjacency_matrix(abs(xl)^2, mode='upper', weighted=T, diag=F)
E(gl)$color[E(g_rawl)$weight > 0] <- 'blue3'
E(gl)$color[E(g_rawl)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.lath <- degree(gl, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.lath <- colrange(max(deg.lath)+1)
col.lath <- col.lath[deg.lath+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gl, layout=legume.mat, edge.width=abs(E(gl)$weight)*7, vertex.color=col.lath, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("D  Lathyrus", adj=0.1, line=-0.5)


# PHASEOLUS NETWORK & CORRELATIONS ####
phas.tot <- merge.acc.all %>% filter(genus=="Phaseolus")
dp <- phas.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dp <- scale(dp, scale=TRUE, center=TRUE)
xp <- cor(dp, use="complete.obs") # allows NA's in matrix, & removes them
resp <- corr.test(dp, adjust="holm", use="complete.obs") # Pearson is the default
# correlation coefficients will be the same as x if "complete.obs" is used for corr.test.
# Keep as use="complete.obs" or the p value will be different! Default is pairwise
colnames(resp$p) <- colnames(xp)
rownames(resp$p) <- rownames(xp)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xp) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xp) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xp, order="original", method="color", col=col,
         p.mat = resp$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Phaseolus", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xp) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matp <- resp$p
screen_matp[screen_matp > 0.05] <- 0
xp[screen_matp == 0] <- 0
g_rawp <- graph_from_adjacency_matrix(xp, mode='upper', weighted=T, diag=F)
gp <- graph_from_adjacency_matrix(abs(xp)^2, mode='upper', weighted=T, diag=F)
E(gp)$color[E(g_rawp)$weight > 0] <- 'blue3'
E(gp)$color[E(g_rawp)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.phas <- degree(gp, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.phas <- colrange(max(deg.phas)+1)
col.phas <- col.phas[deg.phas+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gp, layout=legume.mat, edge.width=abs(E(gp)$weight)*7, vertex.color=col.phas, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("E  Phaseolus", adj=0.1, line=-0.5)


# VICIA NETWORK & CORRELATIONS ####
vicia.tot <- merge.acc.all %>% filter(genus=="Vicia")
dv <- vicia.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dv <- scale(dv, scale=TRUE, center=TRUE)
xv <- cor(dv, use="complete.obs") # allows NA's in matrix, & removes them
resv <- corr.test(dv, adjust="holm", use="complete.obs") # Pearson is the default
# correlation coefficients will be the same as x if "complete.obs" is used for corr.test.
# Keep as use="complete.obs" or the p value will be different! Default is pairwise
colnames(resv$p) <- colnames(xv)
rownames(resv$p) <- rownames(xv)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xv) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xv) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xv, order="original", method="color", col=col,
         p.mat = resv$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Vicia", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xv) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matv <- resv$p
screen_matv[screen_matv > 0.05] <- 0
xv[screen_matv == 0] <- 0
g_rawv <- graph_from_adjacency_matrix(xv, mode='upper', weighted=T, diag=F)
gv <- graph_from_adjacency_matrix(abs(xv)^2, mode='upper', weighted=T, diag=F)
E(gv)$color[E(g_rawv)$weight > 0] <- 'blue3'
E(gv)$color[E(g_rawv)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.vicia <- degree(gv, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.vicia <- colrange(max(deg.vicia)+1)
col.vicia <- col.vicia[deg.vicia+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gv, layout=legume.mat, edge.width=abs(E(gv)$weight)*7, vertex.color=col.vicia, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("F  Vicia", adj=0.1, line=-0.5)

# testing specific relationships
plot <- ggplot(data=lath.tot, aes(x=ppd50.mn, y=nodes.21.mn, color=lifespan))+
  geom_point()+labs(title="annual")
plot

plot <- ggplot(data=merge.acc.all, aes(x=ppd50.mn, y=nodes.35.mn, color=genus))+
  geom_point()+labs(title="all")
plot

merge.acc.all
# Two methods of making a grid with the network plots:

# make a grid
map_base_to_grid <- function(fun) {
  gridGraphics::grid.echo(fun)
  grid::grid.grab()
}

grid1 <- map_base_to_grid(function() plot(g, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid2 <- map_base_to_grid(function() plot(gl, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid3 <- map_base_to_grid(function() plot(gp, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid4 <- map_base_to_grid(function() plot(gv, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid5 <- map_base_to_grid(function() plot(ga, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid6 <- map_base_to_grid(function() plot(gper, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))

ggpubr::ggarrange(grid1, grid2, grid3, grid4, grid5, grid6, nrow=3, ncol=2, widths=c(2,2), heights=c(2,2,2))

plotit <- function(input){plot(input, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
     vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0)}
gridfun.all <- function()plotit(input=g)
gridfun.lath <- function()plotit(input=gl)
gridfun.phas <- function()plotit(input=gp)
gridfun.vic <- function()plotit(input=gv)
gridfun.ann <- function()plotit(input=ga)
gridfun.per <- function()plotit(input=gper)

ggpubr::ggarrange(gridfun.all, gridfun.lath, gridfun.phas, gridfun.vic, gridfun.ann, gridfun.per, nrow=3, ncol=2, widths=c(2,2), heights=c(2,2,2))

# Network For PCs
# merge PCA outputs
merge.pca <- Reduce(function(x,y) merge(x,y, all=T, by=c("accession")), list(seed.acc.pca, veg.acc.pca))
d.pca <- merge.pca %>% dplyr::select(sPC1:vPC4)
d.pca <- scale(d.pca, scale=TRUE, center=TRUE)
x.pca <- cor(d.pca, use="complete.obs") # allows NA's in matrix, & removes them
res.pca <- cor.mtest(d.pca, conf.level = .95)
colnames(res.pca$p) <- colnames(x.pca)
rownames(res.pca$p) <- rownames(x.pca)
screen_mat.pca <- res.pca$p
screen_mat.pca[screen_mat.pca > 0.05] <- 0
x.pca[screen_mat.pca == 0] <- 0
g_raw.pca <- graph_from_adjacency_matrix(x.pca, mode='upper', weighted=T, diag=F)
g.pca <- graph_from_adjacency_matrix(abs(x.pca)^2, mode='upper', weighted=T, diag=F)
E(g.pca)$color[E(g_raw.pca)$weight > 0] <- 'skyblue1'
E(g.pca)$color[E(g_raw.pca)$weight < 0] <- 'salmon1'
plot(g.pca, layout=layout.auto, edge.width=abs(E(g.pca)$weight)*30, 
     vertex.color='white', vertex.size=23, vertex.label.cex=0.7)
plot(g.pca, layout=layout.circle, edge.width=abs(E(g.pca)$weight)*5, vertex.color='white') # circle
plot(g.pca, layout=layout.sphere, edge.width=abs(E(g.pca)$weight)*5, vertex.color='white') # sphere
plot(g.pca, layout=layout.grid, edge.width=abs(E(g.pca)$weight)*5, vertex.color='white') # grid

#######################################################

# PCA - all traits ####
few.traits <- c("W6 2747", "PI 420171", "PI 358829", "PI 250758", "PI 535200", "PI 494749") # drop accessions that lack many trait data
# W6 2747 - no germ rate & proportion is weird (dormant), no veg; 
# PI 420171 - no veg
# PI 358829 - missing most DAP-21 and all growth rate data
# PI 250758 - missing all DAP-21 and all growth rate data
# Others are kept since only missing one trait

# Short names if needed again: colnames(a) <- c("Sm", "Sl", "Sw", "--Sp", "Sa", "Sc", "Sr", "Gt", "Gp", "H21", "H35", "Hgr", "Hrgr", "L21", "L35", "Lgr")

merge.acc.all.d <- merge.acc.all[!merge.acc.all$accession %in% few.traits,]
a <- merge.acc.all.d %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
colnames(a) <- c("S mass", "S length", "S width", "          --S per.", "S area", "S circularity", "S roundness", "G T50", "G proportion", "H DAP-21", "H DAP-35\n", "H AGR", "H RGR", "L DAP-21", "L DAP-35", "L AGR", "L RGR")
estim_ncpPCA(a, scale=T) # estimates # of dimensions needed in the next imputation step (ncp)
a1 <- imputePCA(as.data.frame(a), ncp=3, scale=T) # new df with imputed values in place of NAs
a.pca <- prcomp(a1$completeObs, scale=T)
get_eigenvalue(a.pca) # Eigenvalues for PCs
fviz_eig(a.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(a.pca, col.var = "black") # correlation circle
fviz_cos2(a.pca, choice = "var", axes = 1:2) # Cos2, quality of representation
fviz_contrib(a.pca, choice = "var", axes = 1, top = 10) # Contributions of variables to PC1
fviz_contrib(a.pca, choice = "var", axes = 2, top = 10) # Contributions of variables to PC2
a.pca.loadings <- a.pca$rotation # PC loadings from each variable
sqrt(1/nrow(a.pca$rotation)) # to determine the threshold above which is considered a "strong" loading
#a.pca$rotation[,1] <- -a.pca$rotation[,1] # change the direction of PC1 to flip seed size to be on the right (positive)
#a.pca$x[,1] <- -a.pca$x[,1]

write.xlsx(merge.acc.all.d, file = "C:/Users/SterlingH/Desktop/R_exports/full_names.xlsx") # For names
write.xlsx(a, file = "C:/Users/SterlingH/Desktop/R_exports/full_compare.xlsx")
write.xlsx(a1, file = "C:/Users/SterlingH/Desktop/R_exports/full_impute.xlsx")

# Vairables plot from ggbiplot
# ggbiplot2 is the code modified to allow changes in variable and arrow attributes
# there is a separate R script for the new ggbiplot2 function, since its' very lengthy
# use ggbiplot() for original code

a.pca.varplot <- ggbiplot2(a.pca, circle=T, alpha=0, varname.adjust=1.1, 
                           varname.color="darkblue", arrow.color="black", circle.color="gray20")+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  scale_x_continuous(limits=c(-2.25, 2.25), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2.35, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+ 
  labs(x="PC1 (42.2%)", y="PC2 (24.7%)")+ #title="Variables"
  theme_classic2()
a.pca.varplot

# Make figures with ggplot
a.pca.df <- as.data.frame(a.pca$x)
# add categories to df
a.pca.df$lifespan <- merge.acc.all.d$lifespan
a.pca.df$genus <- merge.acc.all.d$genus
a.pca.df$species.short <- merge.acc.all.d$species.short

# Colored by lifespan
a.pca.lifespan <- ggplot(data=a.pca.df, aes(x=PC1, y=PC2, shape=genus, fill=lifespan, color=lifespan))+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_point(size=2.5, color="black")+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Life span",
                    labels = c("annual", "perennial"),
                    values = c("brown3", "cornflowerblue"))+
  scale_colour_manual(name = "Life span",
                      labels = c("annual", "perennial"),
                      values = c("brown3", "cornflowerblue"))+
  #stat_ellipse(data=a.pca.df, aes(x=PC1, y=PC2, group=lifespan, color=lifespan), level=0.95)+
  labs(x="PC1 (42.2%)", y="PC2 (24.7%)")+ #title="Life span"
  scale_x_continuous(limits=c(-6, 9), breaks=seq(-5, 10, 2.5), expand=c(0, 0))+
  scale_y_continuous(limits=c(-6, 5), breaks=seq(-5, 5, 2.5), expand=c(0, 0))+
  theme_classic2()+
  guides(shape = F)
a.pca.lifespan

# old ellipse code
# stat_conf_ellipse(data=a.pca.df, aes(x=PC1, y=PC2, group=lifespan, color=lifespan))+

# Colored by genus
a.pca.genus <- ggplot(data=a.pca.df, aes(x=PC1, y=PC2, shape=genus, fill=genus, color=genus))+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_point(size=2.5, color="black")+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("orange", "turquoise3", "seashell3"))+
  #stat_ellipse(data=a.pca.df, aes(x=PC1, y=PC2, group=genus, color=genus), level=0.95)+
  labs(x="PC1 (42.2%)", y="PC2 (24.7%)")+ #title="Genus"
  scale_x_continuous(limits=c(-6, 9), breaks=seq(-5, 10, 2.5), expand=c(0, 0))+
  scale_y_continuous(limits=c(-6, 5), breaks=seq(-5, 5, 2.5), expand=c(0, 0))+
  theme_classic2()
a.pca.genus

# old ellipse code
# stat_conf_ellipse(data=a.pca.df, aes(x=PC1, y=PC2, group=genus, color=genus))+

# GENUS SEED PCAs ####

### SEED PCAs
# example for how to change axis direction: ls.pca$rotation[,1] <- -ls.pca$rotation[,1] 

#### LATHYRUS SEED PCA

lath.seed <- comb.seed.clean %>% filter(genus=="Lathyrus")
ls <- lath.seed %>% select(seed.length, seed.width, seed.perim, seed.area, seed.circ, seed.roundness)
colnames(ls) <- c("length", "width", "perimeter", "area\n", "circularity", "roundness")
ls.pca <- prcomp(ls, scale=T)
get_eigenvalue(ls.pca) # Eigenvalues for PCs
fviz_eig(ls.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(ls.pca, col.var = "black", title="Lathyrus seed PCA variables")
ls.pca.loading <- ls.pca$rotation # PC loadings from each variable
sqrt(1/nrow(ls.pca$rotation))
ls.pca$rotation[,1] <- -ls.pca$rotation[,1] # change the direction of PC1 to flip seed size to be on the right (positive)
ls.pca$x[,1] <- -ls.pca$x[,1] # have to also change the PC scores to flip accordingly
ls.pca$rotation[,2] <- -ls.pca$rotation[,2]
ls.pca$x[,2] <- -ls.pca$x[,2]

# ggbiplot
ls.pca.varplot <- ggbiplot2(ls.pca, circle=T, alpha=0, varname.adjust=1.1, 
                            varname.color="darkblue", arrow.color="black", circle.color="gray20")+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  labs(x="PC1 (69.0%)", y="PC2 (23.2%)", title="Seed variables")+
  theme_classic2()
ls.pca.varplot

# LATHYRUS SPECIES SEED PCA
brewer.pal(n=6, name="Reds")
brewer.pal(n=6, name="Blues")

ls.pca.df <- as.data.frame(ls.pca$x)

# std <- function(x) sd(x)/sqrt(length(x)) # standard error function

# summarize by species - STILL use this!
ls.pca.spsum <- ls.pca.df %>%
  group_by(lath.seed$lifespan, lath.seed$specific.epithet) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(ls.pca.spsum)
names(ls.pca.spsum)[1] <- "lifespan" 
names(ls.pca.spsum)[2] <- "species" 
ls.pca.spsum$species <- factor(ls.pca.spsum$species, levels=c("annuus", "aphaca", "cicera", "hirsutus", "odoratus", "japonicus", "latifolius", "pratensis", "sylvestris", "tuberosus"), ordered=T)
  
ls.pca.plot <- ggplot()+
  geom_hline(yintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  scale_shape_manual(name = "Species",
                     values = c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25))+
  scale_colour_manual(name = "Life span",
                    labels = c("annual", "perennial"),
                    values = c("red", "blue"))+
  geom_errorbar(data=ls.pca.spsum, aes(x=PC1.mn, ymin=PC2.mn-PC2.sd, ymax=PC2.mn+PC2.sd, color=lifespan), width=0.1, size=0.35)+
  geom_errorbarh(data=ls.pca.spsum, aes(y=PC2.mn, xmin=PC1.mn-PC1.sd, xmax=PC1.mn+PC1.sd, color=lifespan), height=0.1, size=0.35)+
  geom_point(data=ls.pca.spsum, aes(x=PC1.mn, y=PC2.mn, shape=species, fill=species), color="black", size=3)+
  scale_fill_manual(name = "Species",
                    values = c("#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))+
  labs(x="PC1 (69.0%)", y="PC2 (23.2%)", title="Seed PCA")+
  scale_x_continuous(limits=c(-4, 4), breaks=seq(-4, 4, 2), expand=c(0, 0))+
  scale_y_continuous(limits=c(-4, 2), breaks=seq(-4, 2, 2), expand=c(0, 0))+
  theme_classic2()
ls.pca.plot

# summarize by accession (not using this anymore)
ls.pca.acsum <- ls.pca.df %>%
  group_by(lath.seed$lifespan, lath.seed$specific.epithet, lath.seed$accession) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(ls.pca.acsum)
names(ls.pca.acsum)[1] <- "lifespan" 
names(ls.pca.acsum)[2] <- "species" 
names(ls.pca.acsum)[3] <- "accession" 
ls.pca.acsum$species <- factor(ls.pca.acsum$species, levels=c("annuus", "aphaca", "cicera", "hirsutus", "odoratus", "japonicus", "latifolius", "pratensis", "sylvestris", "tuberosus"), ordered=T)

# If using accession again, put this line at the top under ggplot()+ :  geom_point(data=ls.pca.acsum, aes(x=PC1.mn, y=PC2.mn, shape=species, color=lifespan), size=2)+

# for now, not using ellipse
# stat_conf_ellipse(data=ls.pca.df, aes(x=PC1, y=PC2, group=lath.seed$species.short), linetype="dashed")+

#### PHASEOLUS SEED PCA

phas.seed <- comb.seed.clean %>% filter(genus=="Phaseolus")
ps <- phas.seed %>% select(seed.length, seed.width, seed.perim, seed.area, seed.circ, seed.roundness) 
colnames(ps) <- c("length", "width", "perimeter", "area\n", "circularity", "roundness")
ps.pca <- prcomp(ps, scale=T)
get_eigenvalue(ps.pca) # Eigenvalues for PCs
fviz_eig(ps.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(ps.pca, col.var = "black", title="Phaseolus seed PCA variables")
ps.pca.loading <- ps.pca$rotation
# Needs to be flipped on both PC1 & PC2
ps.pca$rotation[,1] <- -ps.pca$rotation[,1] # change the direction of PC1 to flip seed size to be on the right (positive)
ps.pca$x[,1] <- -ps.pca$x[,1] # have to also change the PC scores to flip accordingly
ps.pca$rotation[,2] <- -ps.pca$rotation[,2]
ps.pca$x[,2] <- -ps.pca$x[,2]

# ggbiplot
ps.pca.varplot <- ggbiplot2(ps.pca, circle=T, alpha=0, varname.adjust=1.1, 
                            varname.color="darkblue", arrow.color="black", circle.color="gray20")+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  labs(x="PC1 (67.1%)", y="PC2 (23.3%%)", title="Seed variables")+
  theme_classic2()
ps.pca.varplot

# PHASEOLUS SPECIES SEED PCA

ps.pca.df <- as.data.frame(ps.pca$x)

# summarize by species - USE
ps.pca.spsum <- ps.pca.df %>%
  group_by(phas.seed$lifespan, phas.seed$specific.epithet) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(ps.pca.spsum)
names(ps.pca.spsum)[1] <- "lifespan" 
names(ps.pca.spsum)[2] <- "species"
ps.pca.spsum$species <- factor(ps.pca.spsum$species, levels=c("acutifolius", "filiformis", "lunatus", "vulgaris", "angustissimus", "leptostachyus", "maculatus", "parvulus", "polystachios"), ordered=T)

ps.pca.plot <- ggplot()+
  geom_hline(yintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  scale_shape_manual(name = "Species",
                     values = c(21, 22, 23, 24, 21, 22, 23, 24,25))+
  scale_colour_manual(name = "Life span",
                      labels = c("annual", "perennial"),
                      values = c("red", "blue"))+
  geom_errorbar(data=ps.pca.spsum, aes(x=PC1.mn, ymin=PC2.mn-PC2.sd, ymax=PC2.mn+PC2.sd, color=lifespan), width=0.1, size=0.35)+
  geom_errorbarh(data=ps.pca.spsum, aes(y=PC2.mn, xmin=PC1.mn-PC1.sd, xmax=PC1.mn+PC1.sd, color=lifespan), height=0.1, size=0.35)+
  geom_point(data=ps.pca.spsum, aes(x=PC1.mn, y=PC2.mn, shape=species, fill=species), color="black", size=3)+
  scale_fill_manual(name = "Species",
                    values = c("#FC9272", "#FB6A4A", "#DE2D26", "#A50F15", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))+
  labs(x="PC1 (67.1%)", y="PC2 (23.3%)", title="Seed PCA")+
  scale_x_continuous(limits=c(-3, 5), breaks=seq(-2, 4, 2), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2, 2.5), breaks=seq(-2, 2, 2), expand=c(0, 0))+
  theme_classic2()
ps.pca.plot

# summarize by accession - NOT USING
ps.pca.acsum <- ps.pca.df %>%
  group_by(phas.seed$lifespan, phas.seed$specific.epithet, phas.seed$accession) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(ps.pca.acsum)
names(ps.pca.acsum)[1] <- "lifespan" 
names(ps.pca.acsum)[2] <- "species" 
names(ps.pca.acsum)[3] <- "accession"
ps.pca.acsum$species <- factor(ps.pca.acsum$species, levels=c("acutifolius", "filiformis", "lunatus", "vulgaris", "angustissimus", "leptostachyus", "maculatus", "parvulus", "polystachios"), ordered=T)

# stat_conf_ellipse(data=ps.pca.df, aes(x=PC1, y=PC2, group=phas.seed$species.short), linetype="dashed")+
# If using accession again, put this line at the top under ggplot()+ :  geom_point(data=ps.pca.acsum, aes(x=PC1.mn, y=PC2.mn, shape=species, color=lifespan), size=2)+

#### VICIA SEED PCA

vicia.seed <- comb.seed.clean %>% filter(genus=="Vicia")
vs <- vicia.seed %>% select(seed.length, seed.width, seed.perim, seed.area, seed.circ, seed.roundness) # changed order to make PC2 orientation comparable to Lathyrus PC2
colnames(vs) <- c("length", "width", "perimeter", "area\n", "circularity", "roundness")
vs.pca <- prcomp(vs, scale=T)
get_eigenvalue(vs.pca) # Eigenvalues for PCs
fviz_eig(vs.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(vs.pca, col.var = "black", title="Vicia seed PCA variables")
vs.pca.loading <- vs.pca$rotation
# Needs to be flipped on both PC1 & PC2
vs.pca$rotation[,2] <- -vs.pca$rotation[,2]
vs.pca$x[,2] <- -vs.pca$x[,2]

# ggbiplot
vs.pca.varplot <- ggbiplot2(vs.pca, circle=T, alpha=0, varname.adjust=1.1, 
                            varname.color="darkblue", arrow.color="black", circle.color="gray20")+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  labs(x="PC1 (69.9%)", y="PC2 (22.4%%)", title="Seed variables")+
  theme_classic2()
vs.pca.varplot

# VICIA SPECIES SEED PCA

vs.pca.df <- as.data.frame(vs.pca$x)

# summarize by species - USE
vs.pca.spsum <- vs.pca.df %>%
  group_by(vicia.seed$lifespan, vicia.seed$specific.epithet) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(vs.pca.spsum)
names(vs.pca.spsum)[1] <- "lifespan" 
names(vs.pca.spsum)[2] <- "species" 
vs.pca.spsum$species <- factor(vs.pca.spsum$species, levels=c("benghalensis", "ervilia", "hirsuta", "sativa", "villosa", "americana", "cracca", "dumetorum", "sepium", "tenuifolia"), ordered=T)

vs.pca.plot <- ggplot()+
  geom_hline(yintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  scale_shape_manual(name = "Species",
                     values = c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25))+
  scale_colour_manual(name = "Life span",
                      labels = c("annual", "perennial"),
                      values = c("red", "blue"))+
  geom_errorbar(data=vs.pca.spsum, aes(x=PC1.mn, ymin=PC2.mn-PC2.sd, ymax=PC2.mn+PC2.sd, color=lifespan), width=0.1, size=0.35)+
  geom_errorbarh(data=vs.pca.spsum, aes(y=PC2.mn, xmin=PC1.mn-PC1.sd, xmax=PC1.mn+PC1.sd, color=lifespan), height=0.1, size=0.35)+
  geom_point(data=vs.pca.spsum, aes(x=PC1.mn, y=PC2.mn, shape=species, fill=species), color="black", size=3)+
  scale_fill_manual(name = "Species",
                    values = c("#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))+
  labs(x="PC1 (69.9%)", y="PC2 (22.4%)", title="Seed PCA")+
  scale_x_continuous(limits=c(-4, 4), breaks=seq(-2, 4, 2), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2, 2), breaks=seq(-2, 2, 2), expand=c(0, 0))+
  theme_classic2()
vs.pca.plot

# stat_conf_ellipse(data=vs.pca.df, aes(x=PC1, y=PC2, group=vicia.seed$species.short), linetype="dashed")+
# If using accession again, put this line at the top under ggplot()+ :  geom_point(data=vs.pca.acsum, aes(x=PC1.mn, y=PC2.mn, shape=species, color=lifespan), size=2)+

# summarize by accession - NOT USING
vs.pca.acsum <- vs.pca.df %>%
  group_by(vicia.seed$lifespan, vicia.seed$specific.epithet, vicia.seed$accession) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(vs.pca.acsum)
names(vs.pca.acsum)[1] <- "lifespan"
names(vs.pca.acsum)[2] <- "species" 
names(vs.pca.acsum)[3] <- "accession" 
vs.pca.acsum$species <- factor(vs.pca.acsum$species, levels=c("benghalensis", "ervilia", "hirsuta", "sativa", "villosa", "americana", "cracca", "dumetorum", "sepium", "tenuifolia"), ordered=T)

# GENUS VEGETATIVE PCAs ####

#### LATHYRUS VEGETATIVE PCA
lath.bad <- c("PI 358829", "PI 250758") # accessions recently removed due to outliers & low sampling
lath.veg <- comb.plant.clean.new %>% filter(genus=="Lathyrus")
lath.veg <- lath.veg[!is.na(lath.veg$height.perday.21.35),] # remove all rows where height AGR = NA, since this only co-occurs when one date has no data at all (except W6 4872, which is gone)
lath.veg <- lath.veg[!lath.veg$accession %in% lath.bad,]
lv <- lath.veg %>% select(height.21.new, height.35.new, height.perday.21.35, height.perday.21.35.ln, leaves.21.new, leaves.35.new, leaves.perday.21.35, leaves.perday.21.35.ln)
colnames(lv) <- c("H DAP-21", "H DAP-35", "H AGR", "H RGR", "L DAP-21", "L DAP-35", "L AGR", "L RGR")
estim_ncpPCA(lv, scale=T) # estimates # of dimenstions needed in the next imputation step (ncp)
lv1 <- imputePCA(as.data.frame(lv), ncp=5, scale=T) # imputes missing values (NAs)
lv.pca <- prcomp(lv1$completeObs, scale=T)
get_eigenvalue(lv.pca) # Eigenvalues for PCs
fviz_eig(lv.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(lv.pca, col.var = "black", title="Lathyrus vegetative PCA variables")
lv.pca.loading <- lv.pca$rotation
sqrt(1/nrow(lv.pca$rotation))
lv.pca$rotation[,1] <- -lv.pca$rotation[,1] # flip to have PC space direction consistent between genera
lv.pca$x[,1] <- -lv.pca$x[,1]

# check imputation
write.xlsx(lath.veg, file = "C:/Users/SterlingH/Desktop/R_exports/pre_lath_veg.xlsx") # For names
write.xlsx(lath.veg, file = "C:/Users/SterlingH/Desktop/R_exports/lath_veg_names1.xlsx") # For names
write.xlsx(lv, file = "C:/Users/SterlingH/Desktop/R_exports/lath_veg_compare.xlsx")
write.xlsx(lv1, file = "C:/Users/SterlingH/Desktop/R_exports/lath_veg_impute.xlsx")

# ggbiplot
lv.pca.varplot <- ggbiplot2(lv.pca, circle=T, alpha=0, varname.adjust=1.1, 
                            varname.color="darkblue", arrow.color="black", circle.color="gray20")+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  labs(x="PC1 (35.4%)", y="PC2 (34.3%)", title="Vegetative variables")+
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  theme_classic2()
lv.pca.varplot

# LATHYRUS SPECIES VEGETATIVE PCA

lv.pca.df <- as.data.frame(lv.pca$x)

# std <- function(x) sd(x)/sqrt(length(x)) # standard error function

# summarize by species - USE
lv.pca.spsum <- lv.pca.df %>%
  group_by(lath.veg$lifespan, lath.veg$specific.epithet) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(lv.pca.spsum)
names(lv.pca.spsum)[1] <- "lifespan" 
names(lv.pca.spsum)[2] <- "species" 
lv.pca.spsum$species <- factor(lv.pca.spsum$species, levels=c("aphaca", "cicera", "hirsutus", "odoratus", "japonicus", "latifolius", "pratensis", "sylvestris", "tuberosus"), ordered=T)

lv.pca.plot <- ggplot()+
  geom_hline(yintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  scale_shape_manual(name = "Species",
                     values = c(22, 23, 24, 25, 21, 22, 23, 24, 25))+
  scale_colour_manual(name = "Life span",
                      labels = c("annual", "perennial"),
                      values = c("red", "blue"))+
  geom_errorbar(data=lv.pca.spsum, aes(x=PC1.mn, ymin=PC2.mn-PC2.sd, ymax=PC2.mn+PC2.sd, color=lifespan), width=0.1, size=0.35)+
  geom_errorbarh(data=lv.pca.spsum, aes(y=PC2.mn, xmin=PC1.mn-PC1.sd, xmax=PC1.mn+PC1.sd, color=lifespan), height=0.1, size=0.35)+
  geom_point(data=lv.pca.spsum, aes(x=PC1.mn, y=PC2.mn, shape=species, fill=species), color="black", size=3)+
  scale_fill_manual(name = "Species",
                    values = c("#FC9272", "#FB6A4A", "#DE2D26", "#A50F15", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))+
  labs(x="PC1 (35.4%)", y="PC2 (34.3%)", title="Vegetative PCA")+
  scale_x_continuous(limits=c(-3.5, 4), breaks=seq(-2, 4, 2), expand=c(0, 0))+
  scale_y_continuous(limits=c(-4, 3), breaks=seq(-4, 2, 2), expand=c(0, 0))+
  theme_classic2()
lv.pca.plot

# stat_conf_ellipse(data=lv.pca.df, aes(x=PC1, y=PC2, group=lath.veg$species.short), linetype="dashed")+

# summarize by accession - NOT USING
lv.pca.acsum <- lv.pca.df %>%
  group_by(lath.veg$lifespan, lath.veg$specific.epithet, lath.veg$accession) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(lv.pca.acsum)
names(lv.pca.acsum)[1] <- "lifespan" 
names(lv.pca.acsum)[2] <- "species"
names(lv.pca.acsum)[3] <- "accession"
lv.pca.acsum$species <- factor(lv.pca.acsum$species, levels=c("annuus", "aphaca", "cicera", "hirsutus", "odoratus", "japonicus", "latifolius", "pratensis", "sylvestris", "tuberosus"), ordered=T)

#### PHASEOLUS VEGETATIVE PCA

phas.veg <- comb.plant.clean.new %>% filter(genus=="Phaseolus")
phas.veg <- phas.veg[!is.na(phas.veg$height.perday.21.35),] # remove all rows where height AGR = NA, since this only co-occurs when one date has no data at all or all veg data are NA
pv <- phas.veg %>% select(height.21.new, height.35.new, height.perday.21.35, height.perday.21.35.ln, leaves.21.new, leaves.35.new, leaves.perday.21.35, leaves.perday.21.35.ln)
colnames(pv) <- c("H DAP-21", "H DAP-35", "H AGR", "H RGR", "L DAP-21", "L DAP-35", "L AGR", "L RGR")
# estim_ncpPCA(pv, scale=T) # Phaseolus has no NAs after dropping rows with a lot of NAs, so not necessary to
# pv1 <- imputePCA(as.data.frame(pv), ncp=2, scale=T)
# There are no missing values now after the data reduction above.
pv.pca <- prcomp(pv, scale=T)
get_eigenvalue(pv.pca) # Eigenvalues for PCs
fviz_eig(pv.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(pv.pca, col.var = "black", title="Phaseolus vegetative PCA variables")
pv.pca.loading <- pv.pca$rotation
pv.pca$rotation[,2] <- -pv.pca$rotation[,2] # flip to have PC space direction consistent between genera
pv.pca$x[,2] <- -pv.pca$x[,2]

# check imputation
write.xlsx(phas.veg, file = "C:/Users/SterlingH/Desktop/R_exports/pre_phas_veg_names.xlsx") # For names
write.xlsx(phas.veg, file = "C:/Users/SterlingH/Desktop/R_exports/phas_veg_names.xlsx") # For names
write.xlsx(pv, file = "C:/Users/SterlingH/Desktop/R_exports/phas_veg_compare.xlsx")
write.xlsx(pv1, file = "C:/Users/SterlingH/Desktop/R_exports/phas_veg_impute.xlsx")

# ggbiplot
pv.pca.varplot <- ggbiplot2(pv.pca, circle=T, alpha=0, varname.adjust=1.1, 
                            varname.color="darkblue", arrow.color="black", circle.color="gray20")+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  labs(x="PC1 (48.7%)", y="PC2 (37.9%)", title="Vegetative variables")+
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  theme_classic2()
pv.pca.varplot

# PHASEOLUS SPECIES VEGETATIVE PCA

pv.pca.df <- as.data.frame(pv.pca$x)

# summarize by species - USE
pv.pca.spsum <- pv.pca.df %>%
  group_by(phas.veg$lifespan, phas.veg$specific.epithet) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(pv.pca.spsum)
names(pv.pca.spsum)[1] <- "lifespan" 
names(pv.pca.spsum)[2] <- "species" 
pv.pca.spsum$species <- factor(pv.pca.spsum$species, levels=c("acutifolius", "filiformis", "lunatus", "vulgaris", "angustissimus", "leptostachyus", "maculatus", "parvulus", "polystachios"), ordered=T)

pv.pca.plot <- ggplot()+
  geom_hline(yintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  scale_shape_manual(name = "Species",
                     values = c(21, 22, 23, 24, 21, 22, 23, 24,25))+
  scale_colour_manual(name = "Life span",
                      labels = c("annual", "perennial"),
                      values = c("red", "blue"))+
  geom_errorbar(data=pv.pca.spsum, aes(x=PC1.mn, ymin=PC2.mn-PC2.sd, ymax=PC2.mn+PC2.sd, color=lifespan), width=0.1, size=0.35)+
  geom_errorbarh(data=pv.pca.spsum, aes(y=PC2.mn, xmin=PC1.mn-PC1.sd, xmax=PC1.mn+PC1.sd, color=lifespan), height=0.1, size=0.35)+
  geom_point(data=pv.pca.spsum, aes(x=PC1.mn, y=PC2.mn, shape=species, fill=species), color="black", size=3)+
  scale_fill_manual(name = "Species",
                    values = c("#FC9272", "#FB6A4A", "#DE2D26", "#A50F15", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))+
  labs(x="PC1 (48.7%)", y="PC2 (37.9%)", title="Vegetative PCA")+
  scale_x_continuous(limits=c(-5, 4), breaks=seq(-4, 4, 2), expand=c(0, 0))+
  scale_y_continuous(limits=c(-3, 3.5), breaks=seq(-2, 2, 2), expand=c(0, 0))+
  theme_classic2()
pv.pca.plot

# stat_conf_ellipse(data=pv.pca.df, aes(x=PC1, y=PC2, group=phas.veg$species.short), linetype="dashed")+
# summarize by accession - NOT USING
pv.pca.acsum <- pv.pca.df %>%
  group_by(phas.veg$lifespan, phas.veg$specific.epithet, phas.veg$accession) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(pv.pca.acsum)
names(pv.pca.acsum)[1] <- "lifespan"
names(pv.pca.acsum)[2] <- "species"
names(pv.pca.acsum)[3] <- "accession" 
pv.pca.acsum$species <- factor(pv.pca.acsum$species, levels=c("acutifolius", "filiformis", "lunatus", "vulgaris", "angustissimus", "leptostachyus", "maculatus", "parvulus", "polystachios"), ordered=T)

#### VICIA VEGETATIVE PCA

vicia.veg <- comb.plant.clean.new %>% filter(genus=="Vicia")
vicia.veg <- vicia.veg[!is.na(vicia.veg$height.perday.21.35),] # remove all rows where height AGR = NA, since this only co-occurs when one date has no data at all or all veg data are NA
vv <- vicia.veg %>% select(height.21.new, height.35.new, height.perday.21.35, height.perday.21.35.ln, leaves.21.new, leaves.35.new, leaves.perday.21.35, leaves.perday.21.35.ln)
colnames(vv) <- c("H DAP-21", "H DAP-35", "H AGR", "H RGR", "L DAP-21", "L DAP-35", "L AGR", "L RGR")
estim_ncpPCA(vv, scale=T) # estimates # of dimenstions needed in the next imputation step (ncp)
vv1 <- imputePCA(as.data.frame(vv), ncp=5, scale=T) # imputes missing values (NAs)
vv.pca <- prcomp(vv1$completeObs, scale=T)
get_eigenvalue(vv.pca) # Eigenvalues for PCs
fviz_eig(vv.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(vv.pca, col.var = "black", title="Vicia vegetative PCA variables")
vv.pca.loading <- vv.pca$rotation
vv.pca$rotation[,2] <- -vv.pca$rotation[,2] # flip to have PC space direction consistent between genera
vv.pca$x[,2] <- -vv.pca$x[,2]

# check imputation
write.xlsx(vicia.veg, file = "C:/Users/SterlingH/Desktop/R_exports/pre_vicia_veg_names1.xlsx") # For names
write.xlsx(vicia.veg, file = "C:/Users/SterlingH/Desktop/R_exports/vicia_veg_names1.xlsx") # For names
write.xlsx(vv, file = "C:/Users/SterlingH/Desktop/R_exports/vicia_veg_compare.xlsx")
write.xlsx(vv1, file = "C:/Users/SterlingH/Desktop/R_exports/vicia_veg_impute.xlsx")

# ggbiplot
vv.pca.varplot <- ggbiplot2(vv.pca, circle=T, alpha=0, varname.adjust=1.1,
                            varname.color="darkblue", arrow.color="black", circle.color="gray20")+
  geom_hline(yintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "gray", size=0.5)+
  labs(x="PC1 (53.6%)", y="PC2 (30.4%)", title="Vegetative variables")+
  scale_x_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  scale_y_continuous(limits=c(-2.5, 2.5), breaks=seq(-2, 2, 1), expand=c(0, 0))+
  theme_classic2()
vv.pca.varplot

# VICIA SPECIES VEGETATIVE PCA

vv.pca.df <- as.data.frame(vv.pca$x)

# summarize by species - USE
vv.pca.spsum <- vv.pca.df %>%
  group_by(vicia.veg$lifespan, vicia.veg$specific.epithet) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(vv.pca.spsum)
names(vv.pca.spsum)[1] <- "lifespan"
names(vv.pca.spsum)[2] <- "species"
vv.pca.spsum$species <- factor(vv.pca.spsum$species, levels=c("benghalensis", "ervilia", "hirsuta", "sativa", "villosa", "americana", "cracca", "dumetorum", "sepium", "tenuifolia"), ordered=T)

vv.pca.plot <- ggplot()+
  geom_hline(yintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  geom_vline(xintercept=0, linetype="dashed",
             color = "grey", size=0.5)+
  scale_shape_manual(name = "Species",
                     values = c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25))+
  scale_colour_manual(name = "Life span",
                      labels = c("annual", "perennial"),
                      values = c("red", "blue"))+
  geom_errorbar(data=vv.pca.spsum, aes(x=PC1.mn, ymin=PC2.mn-PC2.sd, ymax=PC2.mn+PC2.sd, color=lifespan), width=0.1, size=0.35)+
  geom_errorbarh(data=vv.pca.spsum, aes(y=PC2.mn, xmin=PC1.mn-PC1.sd, xmax=PC1.mn+PC1.sd, color=lifespan), height=0.1, size=0.35)+
  geom_point(data=vv.pca.spsum, aes(x=PC1.mn, y=PC2.mn, shape=species, fill=species), color="black", size=3)+
  scale_fill_manual(name = "Species",
                    values = c("#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))+
  labs(x="PC1 (53.6%)", y="PC2 (30.4%)", title="Vegetative PCA")+
  scale_x_continuous(limits=c(-3.5, 8), breaks=seq(-2, 8, 2), expand=c(0, 0))+
  scale_y_continuous(limits=c(-3.5, 2.5), breaks=seq(-2, 2, 2), expand=c(0, 0))+
  theme_classic2()
vv.pca.plot

#   stat_conf_ellipse(data=vv.pca.df, aes(x=PC1, y=PC2, group=vicia.veg$species.short), linetype="dashed")+
std <- function(x) sd(x)/sqrt(length(x)) # standard error function

# summarize by accession - NOT USING
vv.pca.acsum <- vv.pca.df %>%
  group_by(vicia.veg$lifespan, vicia.veg$specific.epithet, vicia.veg$accession) %>%
  summarise(
    PC1.mn = mean(PC1, na.rm=T),
    PC1.sd = sd(PC1, na.rm=T),
    PC2.mn = mean(PC2, na.rm=T),
    PC2.sd = sd(PC2, na.rm=T),
  ) %>% ungroup ()
View(vv.pca.acsum)
names(vv.pca.acsum)[1] <- "lifespan" 
names(vv.pca.acsum)[2] <- "species"
names(vv.pca.acsum)[3] <- "accession"
vv.pca.acsum$species <- factor(vv.pca.acsum$species, levels=c("benghalensis", "ervilia", "hirsuta", "sativa", "villosa", "americana", "cracca", "dumetorum", "sepium", "tenuifolia"), ordered=T)

# PCA PANELS ####

# Panel for Figure 1
ggpubr::ggarrange(a.pca.varplot, a.pca.genus, a.pca.lifespan, nrow=1, ncol=3, widths=c(2,2,2), heights=2,
                  legend="bottom", labels=c("A","B","C"), common.legend=T)
#labels=c("A  Variables","B  Genus","C  Life span")

# Seed & vegetative panels by genus
# Lathyrus
lath.pca <- ggpubr::ggarrange(ls.pca.varplot, ls.pca.plot, lv.pca.varplot, lv.pca.plot, nrow=2, ncol=2, widths=c(2,2), heights=c(2,2),
                  labels=c("A","B","C","D"), legend="bottom", common.legend=T)
annotate_figure(p=lath.pca, top=text_grob("Lathyrus", color ="black", face ="bold", size=15, hjust=5))

# Phaseolus
phas.pca <- ggpubr::ggarrange(ps.pca.varplot, ps.pca.plot, pv.pca.varplot, pv.pca.plot, nrow=2, ncol=2, widths=c(2,2), heights=c(2,2),
                  labels=c("A","B","C","D"), legend="bottom", common.legend=T)
annotate_figure(p=phas.pca, top=text_grob("Phaseolus", color ="black", face ="bold", size=15, hjust=4.2))

# Vicia
vic.pca <- ggpubr::ggarrange(vs.pca.varplot, vs.pca.plot, vv.pca.varplot, vv.pca.plot, nrow=2, ncol=2, widths=c(2,2), heights=c(2,2),
                  labels=c("A","B","C","D"), legend="bottom", common.legend=T)
annotate_figure(p=vic.pca, top=text_grob("Vicia", color ="black", face ="bold", size=15, hjust=9))

# OLD

# Panel of seed & vegetative PCA original
ggpubr::ggarrange(ls.pca.plot, ps.pca.plot, vs.pca.plot, lv.pca.plot, pv.pca.plot, vs.pca.plot, 
                  nrow=2, ncol=3, widths=c(2,2,2), heights=c(2,2),
                  labels=c("A","B","C","D","E","F"), legend="right")

# Grid for additional PCs
ggpubr::ggarrange(pc13, pc14, pc23, pc24, pc34, nrow=3, ncol=2, widths=c(2,2), heights=c(2,2,2),
                  labels=c("A","B","C","D","E"), legend="bottom", common.legend=T)

#######################################################

# Raw species violin plots - trait by trait ####

# Seed data
lath.seed.raw <- comb.seed.clean %>% filter(genus=="Lathyrus")
phas.seed.raw <- comb.seed.clean %>% filter(genus=="Phaseolus")
vicia.seed.raw <- comb.seed.clean %>% filter(genus=="Vicia")

# Germination data
lath.germ.raw <- comb.germ.clean %>% filter(genus=="Lathyrus")
phas.germ.raw <- comb.germ.clean %>% filter(genus=="Phaseolus")
vicia.germ.raw <- comb.germ.clean %>% filter(genus=="Vicia")

# Vegetative growth data
lath.veg.raw <- comb.plant.clean.new %>% filter(genus=="Lathyrus")
phas.veg.raw <- comb.plant.clean.new %>% filter(genus=="Phaseolus")
vicia.veg.raw <- comb.plant.clean.new %>% filter(genus=="Vicia")

# Accession level data
lath.acc1 <- merge.acc.all %>% filter(genus=="Lathyrus")
phas.acc1 <- merge.acc.all %>% filter(genus=="Phaseolus")
vicia.acc1 <- merge.acc.all %>% filter(genus=="Vicia")
ann.acc1 <- merge.acc.all %>% filter(lifespan=="annual")
per.acc1 <- merge.acc.all %>% filter(lifespan=="perennial")

# order species by lifespan alphabetically

lath.seed.raw$specific.epithet <- factor(lath.seed.raw$specific.epithet, levels=c("annuus","aphaca","cicera","hirsutus","odoratus","japonicus","latifolius","pratensis","sylvestris","tuberosus"), ordered=T)
phas.seed.raw$specific.epithet <- factor(phas.seed.raw$specific.epithet, levels=c("acutifolius","filiformis","lunatus","vulgaris","angustissimus","leptostachyus","maculatus","parvulus","polystachios"), ordered=T)
vicia.seed.raw$specific.epithet <- factor(vicia.seed.raw$specific.epithet, levels=c("benghalensis","ervilia","hirsuta","sativa","villosa","americana","cracca","dumetorum","sepium","tenuifolia"), ordered=T)


lath.germ.raw$specific.epithet <- factor(lath.germ.raw$specific.epithet, levels=c("annuus","aphaca","cicera","hirsutus","odoratus","japonicus","latifolius","pratensis","sylvestris","tuberosus"), ordered=T)
phas.germ.raw$specific.epithet <- factor(phas.germ.raw$specific.epithet, levels=c("acutifolius","filiformis","lunatus","vulgaris","angustissimus","leptostachyus","maculatus","parvulus","polystachios"), ordered=T)
vicia.germ.raw$specific.epithet <- factor(vicia.germ.raw$specific.epithet, levels=c("benghalensis","ervilia","hirsuta","sativa","villosa","americana","cracca","dumetorum","sepium","tenuifolia"), ordered=T)


lath.veg.raw$specific.epithet <- factor(lath.veg.raw$specific.epithet, levels=c("annuus","aphaca","cicera","hirsutus","odoratus","japonicus","latifolius","pratensis","sylvestris","tuberosus"), ordered=T)
phas.veg.raw$specific.epithet <- factor(phas.veg.raw$specific.epithet, levels=c("acutifolius","filiformis","lunatus","vulgaris","angustissimus","leptostachyus","maculatus","parvulus","polystachios"), ordered=T)
vicia.veg.raw$specific.epithet <- factor(vicia.veg.raw$specific.epithet, levels=c("benghalensis","ervilia","hirsuta","sativa","villosa","americana","cracca","dumetorum","sepium","tenuifolia"), ordered=T)


lath.acc1$specific.epithet <- factor(lath.acc1$specific.epithet, levels=c("annuus","aphaca","cicera","hirsutus","odoratus","japonicus","latifolius","pratensis","sylvestris","tuberosus"), ordered=T)
phas.acc1$specific.epithet <- factor(phas.acc1$specific.epithet, levels=c("acutifolius","filiformis","lunatus","vulgaris","angustissimus","leptostachyus","maculatus","parvulus","polystachios"), ordered=T)
vicia.acc1$specific.epithet <- factor(vicia.acc1$specific.epithet, levels=c("benghalensis","ervilia","hirsuta","sativa","villosa","americana","cracca","dumetorum","sepium","tenuifolia"), ordered=T)

# to group by mathematical value
# x=reorder(specific.epithet, seed.wt.avg, FUN=median)

# Seed mass

# Lathyrus
lath.sw.boxplot <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=specific.epithet, y=seed.wt.avg, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=seed.wt.avg, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Averaged seed mass (mg)", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.sw.boxplot

lath.sw.boxplot.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=seed.wt.avg, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=seed.wt.avg, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Averaged seed mass (mg)", title="Lathyrus")+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.sw.boxplot.ap

# Phaseolus
phas.sw.boxplot <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=specific.epithet, y=seed.wt.avg, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=seed.wt.avg, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Averaged seed mass (mg)", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.sw.boxplot

phas.sw.boxplot.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=seed.wt.avg, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=seed.wt.avg, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Averaged seed mass (mg)", title="Phaseolus")+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.sw.boxplot.ap

# Vicia
vicia.sw.boxplot <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=specific.epithet, y=seed.wt.avg, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=seed.wt.avg, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Averaged seed mass (mg)", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.sw.boxplot

vicia.sw.boxplot.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=seed.wt.avg, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=seed.wt.avg, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Averaged seed mass (mg)", title="Vicia")+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.sw.boxplot.ap

# Seed length
# Lathyrus
lath.sl.violin <- ggplot()+
  geom_violin(data=lath.seed.raw, aes(x=specific.epithet, y=seed.length, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.seed.raw, aes(x=specific.epithet, y=seed.length), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=seed.length.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed length (mm)", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.sl.violin

lath.sl.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=seed.length.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=seed.length.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.sl.ap

# Phaseolus
phas.sl.violin <- ggplot()+
  geom_violin(data=phas.seed.raw, aes(x=specific.epithet, y=seed.length, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.seed.raw, aes(x=specific.epithet, y=seed.length), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=seed.length.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed length (mm)", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.sl.violin

phas.sl.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=seed.length.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=seed.length.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.sl.ap

# Vicia
vicia.sl.violin <- ggplot()+
  geom_violin(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.length, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.length), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=seed.length.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed length (mm)", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.sl.violin

vicia.sl.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=seed.length.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=seed.length.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.sl.ap

# Seed width
# Lathyrus
lath.sw.violin <- ggplot()+
  geom_violin(data=lath.seed.raw, aes(x=specific.epithet, y=seed.width, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.seed.raw, aes(x=specific.epithet, y=seed.width), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=seed.width.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed width (mm)", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.sw.violin

lath.sw.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=seed.width.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=seed.width.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.sw.ap

# Phaseolus
phas.sw.violin <- ggplot()+
  geom_violin(data=phas.seed.raw, aes(x=specific.epithet, y=seed.width, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.seed.raw, aes(x=specific.epithet, y=seed.width), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=seed.width.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed width (mm)", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.sw.violin

phas.sw.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=seed.width.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=seed.width.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.sw.ap

# Vicia
vicia.sw.violin <- ggplot()+
  geom_violin(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.width, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.width), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=seed.width.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed width (mm)", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.sw.violin

vicia.sw.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=seed.width.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=seed.width.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.sw.ap

# Seed perimeter
# Lathyrus
lath.sp.violin <- ggplot()+
  geom_violin(data=lath.seed.raw, aes(x=specific.epithet, y=seed.perim, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.seed.raw, aes(x=specific.epithet, y=seed.perim), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=seed.perim.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed perimeter (mm)", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.sp.violin

lath.sp.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=seed.perim.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=seed.perim.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.sp.ap

# Phaseolus
phas.sp.violin <- ggplot()+
  geom_violin(data=phas.seed.raw, aes(x=specific.epithet, y=seed.perim, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.seed.raw, aes(x=specific.epithet, y=seed.perim), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=seed.perim.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed perimeter (mm)", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.sp.violin

phas.sp.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=seed.perim.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=seed.perim.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.sp.ap

# Vicia
vicia.sp.violin <- ggplot()+
  geom_violin(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.perim, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.perim), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=seed.perim.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed perimeter (mm)", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.sp.violin

vicia.sp.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=seed.perim.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=seed.perim.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.sp.ap

# Seed area
# Lathyrus
lath.sa.violin <- ggplot()+
  geom_violin(data=lath.seed.raw, aes(x=specific.epithet, y=seed.area, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.seed.raw, aes(x=specific.epithet, y=seed.area), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=seed.area.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed area (mm^2)", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.sa.violin

lath.sa.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=seed.area.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=seed.area.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.sa.ap

# Phaseolus
phas.sa.violin <- ggplot()+
  geom_violin(data=phas.seed.raw, aes(x=specific.epithet, y=seed.area, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.seed.raw, aes(x=specific.epithet, y=seed.area), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=seed.area.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed area (mm^2)", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.sa.violin

phas.sa.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=seed.area.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=seed.area.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.sa.ap

# Vicia
vicia.sa.violin <- ggplot()+
  geom_violin(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.area, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.area), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=seed.area.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed area (mm^2)", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.sa.violin

vicia.sa.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=seed.area.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=seed.area.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.sa.ap

# Seed circularity
# Lathyrus
lath.sc.violin <- ggplot()+
  geom_violin(data=lath.seed.raw, aes(x=specific.epithet, y=seed.circ, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.seed.raw, aes(x=specific.epithet, y=seed.circ), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=seed.circ.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed circularity", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.sc.violin

lath.sc.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=seed.circ.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=seed.circ.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.sc.ap

# Phaseolus
phas.sc.violin <- ggplot()+
  geom_violin(data=phas.seed.raw, aes(x=specific.epithet, y=seed.circ, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.seed.raw, aes(x=specific.epithet, y=seed.circ), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=seed.circ.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed circularity", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.sc.violin

phas.sc.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=seed.circ.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=seed.circ.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.sc.ap

# Vicia
vicia.sc.violin <- ggplot()+
  geom_violin(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.circ, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.circ), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=seed.circ.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed circularity", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.sc.violin

vicia.sc.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=seed.circ.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=seed.circ.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.sc.ap

# Seed roundness

# Lathyrus
lath.sr.violin <- ggplot()+
  geom_violin(data=lath.seed.raw, aes(x=specific.epithet, y=seed.roundness, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.seed.raw, aes(x=specific.epithet, y=seed.roundness), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=seed.round.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed roundness", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.sr.violin

lath.sr.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=seed.round.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=seed.round.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.sr.ap

# Phaseolus
phas.sr.violin <- ggplot()+
  geom_violin(data=phas.seed.raw, aes(x=specific.epithet, y=seed.roundness, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.seed.raw, aes(x=specific.epithet, y=seed.roundness), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=seed.round.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed roundness", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.sr.violin

phas.sr.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=seed.round.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=seed.round.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.sr.ap

# Vicia
vicia.sr.violin <- ggplot()+
  geom_violin(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.roundness, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.seed.raw, aes(x=specific.epithet, y=seed.roundness), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=seed.round.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Seed roundness", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.sr.violin

vicia.sr.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=seed.round.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=seed.round.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.sr.ap

# Germination T50

# Lathyrus
lath.gt.boxplot <- ggplot()+
  geom_boxplot(data=lath.germ.raw, aes(x=specific.epithet, y=ppd50, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.germ.raw, aes(x=specific.epithet, y=ppd50, fill=accession), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  labs(x="Species",y="Germination T50", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.gt.boxplot

lath.gt.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=ppd50.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=ppd50.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.gt.ap

# Phaseolus
phas.gt.boxplot <- ggplot()+
  geom_boxplot(data=phas.germ.raw, aes(x=specific.epithet, y=ppd50, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.germ.raw, aes(x=specific.epithet, y=ppd50, fill=accession), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  labs(x="Species",y="Germination T50", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.gt.boxplot

phas.gt.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=ppd50.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=ppd50.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.gt.ap

# Vicia
vicia.gt.boxplot <- ggplot()+
  geom_boxplot(data=vicia.germ.raw, aes(x=specific.epithet, y=ppd50, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.germ.raw, aes(x=specific.epithet, y=ppd50, fill=accession), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  labs(x="Species",y="Germination T50", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.gt.boxplot

vicia.gt.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=ppd50.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=ppd50.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.gt.ap
  
# Germination proportion

# Lathyrus
lath.gp.boxplot <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=specific.epithet, y=FINAL.germ.prop, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=FINAL.germ.prop, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Germination proportion", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.gp.boxplot

lath.gp.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=FINAL.germ.prop, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=FINAL.germ.prop, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.gp.ap

# Phaseolus
phas.gp.boxplot <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=specific.epithet, y=FINAL.germ.prop, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=FINAL.germ.prop, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Germination proportion", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.gp.boxplot

phas.gp.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=FINAL.germ.prop, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=FINAL.germ.prop, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.gp.ap

# Vicia
vicia.gp.boxplot <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=specific.epithet, y=FINAL.germ.prop, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=FINAL.germ.prop, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Germination proportion", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.gp.boxplot

vicia.gp.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=FINAL.germ.prop, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=FINAL.germ.prop, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.gp.ap

# Height raw growth

# Lathyrus
lath.hg.violin <- ggplot()+
  geom_violin(data=lath.veg.raw, aes(x=specific.epithet, y=height.perday.21.35, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.veg.raw, aes(x=specific.epithet, y=height.perday.21.35), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=height.perday.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height raw growth rate (mm/day)", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.hg.violin

lath.agr.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=height.perday.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=height.perday.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.agr.ap

# Phaseolus
phas.hg.violin <- ggplot()+
  geom_violin(data=phas.veg.raw, aes(x=specific.epithet, y=height.perday.21.35, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.veg.raw, aes(x=specific.epithet, y=height.perday.21.35), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=height.perday.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height raw growth rate (mm/day)", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.hg.violin

phas.agr.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=height.perday.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=height.perday.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.agr.ap

# Vicia
vicia.hg.violin <- ggplot()+
  geom_violin(data=vicia.veg.raw, aes(x=specific.epithet, y=height.perday.21.35, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.veg.raw, aes(x=specific.epithet, y=height.perday.21.35), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=height.perday.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height raw growth rate (mm/day)", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.hg.violin

vicia.agr.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=height.perday.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=height.perday.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.agr.ap

# Height RGR

# Lathyrus
lath.hrgr.violin <- ggplot()+
  geom_violin(data=lath.veg.raw, aes(x=specific.epithet, y=height.perday.21.35.ln, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.veg.raw, aes(x=specific.epithet, y=height.perday.21.35.ln), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=height.rgr.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height relative growth rate", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.hrgr.violin

lath.rgr.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=height.rgr.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=height.rgr.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.rgr.ap

# Phaseolus
phas.hrgr.violin <- ggplot()+
  geom_violin(data=phas.veg.raw, aes(x=specific.epithet, y=height.perday.21.35.ln, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.veg.raw, aes(x=specific.epithet, y=height.perday.21.35.ln), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=height.rgr.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height relative growth rate", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.hrgr.violin

phas.rgr.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=height.rgr.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=height.rgr.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.rgr.ap

# Vicia
vicia.hrgr.violin <- ggplot()+
  geom_violin(data=vicia.veg.raw, aes(x=specific.epithet, y=height.perday.21.35.ln, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.veg.raw, aes(x=specific.epithet, y=height.perday.21.35.ln), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=height.rgr.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height relative growth rate", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.hrgr.violin

vicia.rgr.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=height.rgr.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=height.rgr.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.rgr.ap

# Leaf raw growth

# Lathyrus
lath.lg.violin <- ggplot()+
  geom_violin(data=lath.veg.raw, aes(x=specific.epithet, y=leaves.perday.21.35, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.veg.raw, aes(x=specific.epithet, y=leaves.perday.21.35), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=leaves.perday.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaf growth rate (leaves/day)", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.lg.violin

lath.lagr.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=leaves.perday.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=leaves.perday.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.lagr.ap

# Phaseolus
phas.lg.violin <- ggplot()+
  geom_violin(data=phas.veg.raw, aes(x=specific.epithet, y=leaves.perday.21.35, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.veg.raw, aes(x=specific.epithet, y=leaves.perday.21.35), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=leaves.perday.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaf growth rate (leaves/day)", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.lg.violin

phas.lagr.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=leaves.perday.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=leaves.perday.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.lagr.ap

# Vicia
vicia.lg.violin <- ggplot()+
  geom_violin(data=vicia.veg.raw, aes(x=specific.epithet, y=leaves.perday.21.35, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.veg.raw, aes(x=specific.epithet, y=leaves.perday.21.35), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=leaves.perday.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaf growth rate (leaves/day)", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.lg.violin

vicia.lagr.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=leaves.perday.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=leaves.perday.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.lagr.ap

# Height day 21

# Lathyrus
lath.h21.violin <- ggplot()+
  geom_violin(data=lath.veg.raw, aes(x=specific.epithet, y=height.21.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.veg.raw, aes(x=specific.epithet, y=height.21.new), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=height.21.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height DAP-21", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.h21.violin

lath.h21.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=height.21.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=height.21.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.h21.ap

# Phaseolus
phas.h21.violin <- ggplot()+
  geom_violin(data=phas.veg.raw, aes(x=specific.epithet, y=height.21.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.veg.raw, aes(x=specific.epithet, y=height.21.new), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=height.21.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height DAP-21", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.h21.violin

phas.h21.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=height.21.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=height.21.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.h21.ap

# Vicia
vicia.h21.violin <- ggplot()+
  geom_violin(data=vicia.veg.raw, aes(x=specific.epithet, y=height.21.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.veg.raw, aes(x=specific.epithet, y=height.21.new), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=height.21.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height DAP-21", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.h21.violin

vicia.h21.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=height.21.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=height.21.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.h21.ap

# Height day 35

# Lathyrus
lath.h35.violin <- ggplot()+
  geom_violin(data=lath.veg.raw, aes(x=specific.epithet, y=height.35.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.veg.raw, aes(x=specific.epithet, y=height.35.new), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=height.35.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height DAP-35", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.h35.violin

lath.h35.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=height.35.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=height.35.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.h35.ap

# Phaseolus
phas.h35.violin <- ggplot()+
  geom_violin(data=phas.veg.raw, aes(x=specific.epithet, y=height.35.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.veg.raw, aes(x=specific.epithet, y=height.35.new), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=height.35.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height DAP-35", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.h35.violin

phas.h35.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=height.35.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=height.35.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.h35.ap

# Vicia
vicia.h35.violin <- ggplot()+
  geom_violin(data=vicia.veg.raw, aes(x=specific.epithet, y=height.35.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.veg.raw, aes(x=specific.epithet, y=height.35.new), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=height.35.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Height DAP-35", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.h35.violin

vicia.h35.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=height.35.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=height.35.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.h35.ap

# Leaves day 21

# Lathyrus
lath.l21.violin <- ggplot()+
  geom_violin(data=lath.veg.raw, aes(x=specific.epithet, y=leaves.21.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.veg.raw, aes(x=specific.epithet, y=leaves.21.new), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=leaves.21.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaves DAP-21", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.l21.violin

lath.l21.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=leaves.21.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=leaves.21.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.l21.ap

# Phaseolus
phas.l21.violin <- ggplot()+
  geom_violin(data=phas.veg.raw, aes(x=specific.epithet, y=leaves.21.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.veg.raw, aes(x=specific.epithet, y=leaves.21.new), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=leaves.21.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaves DAP-21", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.l21.violin

phas.l21.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=leaves.21.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=leaves.21.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.l21.ap

# Vicia
vicia.l21.violin <- ggplot()+
  geom_violin(data=vicia.veg.raw, aes(x=specific.epithet, y=leaves.21.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.veg.raw, aes(x=specific.epithet, y=leaves.21.new), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=leaves.21.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaves DAP-21", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.l21.violin

vicia.l21.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=leaves.21.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=leaves.21.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.l21.ap

# Leaves day 35

# Lathyrus
lath.l35.violin <- ggplot()+
  geom_violin(data=lath.veg.raw, aes(x=specific.epithet, y=leaves.35.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=lath.veg.raw, aes(x=specific.epithet, y=leaves.35.new), color="black", fill="white", width=0.1)+
  geom_point(data=lath.acc1, aes(x=specific.epithet, y=leaves.35.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaves DAP-35", title="Lathyrus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
lath.l35.violin

lath.l35.ap <- ggplot()+
  geom_boxplot(data=lath.acc1, aes(x=lifespan, y=leaves.35.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=lath.acc1, aes(x=lifespan, y=leaves.35.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
lath.l35.ap

# Phaseolus
phas.l35.violin <- ggplot()+
  geom_violin(data=phas.veg.raw, aes(x=specific.epithet, y=leaves.35.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=phas.veg.raw, aes(x=specific.epithet, y=leaves.35.new), color="black", fill="white", width=0.1)+
  geom_point(data=phas.acc1, aes(x=specific.epithet, y=leaves.35.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaves DAP-35", title="Phaseolus")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
phas.l35.violin

phas.l35.ap <- ggplot()+
  geom_boxplot(data=phas.acc1, aes(x=lifespan, y=leaves.35.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=phas.acc1, aes(x=lifespan, y=leaves.35.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
phas.l35.ap

# Vicia
vicia.l35.violin <- ggplot()+
  geom_violin(data=vicia.veg.raw, aes(x=specific.epithet, y=leaves.35.new, fill=lifespan), color="black")+
  scale_fill_manual(values = c("firebrick3","dodgerblue4"))+
  geom_boxplot(data=vicia.veg.raw, aes(x=specific.epithet, y=leaves.35.new), color="black", fill="white", width=0.1)+
  geom_point(data=vicia.acc1, aes(x=specific.epithet, y=leaves.35.mn, color=lifespan),
             pch=16, size=2, alpha=0.8, show.legend = F, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))+
  scale_colour_manual(values = c("palevioletred1","deepskyblue"))+
  labs(x="Species",y="Leaves DAP-35", title="Vicia")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) # has to come after theme_classic or the changes are overridden
vicia.l35.violin

vicia.l35.ap <- ggplot()+
  geom_boxplot(data=vicia.acc1, aes(x=lifespan, y=leaves.35.mn, color=lifespan))+
  scale_colour_manual(values = c("firebrick3","dodgerblue4"))+
  geom_point(data=vicia.acc1, aes(x=lifespan, y=leaves.35.mn, fill=lifespan), position=position_jitterdodge(jitter.width=0.5, dodge.width = 0),
             pch=21, size=2, show.legend = F, alpha=0.6)+
  scale_fill_manual(values = c("palevioletred1","deepskyblue"))+
  theme_classic2() # has to come after theme_classic or the changes are overridden
vicia.l35.ap

#######################################################

# LM PREP ####

# Make sure your contrasts are set correctly!!
options("contrasts") # check
options(contrasts = c("contr.treament", "contr.poly")) # this should be the default; required before running a Type III SS with Anova function
options(contrasts = c("contr.sum", "contr.poly"))

# Create vectors to represent each mean to be compared for Tukey comparisons
# Make sure this list is in the same order as the output of mod.emm (emmeans table) for each model
LathA <- c(1,0,0,0,0,0)
LathP <- c(0,1,0,0,0,0)
PhasA <- c(0,0,1,0,0,0)
PhasP <- c(0,0,0,1,0,0)
ViciA <- c(0,0,0,0,1,0)
ViciP <- c(0,0,0,0,0,1)

# Remove species for which there are only one accession:
no.single <- comb.acc.clean[!comb.acc.clean$species.short %in% names(which(table(comb.acc.clean$species.short) == 1)), ]

# dropping data
merge.lunatus <- merge.acc.all[!merge.acc.all$species.short == "P. lunatus",] # for seed mass and germination prop

# LM - PCs ####

# PC1 model
pc1.model <- lm(data=a.pca.df, PC1~genus + lifespan + genus*lifespan + genus*lifespan/species.short)

anova(pc1.model) # default is type I for lm
Anova(pc1.model, type="III", singular.ok=T)
pc1.mod.emm <- emmeans(object=pc1.model, specs=c("lifespan","genus"))
pc1.df <- as.data.frame(pc1.mod.emm)

contrast(pc1.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# PC2 model
pc2.model <- lm(data=a.pca.df, PC2~genus + lifespan + genus*lifespan + genus*lifespan/species.short)

anova(pc2.model) # default is type I for lm
Anova(pc2.model, type="II")
pc2.mod.emm <- emmeans(object=pc2.model, specs=c("lifespan","genus"))
pc2.df <- as.data.frame(pc2.mod.emm)

contrast(pc2.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# LM - seed ####

# dropping data
seed.lunatus <- comb.seed.clean[!comb.seed.clean$species.short == "P. lunatus",]

# SEED mass MODEL

# test data
test.all <- read_excel("C:/Users/SterlingH/Desktop/test_lm.xlsx", 
                       na = "NA", sheet="full")
test.nosingle <- read_excel("C:/Users/SterlingH/Desktop/test_lm.xlsx", 
                       na = "NA", sheet="singles.removed")

sdwt.model <- lm(data=merge.acc.all, seed.wt.avg~genus + lifespan + genus*lifespan + genus*lifespan/species.short)

anova(sdwt.model)
summary(sdwt.model)
summary(sdwt.model)$fstatistic
Anova(sdwt.model, type=3, singular.ok = T)
alias <- alias(sdwt.model)
sdwt.mod.emm <- emmeans(object=sdwt.model, specs=c("lifespan","genus"))
sdwt.df <- as.data.frame(sdwt.mod.emm)

contrast(sdwt.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# test which variables are linearly dependent (aliased coefficients)
alias(model) # check aliased variables

Anova(lm(data=comb.acc.clean, seed.wt.avg~genus + lifespan
         + genus*lifespan + genus*lifespan/species.short), type="III", singular.ok=T)

Anova(lm(data=comb.acc.clean, seed.wt.avg~genus + lifespan + genus*lifespan),type="III",singular.ok=T)

Anova(lm(data=comb.acc.clean, seed.wt.avg~genus*lifespan, contrasts=list(genus=contr.sum, lifespan=contr.sum),type="III"))

Anova(lm(data=comb.acc.clean, seed.wt.avg~genus + lifespan + genus*lifespan
   + genus/lifespan/species.short))

# SEED mass Tukey
av <- aov(lm(data=comb.acc.clean, seed.wt.avg~genus + lifespan + genus*lifespan
             + genus*lifespan/species.short))
TukeyHSD(av, which="genus:lifespan")

# SEED LENGTH MODEL
sdlg.model <- lmer(data=comb.seed.clean, seed.length~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession))

r <- rand(sdlg.model)
anova(sdlg.model)
sdlg.mod.emm <- emmeans(object=sdlg.model, specs=c("lifespan","genus"))
sdlg.df <- as.data.frame(sdlg.mod.emm)


contrast(sdlg.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# SEED WIDTH MODEL
sdwd.model <- lmer(data=comb.seed.clean, seed.width~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession))

r <- rand(sdwd.model)
anova(sdwd.model)
sdwd.mod.emm <- emmeans(object=sdwd.model, specs=c("lifespan","genus"))
sdwd.df <- as.data.frame(sdwd.mod.emm)

contrast(sdwd.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# SEED PERIMETER MODEL
sdpm.model <- lmer(data=comb.seed.clean, seed.perim~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession))

r <- rand(sdpm.model)
anova(sdpm.model)
sdpm.mod.emm <- emmeans(object=sdpm.model, specs=c("lifespan","genus"))
sdpm.df <- as.data.frame(sdpm.mod.emm)

contrast(sdpm.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# SEED AREA MODEL
sdar.model <- lmer(data=comb.seed.clean, seed.area~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession))

r <- rand(sdar.model)
anova(sdar.model)
sdar.mod.emm <- emmeans(object=sdar.model, specs=c("lifespan","genus"))
sdar.df <- as.data.frame(sdar.mod.emm)

contrast(sdar.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# SEED CIRCULARITY
sdcr.model <- lmer(data=comb.seed.clean, seed.circ~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession))

r <- rand(sdcr.model)
anova(sdcr.model)
sdcr.mod.emm <- emmeans(object=sdcr.model, specs=c("lifespan","genus"))
sdcr.df <- as.data.frame(sdcr.mod.emm)

contrast(sdcr.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")


# SEED ROUNDNESS
sdrd.model <- lmer(data=comb.seed.clean, seed.roundness~genus + lifespan + genus*lifespan + genus*lifespan/species.short
              + (1 | species.short:accession))

r <- rand(sdrd.model)
anova(sdrd.model)
sdrd.mod.emm <- emmeans(object=sdrd.model, specs=c("lifespan","genus"))
sdrd.df <- as.data.frame(sdrd.mod.emm)

contrast(sdrd.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# LM - germination ####

# dropping accession where half the seeds were destroyed from bleach
germrep.vhi <- comb.germ.clean[!comb.germ.clean$accession == "PI 219631",] 
germacc.vhi <- merge.acc.all[!merge.acc.all$accession == "PI 219631",]
# excluding accessions with a lot of ungerminated, imbibed, viable seeds
# viab <- c("W6 2427", "PI 494749")
# germrep.viab <- comb.germ.clean[!comb.germ.clean$accession %in% viab,]
# germacc.viab <- merge.acc.all[!merge.acc.all$accession %in% viab,]
germrt.lunatus <- comb.germ.clean[!comb.germ.clean$species.short == "P. lunatus",]

### GERMINATION T50
# normal models use comb.germ.clean (T50)
# merge.acc.all (from acc.germ.clean) (proportion)

# Germination T50
# ORIGINAL
model <- lmer(data=comb.germ.clean, ppd50~genus + lifespan + genus*lifespan + genus*lifespan/species.short
              + (1 | species.short:accession) + (1 | age))
# FINAL
# Nonsignificant: age (marginal), sandpaper
gmrt.modreduced <- lmer(data=comb.germ.clean, ppd50~genus + lifespan + genus*lifespan + genus*lifespan/species.short
                             + (1 | species.short:accession))
# put in age as a random effect if using germrep.vhi

rand(model)
r <- rand(gmrt.modreduced)
anova(gmrt.modreduced)
gmrt.mod.emm <- emmeans(object=gmrt.modreduced, specs=c("lifespan","genus"))
gmrt.df <- as.data.frame(gmrt.mod.emm)

contrast(gmrt.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                                   "PhasA-PhasP" = PhasA-PhasP,
                                   "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")


### GERMINATION PROPORTION

# ORIGINAL
model <- lmer(data=merge.acc.all, FINAL.germ.prop~genus + lifespan + genus*lifespan + genus*lifespan/species.short
              + (1 | age))

# FINAL 
# Nonsignificant: age 
gmprop.modreduced <- lm(data=merge.acc.all, FINAL.germ.prop~genus + lifespan + genus*lifespan + genus*lifespan/species.short)

rand(model)
rand(gmprop.modreduced)
anova(gmprop.modreduced)
gmprop.mod.emm <- emmeans(object=gmprop.modreduced, specs=c("lifespan","genus"))
gmprop.df <- as.data.frame(gmprop.mod.emm)

contrast(gmprop.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                                   "PhasA-PhasP" = PhasA-PhasP,
                                   "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# LM - vegetative ####

# Dropping data

comb.plant.fertzero <- comb.plant.clean.new[!comb.plant.clean.new$fert.zero == "yes",]
comb.plant.fert13 <- comb.plant.clean.new[!comb.plant.clean.new$fert.13 == "yes",]
comb.plant.grnhs <- comb.plant.clean.new[!comb.plant.clean.new$grnhs.2 == "yes",]
comb.plant.lunatus <- comb.plant.clean.new[!comb.plant.clean.new$species.short == "P. lunatus",]

# taking out days from fertilization from all models!
# putting in height.21.new as a control for size as opposed to days from sowing for growth rate; still days from sowing for static growth

# HEIGHT AGR
# ORIGINAL

model <- lmer(data=comb.plant.clean.new, height.perday.21.35~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.35.new) + (1 | repro.21.35.new) + (1 | height.21.new))
# FINAL
# nonsig: height
htagr.model.reduced <- lmer(data=comb.plant.clean.new, height.perday.21.35~genus + lifespan + genus*lifespan + genus*lifespan/species.short
          + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.35.new) + (1 | repro.21.35.new))

rand(model)
r <- rand(htagr.model.reduced)
anova(htagr.model.reduced)
htagr.mod.emm <- emmeans(object=htagr.model.reduced, specs=c("lifespan","genus"))
htagr.df <- as.data.frame(htagr.mod.emm)

contrast(htagr.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# HEIGHT RGR
# ORIGINAL

model <- lmer(data=comb.plant.clean.new, height.perday.21.35.ln~genus + lifespan + genus*lifespan + genus*lifespan/species.short
              + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.35.new) + (1 | repro.21.35.new) + (1 | height.21.new))
# FINAL
# nonsig: vigor, repro
htrgr.model.reduced <- lmer(data=comb.plant.clean.new, height.perday.21.35.ln~genus + lifespan + genus*lifespan + genus*lifespan/species.short
                           + (1 | species.short:accession) + (1 | accession:replicate) + (1 | height.21.new))

rand(model)
r <- rand(htrgr.model.reduced)
anova(htrgr.model.reduced)
htrgr.mod.emm <- emmeans(object=htrgr.model.reduced, specs=c("lifespan","genus"))
htrgr.df <- as.data.frame(htrgr.mod.emm)

contrast(htrgr.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                                   "PhasA-PhasP" = PhasA-PhasP,
                                   "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")


# LEAF AGR
# ORIGINAL
model <- lmer(data=comb.plant.clean.new, leaves.perday.21.35~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.35.new) + (1 | repro.21.35.new) + (1 | height.21.new))

# FINAL
# Nonsignificant: height
lfagr.model.reduced <- lmer(data=comb.plant.clean.new, leaves.perday.21.35~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.35.new) + (1 | repro.21.35.new))

rand(model)
r <- rand(lfagr.model.reduced)
anova(lfagr.model.reduced)
lfagr.mod.emm <- emmeans(object=lfagr.model.reduced, specs=c("lifespan","genus"))
lfagr.df <- as.data.frame(lfagr.mod.emm)

contrast(lfagr.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# LEAF RGR
# ORIGINAL

model <- lmer(data=comb.plant.clean.new, leaves.perday.21.35.ln~genus + lifespan + genus*lifespan + genus*lifespan/species.short
              + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.35.new) + (1 | repro.21.35.new) + (1 | height.21.new))
# FINAL
# nonsig: vigor, repro, height 
lfrgr.model.reduced <- lmer(data=comb.plant.clean.new, leaves.perday.21.35.ln~genus + lifespan + genus*lifespan + genus*lifespan/species.short
                            + (1 | species.short:accession) + (1 | accession:replicate))

rand(model)
r <- rand(lfrgr.model.reduced)
anova(lfrgr.model.reduced)
lfrgr.mod.emm <- emmeans(object=lfrgr.model.reduced, specs=c("lifespan","genus"))
lfrgr.df <- as.data.frame(lfrgr.mod.emm)

contrast(lfrgr.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                                    "PhasA-PhasP" = PhasA-PhasP,
                                    "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# HEIGHT at 21 DAYS
# ORIGINAL
# originally had random effect + (1 | sowing.to.planting), but this was not able to be determined accurately

model <- lmer(data=comb.plant.clean.new, height.21.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
     + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.new) + (1 | repro.21.new))

# FINAL
# Nonsignificant: repro
ht21.model.reduced <- lmer(data=comb.plant.clean.new, height.21.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
                      + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.new))

rand(model)
r <- rand(ht21.model.reduced)
anova(ht21.model.reduced)
ht21.mod.emm <- emmeans(object=ht21.model.reduced, specs=c("lifespan","genus"))
ht21.df <- as.data.frame(ht21.mod.emm)

contrast(ht21.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
                              adjust="bonferroni")

# LEAVES at 21 DAYS
# ORIGINAL
model <- lmer(data=comb.plant.clean.new, leaves.21.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
              + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.new) + (1 | repro.21.new))

# FINAL
# Nonsignificant: repro
lf21.model.reduced <- lmer(data=comb.plant.clean.new, leaves.21.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
                      + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.21.new))

rand(model)
r <- rand(lf21.model.reduced)
anova(lf21.model.reduced)
lf21.mod.emm <- emmeans(object=lf21.model.reduced, specs=c("lifespan","genus"))
lf21.df <- as.data.frame(lf21.mod.emm)

contrast(lf21.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# HEIGHT at 35 DAYS

# ORIGINAL
model <- lmer(data=comb.plant.clean.new, height.35.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.35.new) + (1 | repro.35.new))

# FINAL
# Nonsignificant: all significant.
ht35.model.reduced <- lmer(data=comb.plant.clean.new, height.35.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
     + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.35.new) + (1 | repro.35.new))

rand(model)
r <- rand(ht35.model.reduced)
anova(ht35.model.reduced)
ht35.mod.emm <- emmeans(object=ht35.model.reduced, specs=c("lifespan","genus"))
ht35.df <- as.data.frame(ht35.mod.emm)

contrast(ht35.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# LEAVES at 35 DAYS
# ORIGINAL
model <- lmer(data=comb.plant.clean.new, leaves.35.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
           + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.35.new) + (1 | repro.35.new))

# FINAL
# Nonsignificant: all significant.
lf35.model.reduced <- lmer(data=comb.plant.clean.new, leaves.35.new~genus + lifespan + genus*lifespan + genus*lifespan/species.short
     + (1 | species.short:accession) + (1 | accession:replicate) + (1 | vigor.35.new) + (1 | repro.35.new))

rand(model)
r <- rand(lf35.model.reduced)
anova(lf35.model.reduced)
lf35.mod.emm <- emmeans(object=lf35.model.reduced, specs=c("lifespan","genus"))
lf35.df <- as.data.frame(lf35.mod.emm)

contrast(lf35.mod.emm, method=list("LathA-LathP" = LathA-LathP,
                              "PhasA-PhasP" = PhasA-PhasP,
                              "ViciA-ViciP" = ViciA-ViciP),
         adjust="bonferroni")

# Subsets of data
lath.plant <- comb.plant.clean.new[comb.plant.clean.new$genus=="Lathyrus",]
phas.plant <- comb.plant.clean.new[comb.plant.clean.new$genus=="Phaseolus",]
vic.plant <- comb.plant.clean.new[comb.plant.clean.new$genus=="Vicia",]

comb.plant.clean.fert <- comb.plant.clean.new[!comb.plant.clean.new$days.since.fert.21 == 0,] 

ggplot(data=phas.plant, aes(x=growth.ratio, y=height.perday.21.35, color = lifespan))+
  geom_point(alpha=0.5)

ggplot(data=phas.plant, aes(x=days.since.fert.21, y=height.21.new, color = lifespan))+
  geom_point(alpha=0.5)

ggplot(data=phas.plant, aes(x=days.since.fert.35, y=height.35.new, color = lifespan))+
  geom_point(alpha=0.5)

summary(lm(data=phas.plant, growth.ratio~height.perday.21.35))

summary(lm(data=phas.plant, days.since.fert.21~height.21.new))

summary(lm(data=phas.plant, days.since.fert.35~height.35.new))

ggplot(data=comb.plant.clean.new, aes(x=pesticide, y=height.perday.21.35, color = genus))+
  geom_boxplot()

ggplot(data=comb.plant.clean.new, aes(x=pesticide.21, y=nodes.21.new, color = genus))+
  geom_boxplot()

ggplot(data=comb.plant.clean.new, aes(x=pesticide.35, y=nodes.35.new, color = genus))+
  geom_boxplot()

#######################################################


# TEST PLOTS & OUTLIERS ####
lath.merge <- merge.acc.all[merge.acc.all$genus=="Lathyrus",]
phas.merge <- merge.acc.all[merge.acc.all$genus=="Phaseolus",]
vic.merge <- merge.acc.all[merge.acc.all$genus=="Vicia",]
ann.merge <- merge.acc.all[merge.acc.all$lifespan=="annual",]
per.merge <- merge.acc.all[merge.acc.all$lifespan=="perennial",]

ggplot(data=lath.merge, aes(x=seed.perim.mn, y=height.perday.mn, 
   shape=genus, color=lifespan))+
  geom_point(size=3)

ggplot(data=phas.merge, aes(x=seed.area.mn, y=height.21.mn, 
                  shape=genus, color=lifespan, label=accession))+
  geom_point(size=3)+
  geom_text(color="black", size=3)

ggplot(data=per.merge, aes(x=lifespan, y=height.21.new, color=accession))+
  geom_point(size=3, position=position_jitterdodge(jitter.width=0.5, dodge.width = 0))

# examine outliers function
outlier.look <- function(d, trait) {
  plot <- ggplot(data=d, aes(x=trait, y=trait, 
                                label=accession))+
    geom_point(size=3, color="cyan")+
    geom_text(color="black", size=3)
  
  Q <- quantile(trait, probs=c(.25, .75), na.rm = T)
  iqr <- IQR(trait,  na.rm = T)
  upper <- Q[2]+1.5*iqr # Upper Range  
  lower <- Q[1]-1.5*iqr # Lower Range
  my_list <- list(plot, upper, lower)
  return(my_list)
}

outlier.look(d = vic.merge, trait = vic.merge$FINAL.germ.prop)

# examining outliers
ggplot(data=lath.germ, aes(x=ppd50, y=ppd50, 
                            color=accession))+
  geom_point(size=3, alpha=0.3, position=position_jitterdodge(jitter.width=0.5))

# test averages and SD
# Seed
lath.seed <- comb.seed.clean %>% filter(genus=="Lathyrus")
phas.seed <- comb.seed.clean %>% filter(genus=="Phaseolus")
vic.seed <- comb.seed.clean %>% filter(genus=="Vicia")
ann.lath.seed <- comb.seed.clean %>% filter(genus=="Lathyrus" & lifespan=="annual")
per.lath.seed <- comb.seed.clean %>% filter(genus=="Lathyrus" & lifespan=="perennial")
ann.phas.seed <- comb.seed.clean %>% filter(genus=="Phaseolus" & lifespan=="annual")
per.phas.seed <- comb.seed.clean %>% filter(genus=="Phaseolus" & lifespan=="perennial")
ann.vic.seed <- comb.seed.clean %>% filter(genus=="Vicia" & lifespan=="annual")
per.vic.seed <- comb.seed.clean %>% filter(genus=="Vicia" & lifespan=="perennial")

# Germination T50
lath.germ <- comb.germ.clean %>% filter(genus=="Lathyrus")
phas.germ <- comb.germ.clean %>% filter(genus=="Phaseolus")
vic.germ <- comb.germ.clean %>% filter(genus=="Vicia")
ann.lath.germ <- comb.germ.clean %>% filter(genus=="Lathyrus" & lifespan=="annual")
per.lath.germ <- comb.germ.clean %>% filter(genus=="Lathyrus" & lifespan=="perennial")
ann.phas.germ <- comb.germ.clean %>% filter(genus=="Phaseolus" & lifespan=="annual")
per.phas.germ <- comb.germ.clean %>% filter(genus=="Phaseolus" & lifespan=="perennial")
ann.vic.germ <- comb.germ.clean %>% filter(genus=="Vicia" & lifespan=="annual")
per.vic.germ <- comb.germ.clean %>% filter(genus=="Vicia" & lifespan=="perennial")

# Vegetative
lath.veg <- comb.plant.clean.new %>% filter(genus=="Lathyrus")
phas.veg <- comb.plant.clean.new %>% filter(genus=="Phaseolus")
vic.veg <- comb.plant.clean.new %>% filter(genus=="Vicia")
ann.lath.veg <- comb.plant.clean.new %>% filter(genus=="Lathyrus" & lifespan=="annual")
per.lath.veg <- comb.plant.clean.new %>% filter(genus=="Lathyrus" & lifespan=="perennial")
ann.phas.veg <- comb.plant.clean.new %>% filter(genus=="Phaseolus" & lifespan=="annual")
per.phas.veg <- comb.plant.clean.new %>% filter(genus=="Phaseolus" & lifespan=="perennial")
ann.vic.veg <- comb.plant.clean.new %>% filter(genus=="Vicia" & lifespan=="annual")
per.vic.veg <- comb.plant.clean.new %>% filter(genus=="Vicia" & lifespan=="perennial")

# Accession level
ann.lath.acc <- merge.acc.all %>% filter(genus=="Lathyrus" & lifespan=="annual")
per.lath.acc <- merge.acc.all %>% filter(genus=="Lathyrus" & lifespan=="perennial")
ann.phas.acc <- merge.acc.all %>% filter(genus=="Phaseolus" & lifespan=="annual")
per.phas.acc <- merge.acc.all %>% filter(genus=="Phaseolus" & lifespan=="perennial")
ann.vic.acc <- merge.acc.all %>% filter(genus=="Vicia" & lifespan=="annual")
per.vic.acc <- merge.acc.all %>% filter(genus=="Vicia" & lifespan=="perennial")

mean(per.vic.acc$leaves.rgr.mn, na.rm=T)
sd(per.vic.acc$leaves.rgr.mn, na.rm=T)
# ann.lath
# per.lath
# ann.phas
# per.phas
# ann.vic
# per.vic

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ OLD ... ####

# data truncating outliers & testing correlations
bad.seed.corr <- c("PI 535340","PI 535348","PI 535202","PI 661810","PI 661811", "PI 640969")
ann.bad.corr <- ann.merge[!ann.merge$accession %in% bad.seed.corr,]
all.bad.corr <- merge.acc.all[!merge.acc.all$accession %in% bad.seed.corr,]
phas.bad.corr <- phas.merge[!phas.merge$accession %in% bad.seed.corr,]

bad.seed.corr1 <- c("PI 642133","W6 39937","PI 535372")
per.bad.corr <- per.merge[!per.merge$accession %in% bad.seed.corr1,]

per.bad.germ <- c("PI 494749", "PI 372552")
per.germ.corr <- per.merge[!per.merge$accession %in% per.bad.germ,]

ann.bad.germ <- c("PI 494702", "PI 219631", "W6 17772") # had the longest T50s for annuals
ann.germ.corr <- ann.merge[!ann.merge$accession %in% ann.bad.germ,]

lath.bad.veg <- c("PI 358829", "W6 18210", "W6 20260", "PI 358865")
lath.bad.corr <- lath.merge[!lath.merge$accession %in% lath.bad.veg,]

da <- phas.bad.corr %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
da <- scale(da, scale=TRUE, center=TRUE)
xa <- cor(da, use="complete.obs") # allows NA's in matrix, & removes them
resa <- cor.mtest(da, conf.level = .95)
colnames(resa$p) <- colnames(xa)
rownames(resa$p) <- rownames(xa)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xa) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xa) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xa, order="original", method="color", col=col,
         p.mat = resa$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, title="Test!", mar=c(0,0,1,0))

# make smaller label names for the network nodes
colnames(xa) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf No.\nDAP-21", "Leaf No.\nDAP-35", "Leaf No.\nAGR", "Leaf No.\nRGR")
screen_mata <- resa$p
screen_mata[screen_mata > 0.05] <- 0
xa[screen_mata == 0] <- 0
g_rawa <- graph_from_adjacency_matrix(xa, mode='upper', weighted=T, diag=F)
ga <- graph_from_adjacency_matrix(abs(xa)^2, mode='upper', weighted=T, diag=F)
E(ga)$color[E(g_rawa)$weight > 0] <- 'blue3'
E(ga)$color[E(g_rawa)$weight < 0] <- 'brown3'
par(mar=c(0,0,1,0))
plot <- plot(ga, layout=legume.mat, edge.width=abs(E(g)$weight)*7, vertex.color=col, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.color="black")

# OLD CORRELATIONS WITHOUT P VALUE ADJUSTMENTS) ####
# FULL DATASET NETWORK & CORRELATIONS ###
# Have to start by making a correlation matrix
d <- merge.acc.all %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
d <- scale(d, scale=TRUE, center=TRUE)
x <- cor(d, use="complete.obs") # allows NA's in matrix, & removes them
res1 <- cor.mtest(d, conf.level = .95) # Pearson is the default
colnames(res1$p) <- colnames(x)
rownames(res1$p) <- rownames(x)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(x) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(x) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(x, order="original", method="color", col=col,
         p.mat = res1$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Full dataset", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(x) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")

# NETWORK!!!
screen_mat <- res1$p
screen_mat[screen_mat > 0.05] <- 0
x[screen_mat == 0] <- 0
g_raw <- graph_from_adjacency_matrix(x, mode='upper', weighted=T, diag=F)
g <- graph_from_adjacency_matrix(abs(x)^2, mode='upper', weighted=T, diag=F) # use r2 for the correlation values, since iGraph can't handle negative corr values.
E(g)$color[E(g_raw)$weight > 0] <- 'blue3'
E(g)$color[E(g_raw)$weight < 0] <- 'brown3'
plot(g, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white') # auto, line thickness increased 5x
plot(g, layout=layout.circle, edge.width=abs(E(g)$weight)*5, vertex.color='white') # circle
plot(g, layout=layout.sphere, edge.width=abs(E(g)$weight)*5, vertex.color='white') # sphere
plot(g, layout=layout.grid, edge.width=abs(E(g)$weight)*5, vertex.color='white') # grid

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg <- degree(g, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col <- colrange(max(deg)+1)
col <- col[deg+1]

# get coordinates of tkplot to replot with the norm plot function
tk <- tkplot(gv, layout=layout_with_fr, edge.width=abs(E(g)$weight)*7, vertex.color=col.vicia, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.color="black", main="Full dataset")

# first move the nodes where you want them and then get the coordinates
# Keep the tk plot open!
layout.legume <- tkplot.getcoords(tk)

# matrix coords for future reference
legume.mat <- cbind(c(130,187,41,221,71,144,68,312,397,272,287,304,272,416,399,388,418), 
                    c(199,3,127,124,2,370,314,384,356,282,183,90,0,281,184,88,0))
legume.mat==layout.legume # check that it's the same

# Set up grid of plots
par(mar=c(0,0,1,0), mfrow=c(2,3))

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(g, layout=legume.mat, edge.width=abs(E(g)$weight)*7, vertex.color=col, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("A  Full dataset", adj=0.1, line=-0.5)

# trying changing color based on betweenness
deg <- betweenness(g, directed=T, weights=NA)
colrange <- colorRampPalette(c("yellow", "red"))
col <- colrange(max(deg)+1)
col <- col[deg+1]

plot(g, layout=layout_with_fr, edge.width=abs(E(g)$weight)*5, vertex.color=col, 
     vertex.size=20, vertex.label.cex=0.45, vertex.label.font=2, vertex.label.color="black", main="Full dataset")

tkplot(g, layout=layout_with_fr, edge.width=abs(E(g)$weight)*5, vertex.color=col, 
       vertex.size=18, vertex.label.cex=0.45, vertex.label.font=2, vertex.label.color="black", main="Full dataset")

# ANNUAL NETWORK & CORRELATIONS ###
ann.tot <- merge.acc.all %>% filter(lifespan=="annual")
da <- ann.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
da <- scale(da, scale=TRUE, center=TRUE)
xa <- cor(da, use="complete.obs") # allows NA's in matrix, & removes them
resa <- cor.mtest(da, conf.level = .95)
colnames(resa$p) <- colnames(xa)
rownames(resa$p) <- rownames(xa)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xa) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xa) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xa, order="original", method="color", col=col,
         p.mat = resa$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Annual", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xa) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_mata <- resa$p
screen_mata[screen_mata > 0.05] <- 0
xa[screen_mata == 0] <- 0
g_rawa <- graph_from_adjacency_matrix(xa, mode='upper', weighted=T, diag=F)
ga <- graph_from_adjacency_matrix(abs(xa)^2, mode='upper', weighted=T, diag=F)
E(ga)$color[E(g_rawa)$weight > 0] <- 'blue3'
E(ga)$color[E(g_rawa)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.ann <- degree(ga, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.ann <- colrange(max(deg.ann)+1)
col.ann <- col.ann[deg.ann+1]

# Keep the tk plot open!
# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(ga, layout=legume.mat, edge.width=abs(E(ga)$weight)*7, vertex.color=col.ann, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("B  Annual", adj=0.1, line=-0.5)

# PERENNIAL NETWORK & CORRELATIONS ###
per.tot <- merge.acc.all %>% filter(lifespan=="perennial")
dper <- per.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dper <- scale(dper, scale=TRUE, center=TRUE)
xper <- cor(dper, use="complete.obs") # allows NA's in matrix, & removes them
resper <- cor.mtest(dper, conf.level = .95)
colnames(resper$p) <- colnames(xper)
rownames(resper$p) <- rownames(xper)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xper) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xper) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xper, order="original", method="color", col=col,
         p.mat = resper$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Perennial", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xper) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matper <- resper$p
screen_matper[screen_matper > 0.05] <- 0
xper[screen_matper == 0] <- 0
g_rawper <- graph_from_adjacency_matrix(xper, mode='upper', weighted=T, diag=F)
gper <- graph_from_adjacency_matrix(abs(xper)^2, mode='upper', weighted=T, diag=F)
E(gper)$color[E(g_rawper)$weight > 0] <- 'blue3'
E(gper)$color[E(g_rawper)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.per <- degree(gper, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.per <- colrange(max(deg.per)+1)
col.per <- col.per[deg.per+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gper, layout=legume.mat, edge.width=abs(E(gper)$weight)*7, vertex.color=col.per, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("C  Perennial", adj=0.1, line=-0.5)

# LATHYRUS NETWORK & CORRELATIONS ###
lath.tot <- merge.acc.all %>% filter(genus=="Lathyrus")
dl <- lath.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dl <- scale(dl, scale=TRUE, center=TRUE)
xl <- cor(dl, use="complete.obs") # allows NA's in matrix, & removes them
resl <- cor.mtest(dl, conf.level = .95)
colnames(resl$p) <- colnames(xl)
rownames(resl$p) <- rownames(xl)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xl) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xl) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xl, order="original", method="color", col=col,
         p.mat = resl$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Lathyrus", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xl) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matl <- resl$p
screen_matl[screen_matl > 0.05] <- 0
xl[screen_matl == 0] <- 0
g_rawl <- graph_from_adjacency_matrix(xl, mode='upper', weighted=T, diag=F)
gl <- graph_from_adjacency_matrix(abs(xl)^2, mode='upper', weighted=T, diag=F)
E(gl)$color[E(g_rawl)$weight > 0] <- 'blue3'
E(gl)$color[E(g_rawl)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.lath <- degree(gl, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.lath <- colrange(max(deg.lath)+1)
col.lath <- col.lath[deg.lath+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gl, layout=legume.mat, edge.width=abs(E(gl)$weight)*7, vertex.color=col.lath, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("D  Lathyrus", adj=0.1, line=-0.5)


# PHASEOLUS NETWORK & CORRELATIONS ###
phas.tot <- merge.acc.all %>% filter(genus=="Phaseolus")
dp <- phas.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dp <- scale(dp, scale=TRUE, center=TRUE)
xp <- cor(dp, use="complete.obs") # allows NA's in matrix, & removes them
resp <- cor.mtest(dp, conf.level = .95)
colnames(resp$p) <- colnames(xp)
rownames(resp$p) <- rownames(xp)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xp) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xp) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xp, order="original", method="color", col=col,
         p.mat = resp$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Phaseolus", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xp) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matp <- resp$p
screen_matp[screen_matp > 0.05] <- 0
xp[screen_matp == 0] <- 0
g_rawp <- graph_from_adjacency_matrix(xp, mode='upper', weighted=T, diag=F)
gp <- graph_from_adjacency_matrix(abs(xp)^2, mode='upper', weighted=T, diag=F)
E(gp)$color[E(g_rawp)$weight > 0] <- 'blue3'
E(gp)$color[E(g_rawp)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.phas <- degree(gp, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.phas <- colrange(max(deg.phas)+1)
col.phas <- col.phas[deg.phas+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gp, layout=legume.mat, edge.width=abs(E(gp)$weight)*7, vertex.color=col.phas, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("E  Phaseolus", adj=0.1, line=-0.5)


# VICIA NETWORK & CORRELATIONS ###
vicia.tot <- merge.acc.all %>% filter(genus=="Vicia")
dv <- vicia.tot %>% dplyr::select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, height.rgr.mn, leaves.21.mn, leaves.35.mn, leaves.perday.mn, leaves.rgr.mn)
dv <- scale(dv, scale=TRUE, center=TRUE)
xv <- cor(dv, use="complete.obs") # allows NA's in matrix, & removes them
resv <- cor.mtest(dv, conf.level = .95)
colnames(resv$p) <- colnames(xv)
rownames(resv$p) <- rownames(xv)
# full corrplot for supplement (have to do before network screen_mat below, or some coefficients = 0)
col<- colorRampPalette(c("firebrick", "white", "dodgerblue4"))(20)
colnames(xv) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
rownames(xv) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity","Seed roundness", "Germination T50","Germination proportion", "Height DAP-21", "Height DAP-35", "Height AGR", "Height RGR" , "Leaf number DAP-21", "Leaf number DAP-35", "Leaf number AGR", "Leaf number RGR")
corrplot(xv, order="original", method="color", col=col,
         p.mat = resv$p, sig.level = 0.05, insig = "blank", addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex=0.5, number.cex=0.4, cl.cex = 0.6, mar=c(0,0,0,0))
mtext("Vicia", at=0, line=1.5, cex=1.25)

# make smaller label names for the network nodes
colnames(xv) <- c("Seed\nmass", "Seed\nlength", "Seed\nwidth", "Seed\nperim.", "Seed\narea", "Seed\ncirc.","Seed\nround.", "Germ.\nT50","Germ.\nprop.", "Height\nDAP-21", "Height\nDAP-35", "Height\nAGR", "Height\nRGR", "Leaf N\nDAP-21", "Leaf N\nDAP-35", "Leaf N\nAGR", "Leaf N\nRGR")
screen_matv <- resv$p
screen_matv[screen_matv > 0.05] <- 0
xv[screen_matv == 0] <- 0
g_rawv <- graph_from_adjacency_matrix(xv, mode='upper', weighted=T, diag=F)
gv <- graph_from_adjacency_matrix(abs(xv)^2, mode='upper', weighted=T, diag=F)
E(gv)$color[E(g_rawv)$weight > 0] <- 'blue3'
E(gv)$color[E(g_rawv)$weight < 0] <- 'brown3'

# Have to set margins so you can get a title and not a ton of white space
par(mar=c(0,0,1,0))

# Change color based on increasing degree
deg.vicia <- degree(gv, mode="all")
colrange <- colorRampPalette(c("yellow", "red"))
col.vicia <- colrange(max(deg.vicia)+1)
col.vicia <- col.vicia[deg.vicia+1]

# still have to put in the usual code, but with the layout = tkplot coordinates
plot <- plot(gv, layout=legume.mat, edge.width=abs(E(gv)$weight)*7, vertex.color=col.vicia, 
             vertex.size=26, vertex.label.cex=0.6, vertex.label.font=2, vertex.label.family="Helvetica", vertex.label.color="black")
title("F  Vicia", adj=0.1, line=-0.5)

# testing specific relationships
plot <- ggplot(data=lath.tot, aes(x=ppd50.mn, y=nodes.21.mn, color=lifespan))+
  geom_point()+labs(title="annual")
plot

plot <- ggplot(data=merge.acc.all, aes(x=ppd50.mn, y=nodes.35.mn, color=genus))+
  geom_point()+labs(title="all")
plot

merge.acc.all
# Two methods of making a grid with the network plots:

# make a grid
map_base_to_grid <- function(fun) {
  gridGraphics::grid.echo(fun)
  grid::grid.grab()
}

grid1 <- map_base_to_grid(function() plot(g, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid2 <- map_base_to_grid(function() plot(gl, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid3 <- map_base_to_grid(function() plot(gp, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid4 <- map_base_to_grid(function() plot(gv, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid5 <- map_base_to_grid(function() plot(ga, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))
grid6 <- map_base_to_grid(function() plot(gper, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                                          vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0))

ggpubr::ggarrange(grid1, grid2, grid3, grid4, grid5, grid6, nrow=3, ncol=2, widths=c(2,2), heights=c(2,2,2))

plotit <- function(input){plot(input, layout=layout.auto, edge.width=abs(E(g)$weight)*5, vertex.color='white', 
                               vertex.size=22, vertex.label.cex=0.45, vertex.label.font=2, margin = c(0,0,0,0), asp=0)}
gridfun.all <- function()plotit(input=g)
gridfun.lath <- function()plotit(input=gl)
gridfun.phas <- function()plotit(input=gp)
gridfun.vic <- function()plotit(input=gv)
gridfun.ann <- function()plotit(input=ga)
gridfun.per <- function()plotit(input=gper)

ggpubr::ggarrange(gridfun.all, gridfun.lath, gridfun.phas, gridfun.vic, gridfun.ann, gridfun.per, nrow=3, ncol=2, widths=c(2,2), heights=c(2,2,2))

# Network For PCs
# merge PCA outputs
merge.pca <- Reduce(function(x,y) merge(x,y, all=T, by=c("accession")), list(seed.acc.pca, veg.acc.pca))
d.pca <- merge.pca %>% dplyr::select(sPC1:vPC4)
d.pca <- scale(d.pca, scale=TRUE, center=TRUE)
x.pca <- cor(d.pca, use="complete.obs") # allows NA's in matrix, & removes them
res.pca <- cor.mtest(d.pca, conf.level = .95)
colnames(res.pca$p) <- colnames(x.pca)
rownames(res.pca$p) <- rownames(x.pca)
screen_mat.pca <- res.pca$p
screen_mat.pca[screen_mat.pca > 0.05] <- 0
x.pca[screen_mat.pca == 0] <- 0
g_raw.pca <- graph_from_adjacency_matrix(x.pca, mode='upper', weighted=T, diag=F)
g.pca <- graph_from_adjacency_matrix(abs(x.pca)^2, mode='upper', weighted=T, diag=F)
E(g.pca)$color[E(g_raw.pca)$weight > 0] <- 'skyblue1'
E(g.pca)$color[E(g_raw.pca)$weight < 0] <- 'salmon1'
plot(g.pca, layout=layout.auto, edge.width=abs(E(g.pca)$weight)*30, 
     vertex.color='white', vertex.size=23, vertex.label.cex=0.7)
plot(g.pca, layout=layout.circle, edge.width=abs(E(g.pca)$weight)*5, vertex.color='white') # circle
plot(g.pca, layout=layout.sphere, edge.width=abs(E(g.pca)$weight)*5, vertex.color='white') # sphere
plot(g.pca, layout=layout.grid, edge.width=abs(E(g.pca)$weight)*5, vertex.color='white') # grid

# SUMMARIES FOR FIGURES (by genus & lifespan) # CURRENTLY OBSOLETE ####
# USE THE FILTERED DATA!! (above)

# FOR FIGURES ImageJ seed data
seed.fig <- comb.seed.clean %>%
  group_by(lifespan, genus) %>%
  summarise(
    seed.area.mn = mean(seed.area, na.rm=T),
    seed.area.sd = sd(seed.area, na.rm=T),
    seed.length.mn = mean(seed.length, na.rm=T),
    seed.length.sd = sd(seed.length, na.rm=T),
    seed.width.mn = mean(seed.width, na.rm=T),
    seed.width.sd = sd(seed.width, na.rm=T),
    seed.circ.mn = mean(seed.circ, na.rm=T),
    seed.circ.sd = sd(seed.circ, na.rm=T),
    seed.perim.mn = mean(seed.perim, na.rm=T),
    seed.perim.sd = sd(seed.perim, na.rm=T)
  ) %>% ungroup ()
View(seed.fig)

# FOR FIGURES Seed mass data
seedwt.fig <- comb.acc.clean %>%
  group_by(lifespan, genus) %>%
  summarise(
    seed.wt.mn = mean(seed.wt.avg, na.rm=T),
    seed.wt.sd = sd(seed.wt.avg, na.rm=T),
  ) %>% ungroup ()
View(seedwt.fig)

# FOR FIGURES Germination proportion data
germ.fig <- acc.germ.clean %>%
  group_by(lifespan, genus) %>%
  summarise(    
    germ.prop.mn = mean(FINAL.germ.prop, na.rm=T),
    germ.prop.sd = sd(FINAL.germ.prop, na.rm=T),
    viab.prop.mn = mean(FINAL.viable.prop, na.rm=T),
    viab.prop.sd = sd(FINAL.viable.prop, na.rm=T),
  ) %>% ungroup ()
View(germ.fig)

# FOR FIGURES Germination rate data
germrate.fig <- comb.germ.clean %>%
  group_by(lifespan, genus) %>%
  summarise(    
    ppd50.fig.mn = mean(ppd50, na.rm=T),
    ppd50.fig.sd = sd(ppd50, na.rm=T),
  ) %>% ungroup ()
View(germrate.fig)

# FOR FIGURES Growth data
growth.fig <- comb.plant.clean.new %>%
  group_by(lifespan, genus) %>%
  summarise(    
    height.perday.mn = mean(height.perday.21.35, na.rm=T),
    height.perday.sd = sd(height.perday.21.35, na.rm=T),
    nodes.perday.mn = mean(nodes.perday.21.35, na.rm=T),
    nodes.perday.sd = sd(nodes.perday.21.35, na.rm=T),
    height.21.mn = mean(height.21.new, na.rm=T),
    height.21.sd = sd(height.21.new, na.rm=T),
    nodes.21.mn = mean(nodes.21.new, na.rm=T),
    nodes.21.sd = sd(nodes.21.new, na.rm=T),
    height.35.mn = mean(height.35.new, na.rm=T),
    height.35.sd = sd(height.35.new, na.rm=T),
    nodes.35.mn = mean(nodes.35.new, na.rm=T),
    nodes.35.sd = sd(nodes.35.new, na.rm=T),
  ) %>% ungroup ()
View(growth.fig)

# FOR FIGURES leaf data
leaf.fig <- comb.leaf.clean %>%
  group_by(lifespan, genus) %>%
  summarise(    
    leaf.fresh.mn = mean(leaf.fresh, na.rm=T),
    leaf.fresh.sd = sd(leaf.fresh, na.rm=T),
    leaf.dry.mn = mean(leaf.dry, na.rm=T),
    leaf.dry.sd = sd(leaf.dry, na.rm=T),
    leaf.area.mn = mean(leaf.area, na.rm=T),
    leaf.area.sd = sd(leaf.area, na.rm=T),
    ldrm.mn = mean(ldmc, na.rm=T),
    ldrm.sd = sd(ldmc, na.rm=T),
    sla.mn = mean(sla, na.rm=T),
    sla.sd = sd(sla, na.rm=T),
  ) %>% ungroup ()
View(leaf.fig)

# Interaction plots ####
# 'reaction norm plot' with points for each lifespan type, and color and shape by genus.
# You will need the summary datasets above for these.
# to remove the legend and axes tick labels:
# add axis.text.x = element_blank(), legend.position = "none"
# to theme()

# Seed mass
p.sw <- ggplot(data=sdwt.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
sw.plot <- p.sw + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Averaged single seed mass (mg)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 60, 20), limits = c(0,60), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sw.plot

# Seed length
p.sl <- ggplot(data=sdlg.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
sl.plot <- p.sl + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Seed length (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 7, 1), limits = c(0,7), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sl.plot

# Seed width
p.swd <- ggplot(data=sdwd.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
swd.plot <- p.swd + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Seed width (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 5, 1), limits = c(0,5), expand = c(0,0))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
swd.plot

# Seed perimeter
p.sp <- ggplot(data=sdpm.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
sp.plot <- p.sp + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y=bquote("Seed perimeter (mm)"))+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0,20), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sp.plot

# Seed area
p.sa <- ggplot(data=sdar.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
sa.plot <- p.sa + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y=bquote("Seed area "~(mm^2)))+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0,25), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sa.plot

# Seed circularity
p.sc <- ggplot(data=sdcr.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
sc.plot <- p.sc + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Seed circularity")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0.82, 0.92, 0.02), limits = c(0.82,0.92), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sc.plot

p.sr <- ggplot(data=sdrd.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
sr.plot <- p.sr + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Seed roundness")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0.75, 0.95, 0.05), limits = c(0.75, 0.95), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sr.plot

# Germination rate (T50)
# Lathyrus T50
genus <- c("Lathyrus", "Lathyrus") # need to add a genus column to get the line between the points
lath.gmrt.df %>% add_column(genus)

lath.gr.plot <- ggplot(data=lath.gmrt.df, aes(x=lifespan, y=emmean, group=genus))+
  geom_line(size=.75, color="midnightblue")+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.1, color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, shape=21, color="black", fill="midnightblue", stat="identity")+
  labs(x="Lifespan",y="T50 (days)", title="Lathyrus")+
  scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0,20), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
lath.gr.plot

# Phaseolus T50
genus <- c("Phaseolus", "Phaseolus") # need to add a genus column to get the line between the points
phas.gmrt.df %>% add_column(genus)

phas.gr.plot <- ggplot(data=phas.gmrt.df, aes(x=lifespan, y=emmean, group=genus))+
  geom_line(size=.75, color="grey50")+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.1, color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, shape=22, color="black", fill="grey50", stat="identity")+
  labs(x="Lifespan",y="T50 (days)", title="Phaseolus")+
  scale_y_continuous(breaks = seq(0, 3, 1), limits = c(0,3), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
phas.gr.plot

# Vicia T50
genus <- c("Vicia", "Vicia") # need to add a genus column to get the line between the points
vic.gmrt.df %>% add_column(genus)

vic.gr.plot <- ggplot(data=vic.gmrt.df, aes(x=lifespan, y=emmean, group=genus))+
  geom_line(size=.75, color="dodgerblue")+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.1, color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, shape=24, color="black", fill="dodgerblue", stat="identity")+
  labs(x="Lifespan",y="T50 (days)", title="Vicia")+
  scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0,20), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
vic.gr.plot

# Germination proportion

# Lathyrus germ prop
genus <- c("Lathyrus", "Lathyrus") # need to add a genus column to get the line between the points
lath.gmprop.df %>% add_column(genus)

lath.gp.plot <- ggplot(data=lath.gmprop.df, aes(x=lifespan, y=emmean, group=genus))+
  geom_line(size=.75, color="midnightblue")+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.1, color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, shape=21, color="black", fill="midnightblue", stat="identity")+
  labs(x="Lifespan",y="Germination proportion")+
  scale_y_continuous(breaks = seq(0, 1.1, 0.25), limits = c(0,1.1), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
lath.gp.plot

# Phaseolus germ prop
genus <- c("Phaseolus", "Phaseolus") # need to add a genus column to get the line between the points
phas.gmprop.df %>% add_column(genus)

phas.gp.plot <- ggplot(data=phas.gmprop.df, aes(x=lifespan, y=emmean, group=genus))+
  geom_line(size=.75, color="grey50")+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.1, color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, shape=22, color="black", fill="grey50", stat="identity")+
  labs(x="Lifespan",y="Germination proportion")+
  scale_y_continuous(breaks = seq(0, 1.1, 0.25), limits = c(0,1.1), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
phas.gp.plot

# Vicia germ prop
genus <- c("Vicia", "Vicia") # need to add a genus column to get the line between the points
vic.gmprop.df %>% add_column(genus)

vic.gp.plot <- ggplot(data=vic.gmprop.df, aes(x=lifespan, y=emmean, group=genus))+
  geom_line(size=.75, color="dodgerblue")+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.1, color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, shape=24, color="black", fill="dodgerblue", stat="identity")+
  labs(x="Lifespan",y="Germination proportion")+
  scale_y_continuous(breaks = seq(0, 1.1, 0.25), limits = c(0,1.1), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
vic.gp.plot

# Viable proportion
p.vp <- ggplot(data=germ.fig, aes(x=lifespan, y=viab.prop.mn, group=genus, shape=genus, color=genus, fill=genus))
vp.plot <- p.vp + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=viab.prop.mn-viab.prop.sd, ymax=viab.prop.mn+viab.prop.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Viable proportion")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
vp.plot

# Height growth
p.ht <- ggplot(data=htagr.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
ht.plot <- p.ht + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Height raw growth rate (mm/day)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 20, 5), limits = c(0,20), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
ht.plot

# Height RGR
p.htrgr <- ggplot(data=htrgr.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
htrgr.plot <- p.htrgr + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Height relative growth rate")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 0.07, 0.01), limits = c(0,0.07), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
htrgr.plot

# Leaf growth
p.lf <- ggplot(data=lfagr.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
lf.plot <- p.lf + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Leaf raw growth rate (leaves/day)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 0.4, 0.1), limits = c(0,0.4), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
lf.plot

# Height 21 days
p.ht.21 <- ggplot(data=ht21.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
ht.21.plot <- p.ht.21 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Height DAP-21 (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 300, 100), limits = c(0,300), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
ht.21.plot

# Leaves 21 days
p.lf.21 <- ggplot(data=lf21.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
lf.21.plot <- p.lf.21 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Leaf number DAP-21")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0,10), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
lf.21.plot

# Height 35 days
p.ht.35 <- ggplot(data=ht35.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
ht.35.plot <- p.ht.35 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Height DAP-35 (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 600, 100), limits = c(0,600), expand = c(0,0))+
  expand_limits(y = 0)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
ht.35.plot

# Leaves 35 days
p.lf.35 <- ggplot(data=lf35.df, aes(x=lifespan, y=emmean, group=genus, shape=genus, color=genus, fill=genus))
lf.35.plot <- p.lf.35 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),
                width=0.2, position = position_dodge(.4), color="black", size=0.35)+
  geom_point(size=2.5, stroke=0.5, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Leaf number DAP-35")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_y_continuous(breaks = seq(0, 14, 2), limits = c(0,14), expand = c(0,0))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
lf.35.plot

# GRID of main interaction plots ####
# ggarrange has to be from the ggpubr package (not egg) to get the shared legend

# Seed plots
# USE THIS ONE ->
ggpubr::ggarrange(sw.plot, sl.plot, swd.plot, sp.plot, sa.plot, sc.plot, sr.plot, nrow=2, ncol=4, widths=c(2,2,2,2), heights=c(2,2),
                  labels=c("A","B","C","D","E","F","G"), common.legend=T, legend="bottom")

# GRID of VEGETATIVE PLOTS
ggpubr::ggarrange(ht.21.plot, ht.35.plot, ht.plot, htrgr.plot, lf.21.plot, lf.35.plot, lf.plot, nrow=2, ncol=4, widths=c(2,2,2,2), heights=c(2,2),
                  labels=c("A","B","C","D","E","F","G"), common.legend=T, legend="bottom")

# GRID of germination plots
ggpubr::ggarrange(lath.gr.plot, phas.gr.plot, vic.gr.plot, lath.gp.plot, phas.gp.plot, vic.gp.plot, nrow=2, ncol=3, widths=c(2,2,2), heights=c(2,2),
                  labels=c("A", "B", "C", "D", "E","F"), common.legend=T, legend="bottom")

# LEAF plots ####
cut.species <- c("L. japonicus", "P. maculatus")

# remove species with only 1-2 points:
comb.leaf.clean <- comb.leaf[!comb.leaf$species.short %in% cut.species,]
comb.leaf.lath <- comb.leaf.clean[comb.leaf.clean$genus=="Lathyrus",]
comb.leaf.phas <- comb.leaf.clean[comb.leaf.clean$genus=="Phaseolus",]
comb.leaf.vic <- comb.leaf.clean[comb.leaf.clean$genus=="Vicia",]

# For dotplot stacked in center:
# geom_dotplot(binaxis='y', stackdir='center', dotsize=1, aes(fill=accession), show.legend=F)+
# to reorder by median: aes(x=reorder(specific.epithet, sla, median, na.rm=T)
# y=bquote("Seed area "~(mm^2))

# SLA
sla.lath.p <- ggplot(data=comb.leaf.lath, aes(x=specific.epithet, y=sla))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("SLA "~(mm^2 ~mg^-1)), title="Lathyrus")+
  theme_classic()
sla.lath.p

sla.phas.p <- ggplot(data=comb.leaf.phas, aes(x=specific.epithet, y=sla))+ 
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("SLA "~(mm^2 ~mg^-1)), title="Phaseolus")+
  theme_classic()
sla.phas.p

sla.vic.p <- ggplot(data=comb.leaf.vic, aes(x=specific.epithet, y=sla))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("SLA "~(mm^2 ~mg^-1)), title="Vicia")+
  theme_classic()
sla.vic.p

# LDMC
ldmc.lath.p <- ggplot(data=comb.leaf.lath, aes(x=specific.epithet, y=ldmc))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("LDMC "~(mg ~g^-1)))+
  theme_classic()
ldmc.lath.p

ldmc.phas.p <- ggplot(data=comb.leaf.phas, aes(x=specific.epithet, y=ldmc))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("LDMC "~(mg ~g^-1)))+
  theme_classic()
ldmc.phas.p

ldmc.vic.p <- ggplot(data=comb.leaf.vic, aes(x=specific.epithet, y=ldmc))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("LDMC "~(mg ~g^-1)))+
  theme_classic()
ldmc.vic.p

# Leaf area
la.lath.p <- ggplot(data=comb.leaf.lath, aes(x=specific.epithet, y=leaf.area))+ 
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("LA "~(mm^2)))+
  theme_classic()
la.lath.p

la.phas.p <- ggplot(data=comb.leaf.phas, aes(x=specific.epithet, y=leaf.area))+ 
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("LA "~(mm^2)))+
  theme_classic()
la.phas.p

la.vic.p <- ggplot(data=comb.leaf.vic, aes(x=specific.epithet, y=leaf.area))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             pch=21, size=2, aes(fill=accession), show.legend = F)+
  labs(x="Species",y=bquote("LA "~(mm^2)))+
  theme_classic()
la.vic.p

# GRID of leaf plots
ggarrange(sla.lath.p, ldmc.lath.p, la.lath.p, sla.phas.p, ldmc.phas.p, la.phas.p, sla.vic.p, ldmc.vic.p, la.vic.p, nrow=3, ncol=3, widths=c(2,2,2), heights=c(2,2,2),
          labels=c("A","B","C","D","E","F","G","H","I"))

# Bivariate plots

sla.ldmc <- ggplot(data=comb.leaf, aes(x=sla, y=ldmc)) +
  geom_point(size=2, aes(fill=genus, shape=genus))+
  labs(x=bquote("SLA "~(mm^2 ~mg^-1)),y=bquote("LDMC "~(mg ~g^-1)))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  theme_classic()
sla.ldmc
summary(lm(data=comb.leaf, sla~ldmc))

la.ldmc <- ggplot(data=comb.leaf, aes(x=leaf.area, y=ldmc)) +
  geom_point(size=2, aes(fill=genus, shape=genus))+
  labs(x=bquote("LA "~(mm^2)),y=bquote("LDMC "~(mg ~g^-1)))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  theme_classic()
la.ldmc
summary(lm(data=comb.leaf, leaf.area~ldmc))

# Bivariate leaf grid
ggarrange(sla.ldmc, la.ldmc, nrow=1, ncol=2, widths=c(2,2), heights=c(2),
          labels=c("A","B"), common.legend=T, legend="bottom")

# Correlation matrices ####

# ALL DATA together
all.traits <- merge.acc.all %>% select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, nodes.21.mn, nodes.35.mn, nodes.perday.mn)

all.cor<-rcorr(as.matrix(all.traits), type="pearson") # all data
# Have to change column AND row names of original corr matrix to change them in corrplot
# For an rcorr object, have to change for both $r and $P
colnames(all.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(all.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
colnames(all.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(all.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(all.cor$r, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex=0.5, #Text label color and rotation
         number.cex=0.4, cl.cex = 0.6,
         # Combine with significance
         p.mat = all.cor$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

# SUBSETS of the data correlations

# ANNUAL correlations
ann.acc <- merge.acc.all[merge.acc.all$lifespan=="annual",]
ann.acc.traits <- ann.acc[, c(16,24,22,26,28,30,79,38,87,93,81,90,96,84)]
ann.cor<-rcorr(as.matrix(ann.acc.traits), type="pearson")
colnames(ann.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(ann.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
colnames(ann.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(ann.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(ann.cor$r, method="color", col=col(200), 
         type="upper", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex=0.5, #Text label color and rotation
         number.cex=0.4, cl.cex = 0.6,
         # Combine with significance
         p.mat = ann.cor$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

# PERENNIAL correlations
per.acc <- merge.acc.all[merge.acc.all$lifespan=="perennial",]
per.acc.traits <- per.acc[, c(16,24,22,26,28,30,79,38,87,93,81,90,96,84)]
per.cor<-rcorr(as.matrix(per.acc.traits), type="pearson")
colnames(per.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(per.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
colnames(per.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(per.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(per.cor$r, method="color", col=col(200),  
         type="upper", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex=0.5, #Text label color and rotation
         number.cex=0.4, cl.cex = 0.6,
         # Combine with significance
         p.mat = per.cor$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

# LATHYRUS correlations
lath.acc <- merge.acc.all[merge.acc.all$genus=="Lathyrus",]
lath.acc.traits <- lath.acc[, c(16,24,22,26,28,30,79,38,87,93,81,90,96,84)]
lath.cor<-rcorr(as.matrix(lath.acc.traits), type="pearson")
colnames(lath.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(lath.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
colnames(lath.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(lath.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(lath.cor$r, method="color", col=col(200),  
         type="upper", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex=0.5, #Text label color and rotation
         number.cex=0.4, cl.cex = 0.6,
         # Combine with significance
         p.mat = lath.cor$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

# PHASEOLUS correlations
phas.acc <- merge.acc.all[merge.acc.all$genus=="Phaseolus",]
phas.acc.traits <- phas.acc[, c(16,24,22,26,28,30,79,38,87,93,81,90,96,84)]
phas.cor<-rcorr(as.matrix(phas.acc.traits), type="pearson")
colnames(phas.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(phas.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
colnames(phas.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(phas.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(phas.cor$r, method="color", col=col(200),  
         type="upper", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex=0.5, #Text label color and rotation
         number.cex=0.4, cl.cex = 0.6,
         # Combine with significance
         p.mat = phas.cor$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

# VICIA correlations
vicia.acc <- merge.acc.all[merge.acc.all$genus=="Vicia",]
vicia.acc.traits <- vicia.acc[, c(16,24,22,26,28,30,79,38,87,93,81,90,96,84)]
vicia.cor<-rcorr(as.matrix(vicia.acc.traits), type="pearson")
colnames(vicia.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(vicia.cor$r) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
colnames(vicia.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
rownames(vicia.cor$P) <- c("Seed mass", "Seed length", "Seed width", "Seed perimeter", "Seed area", "Seed circularity", "Germination T50", "Germination proportion", "Height DAP-21", "Height DAP-35", "Height growth", "Nodes DAP-21", "Nodes DAP-35", "Node growth")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(vicia.cor$r, method="color", col=col(200),  
         type="upper", order="original", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex=0.5, #Text label color and rotation
         number.cex=0.4, cl.cex = 0.6,
         # Combine with significance
         p.mat = vicia.cor$P, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

##### Experimenting

# LINEAR MODEL experimenting ####

# T50 by seed area
all <- lm(data=merge.acc.all, seed.area.mn~ppd50.mn)
summary(all)
lath.ann.t <- lm(data=lath.ann.merge, seed.area.mn~ppd50.mn)
summary(lath.ann.t)
lath.per.t <- lm(data=lath.per.merge, seed.area.mn~ppd50.mn)
summary(lath.per.t)
phas.ann.t <- lm(data=phas.ann.merge, seed.area.mn~ppd50.mn)
summary(phas.ann.t)
phas.per.t <- lm(data=phas.per.merge, seed.area.mn~ppd50.mn)
summary(phas.per.t)
vicia.ann.t <- lm(data=vicia.ann.merge, seed.area.mn~ppd50.mn) #sig
summary(vicia.ann.t)
vicia.per.t <- lm(data=vicia.per.merge, seed.area.mn~ppd50.mn)
summary(vicia.per.t)

plot(data=merge.acc.all, ppd50.mn~seed.circ.mn)

# T50 by seed circularity
all <- lm(data=merge.acc.all, seed.circ.mn~ppd50.mn) # sig
summary(all)
lath.ann.t <- lm(data=lath.ann.merge, seed.area.mn~ppd50.mn)
summary(lath.ann.t)
lath.per.t <- lm(data=lath.per.merge, seed.area.mn~ppd50.mn)
summary(lath.per.t)
phas.ann.t <- lm(data=phas.ann.merge, seed.area.mn~ppd50.mn)
summary(phas.ann.t)
phas.per.t <- lm(data=phas.per.merge, seed.area.mn~ppd50.mn)
summary(phas.per.t)
vicia.ann.t <- lm(data=vicia.ann.merge, seed.area.mn~ppd50.mn) #sig
summary(vicia.ann.t)
vicia.per.t <- lm(data=vicia.per.merge, seed.area.mn~ppd50.mn)
summary(vicia.per.t)

# Viability proportion by seed circularity
all <- lm(data=merge.acc.all, seed.circ.mn~FINAL.viable.prop) #sig
summary(all)
lath.ann.t <- lm(data=lath.ann.merge, seed.circ.mn~FINAL.viable.prop)
summary(lath.ann.t)
lath.per.t <- lm(data=lath.per.merge, seed.circ.mn~FINAL.viable.prop)
summary(lath.per.t)
phas.ann.t <- lm(data=phas.ann.merge, seed.circ.mn~FINAL.viable.prop)
summary(phas.ann.t)
phas.per.t <- lm(data=phas.per.merge, seed.circ.mn~FINAL.viable.prop)
summary(phas.per.t)
vicia.ann.t <- lm(data=vicia.ann.merge, seed.circ.mn~FINAL.viable.prop)
summary(vicia.ann.t)
vicia.per.t <- lm(data=vicia.per.merge, seed.circ.mn~FINAL.viable.prop)
summary(vicia.per.t)

# Germination proportion by seed circularity
all <- lm(data=merge.acc.all, seed.circ.mn~FINAL.germ.prop) #sig
summary(all)
lath.ann.t <- lm(data=lath.ann.merge, seed.circ.mn~FINAL.germ.prop)
summary(lath.ann.t)
lath.per.t <- lm(data=lath.per.merge, seed.circ.mn~FINAL.germ.prop)
summary(lath.per.t)
phas.ann.t <- lm(data=phas.ann.merge, seed.circ.mn~FINAL.germ.prop)
summary(phas.ann.t)
phas.per.t <- lm(data=phas.per.merge, seed.circ.mn~FINAL.germ.prop)
summary(phas.per.t)
vicia.ann.t <- lm(data=vicia.ann.merge, seed.circ.mn~FINAL.germ.prop)
summary(vicia.ann.t)
vicia.per.t <- lm(data=vicia.per.merge, seed.circ.mn~FINAL.germ.prop)
summary(vicia.per.t)

# Viability proportion by seed area
all <- lm(data=merge.acc.all, seed.area.mn~FINAL.viable.prop)
summary(all)
lath.ann.t <- lm(data=lath.ann.merge, seed.area.mn~FINAL.viable.prop)
summary(lath.ann.t)
lath.per.t <- lm(data=lath.per.merge, seed.area.mn~FINAL.viable.prop)
summary(lath.per.t)
phas.ann.t <- lm(data=phas.ann.merge, seed.area.mn~FINAL.viable.prop)
summary(phas.ann.t)
phas.per.t <- lm(data=phas.per.merge, seed.area.mn~FINAL.viable.prop)
summary(phas.per.t)
vicia.ann.t <- lm(data=vicia.ann.merge, seed.area.mn~FINAL.viable.prop)
summary(vicia.ann.t)
vicia.per.t <- lm(data=vicia.per.merge, seed.area.mn~FINAL.viable.prop)
summary(vicia.per.t)

# Germination proportion by seed area
all <- lm(data=merge.acc.all, seed.area.mn~FINAL.germ.prop) #sig
summary(all)
lath.ann.t <- lm(data=lath.ann.merge, seed.area.mn~FINAL.germ.prop)
summary(lath.ann.t)
lath.per.t <- lm(data=lath.per.merge, seed.area.mn~FINAL.germ.prop)
summary(lath.per.t)
phas.ann.t <- lm(data=phas.ann.merge, seed.area.mn~FINAL.germ.prop)
summary(phas.ann.t)
phas.per.t <- lm(data=phas.per.merge, seed.area.mn~FINAL.germ.prop)
summary(phas.per.t)
vicia.ann.t <- lm(data=vicia.ann.merge, seed.area.mn~FINAL.germ.prop)
summary(vicia.ann.t)
vicia.per.t <- lm(data=vicia.per.merge, seed.area.mn~FINAL.germ.prop)
summary(vicia.per.t)
phas.ann.t <- lm(data=phas.ann.merge, seed.area.mn~height.perday.mn)
summary(phas.ann.t)
phas.per.t <- lm(data=phas.per.merge, seed.area.mn~height.perday.mn)
summary(phas.per.t)

####
# Getting regression line text ####
# To get full regression line equation in the plots
lm_eqn <- function(df, y, x){
  formula = as.formula(sprintf('%s ~ %s', y, x))
  m <- lm(formula, data=df);
  # formating the values into a summary string to print out
  # ~ give some space, but equal size and comma need to be quoted
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
                   list(a = format(as.vector(coef(m)[1]), digits = 2), 
                        b = format(as.vector(coef(m)[2]), digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3),
                        # getting the pvalue is painful
                        pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
                   )
  )
  as.character(as.expression(eq));                 
}
# Seed discussion in
# https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph

# regression line with only r2 and p values
lm_eqn1 <- function(df, y, x){
  formula = as.formula(sprintf('%s ~ %s', y, x))
  m <- lm(formula, data=df);
  # formating the values into a summary string to print out
  # ~ give some space, but equal size and comma need to be quoted
  eq <- substitute(~~italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
                   list(r2 = format(summary(m)$r.squared, digits = 3),
                        # getting the pvalue is painful
                        pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
                   )
  )
  as.character(as.expression(eq));   
}

# BIVARIATE plots - raw seed data ####

# Subsetting
lath.comb.seed <- comb.seed.clean[comb.seed.clean$genus=="Lathyrus",]
phas.comb.seed <- comb.seed.clean[comb.seed.clean$genus=="Phaseolus",]
vicia.comb.seed <- comb.seed.clean[comb.seed.clean$genus=="Vicia",]
lath.ann.seed <- comb.seed.clean[with(comb.seed.clean, genus == "Lathyrus" & lifespan == "annual"),]
lath.per.seed <- comb.seed.clean[with(comb.seed.clean, genus == "Lathyrus" & lifespan == "perennial"),]
phas.ann.seed <- comb.seed.clean[with(comb.seed.clean, genus == "Phaseolus" & lifespan == "annual"),]
phas.per.seed <- comb.seed.clean[with(comb.seed.clean, genus == "Phaseolus" & lifespan == "perennial"),]
vicia.ann.seed <- comb.seed.clean[with(comb.seed.clean, genus == "Vicia" & lifespan == "annual"),]
vicia.per.seed <- comb.seed.clean[with(comb.seed.clean, genus == "Vicia" & lifespan == "perennial"),]

ggplot(data=phas.per.seed, aes(x=seed.circ, y=seed.area, color=species.short))+
  geom_point(size=3)+
  labs(title="Perennial Lathyrus")


lath.ann.lm <- lm_eqn1(lath.ann.seed, "seed.circ", "seed.area")
lath.per.lm <- lm_eqn1(lath.per.seed, "seed.circ", "seed.area")
test <- lm(data=lath.ann.seed, seed.circ~seed.area)
p.lath.seed.area.x.seed.circ <- ggplot(data=lath.comb.seed, aes(x=seed.area, y=seed.circ, 
                            group=lifespan, shape=lifespan, fill=lifespan, color=lifespan))+
  geom_smooth(method=lm, se=F, fullrange=T)+
  geom_point(size=3, stroke=1, color="black")+
  labs(x="Seed area",y="Seed circularity", title="Lathyrus")+
  scale_shape_manual(name = "Lifespan",
                     labels = c("Annual","Perennial"),
                     values = c(21, 24))+
  scale_fill_manual(name = "Lifespan",
                    labels = c("Annual", "Perennial"),
                    values = c("firebrick3", "mediumblue")) +
  scale_colour_manual(name = "Lifespan",
                      labels = c("Annual", "Perennial"),
                      values = c("firebrick3", "mediumblue")) +
  theme_classic()
p.lath.seed.area.x.seed.circ

# BIVARIATE plots - accession level ####

lath.merge <- merge.acc.all[merge.acc.all$genus=="Lathyrus",]
phas.merge <- merge.acc.all[merge.acc.all$genus=="Phaseolus",]
vicia.merge <- merge.acc.all[merge.acc.all$genus=="Vicia",]
lath.ann.merge <- merge.acc.all[with(merge.acc.all, genus == "Lathyrus" & lifespan == "annual"),]
lath.per.merge <- merge.acc.all[with(merge.acc.all, genus == "Lathyrus" & lifespan == "perennial"),]
phas.ann.merge <- merge.acc.all[with(merge.acc.all, genus == "Phaseolus" & lifespan == "annual"),]
phas.per.merge <- merge.acc.all[with(merge.acc.all, genus == "Phaseolus" & lifespan == "perennial"),]
vicia.ann.merge <- merge.acc.all[with(merge.acc.all, genus == "Vicia" & lifespan == "annual"),]
vicia.per.merge <-merge.acc.all[with(merge.acc.all, genus == "Vicia" & lifespan == "perennial"),]

# NEW, simple visualizing of linear relationships

# generic bivariate figure code for within group
biv.fig <- ggplot(data=lath.acc1, aes(x=seed.area.mn, y=seed.width.mn))+
  geom_point(aes(stroke=1, color=species.short, shape=species.short), size=2)+
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9,10))+
  geom_smooth(method=lm, se=F, fullrange=T)+
  theme_classic()
biv.fig

# bivariate correlations for annual/perennial dataset
biv.lifespan <- ggplot(data=per.acc1, aes(x=seed.area.mn, y=leaves.rgr.mn))+
  geom_point(aes(color=genus, shape=genus, alpha=0.9), size=3)+
  geom_smooth(method=lm, se=F, fullrange=T)+
  theme_classic()
biv.lifespan

# bivariate correlations for full dataset
biv.full <- ggplot(data=merge.acc.all, aes(x=seed.perim.mn, y=height.21.mn))+
  geom_point(aes(color=genus, shape=lifespan, alpha=0.75), size=3)+
  geom_smooth(method=lm, se=F, fullrange=T)+
  theme_classic()
biv.full

# checking out a bunch of correlations at a time
# need to run correlation analyses above 
library("PerformanceAnalytics")
chart.Correlation(xv, pch=19)


# A few examples with regression stats put in plot via ggplot
# Lathyrus seed area x circularity
lath.ann.lm <- lm_eqn1(lath.ann.merge, "seed.circ.mn", "seed.area.mn")
lath.per.lm <- lm_eqn1(lath.per.merge, "seed.circ.mn", "seed.area.mn")
p.lath.seed.area.x.seed.circ <- ggplot(data=lath.merge, aes(x=seed.area.mn, y=seed.circ.mn, 
                                    group=lifespan, shape=lifespan, fill=lifespan, color=lifespan))+
  geom_smooth(method=lm, se=F, fullrange=T)+
  geom_point(size=3, stroke=1, color="black")+
  geom_text(x=15,y=0.913,label=lath.ann.lm, color='firebrick3', parse=T, size=3)+
  geom_text(x=15,y=0.91,label=lath.per.lm, color='mediumblue', parse=T, size=3)+
  labs(x="Seed area",y="Seed circularity", title="Lathyrus")+
  scale_shape_manual(name = "Lifespan",
                     labels = c("Annual","Perennial"),
                     values = c(21, 24))+
  scale_fill_manual(name = "Lifespan",
                    labels = c("Annual", "Perennial"),
                    values = c("firebrick3", "mediumblue")) +
  scale_colour_manual(name = "Lifespan",
                    labels = c("Annual", "Perennial"),
                    values = c("firebrick3", "mediumblue")) +
  theme_classic()
p.lath.seed.area.x.seed.circ

# Phaseolus seed area x circularity
phas.ann.lm <- lm_eqn1(phas.ann.merge, "seed.circ.mn", "seed.area.mn")
phas.per.lm <- lm_eqn1(phas.per.merge, "seed.circ.mn", "seed.area.mn")
p.phas.seed.area.x.seed.circ <- ggplot(data=phas.merge, aes(x=seed.area.mn, y=seed.circ.mn, 
                              group=lifespan, shape=lifespan, fill=lifespan, color=lifespan))+
  geom_smooth(method=lm, se=F, fullrange=T)+
  geom_point(size=3, stroke=1, color="black")+
  geom_text(x=35,y=0.88,label=phas.ann.lm, color='firebrick3', parse=T, size=3)+
  geom_text(x=35,y=0.877,label=phas.per.lm, color='mediumblue', parse=T, size=3)+
  labs(x="Seed area",y="Seed circularity", title="Phaseolus")+
  scale_shape_manual(name = "Lifespan",
                     labels = c("Annual","Perennial"),
                     values = c(21, 24))+
  scale_fill_manual(name = "Lifespan",
                    labels = c("Annual", "Perennial"),
                    values = c("firebrick3", "mediumblue")) +
  scale_colour_manual(name = "Lifespan",
                      labels = c("Annual", "Perennial"),
                      values = c("firebrick3", "mediumblue")) +
  theme_classic()
p.phas.seed.area.x.seed.circ

# Vicia seed area x circularity
vicia.ann.lm <- lm_eqn1(vicia.ann.merge, "seed.circ.mn", "seed.area.mn")
vicia.per.lm <- lm_eqn1(vicia.per.merge, "seed.circ.mn", "seed.area.mn")
p.vicia.seed.area.x.seed.circ <- ggplot(data=vicia.merge, aes(x=seed.area.mn, y=seed.circ.mn, 
                              group=lifespan, shape=lifespan, fill=lifespan, color=lifespan))+
  geom_smooth(method=lm, se=F, fullrange=T)+
  geom_point(size=3, stroke=1, color="black")+
  geom_text(x=14,y=0.92,label=vicia.ann.lm, color='firebrick3', parse=T, size=3)+
  geom_text(x=14,y=0.917,label=vicia.per.lm, color='mediumblue', parse=T, size=3)+
  labs(x="Seed area",y="Seed circularity", title="Vicia")+
  scale_shape_manual(name = "Lifespan",
                     labels = c("Annual","Perennial"),
                     values = c(21, 24))+
  scale_fill_manual(name = "Lifespan",
                    labels = c("Annual", "Perennial"),
                    values = c("firebrick3", "mediumblue")) +
  scale_colour_manual(name = "Lifespan",
                      labels = c("Annual", "Perennial"),
                      values = c("firebrick3", "mediumblue")) +
  theme_classic()
p.vicia.seed.area.x.seed.circ

# seed size by germination rate
vicia.ann.lm <- lm_eqn1(vicia.ann.merge, "ppd50.mn", "seed.area.mn")
vicia.per.lm <- lm_eqn1(vicia.per.merge, "ppd50.mn", "seed.area.mn")
p.vicia.seed.area.x.ppd50 <- ggplot(data=vicia.merge, aes(x=seed.area.mn, y=ppd50.mn, 
                                                              group=lifespan, shape=lifespan, fill=lifespan, color=lifespan))+
  geom_smooth(method=lm, se=F, fullrange=T)+
  geom_point(size=3, stroke=1, color="black")+
  geom_text(x=14,y=11,label=vicia.ann.lm, color='firebrick3', parse=T, size=3)+
  geom_text(x=14,y=10,label=vicia.per.lm, color='mediumblue', parse=T, size=3)+
  labs(x="Seed area",y="T50 germination", title="Vicia")+
  scale_shape_manual(name = "Lifespan",
                     labels = c("Annual","Perennial"),
                     values = c(21, 24))+
  scale_fill_manual(name = "Lifespan",
                    labels = c("Annual", "Perennial"),
                    values = c("firebrick3", "mediumblue")) +
  scale_colour_manual(name = "Lifespan",
                      labels = c("Annual", "Perennial"),
                      values = c("firebrick3", "mediumblue")) +
  theme_classic()
p.vicia.seed.area.x.ppd50

# Vicia ppd50 by seed circularity
vicia.ann.lm <- lm_eqn1(vicia.ann.merge, "ppd50.mn", "seed.circ.mn")
vicia.per.lm <- lm_eqn1(vicia.per.merge, "ppd50.mn", "seed.circ.mn")
p.vicia.seed.area.x.ppd50 <- ggplot(data=vicia.merge, aes(x=seed.circ.mn, y=ppd50.mn, 
                            group=lifespan, shape=lifespan, fill=lifespan, color=lifespan))+
  geom_smooth(method=lm, se=F, fullrange=T)+
  geom_point(size=3, stroke=1, color="black")+
  geom_text(x=14,y=11,label=vicia.ann.lm, color='firebrick3', parse=T, size=3)+
  geom_text(x=14,y=10,label=vicia.per.lm, color='mediumblue', parse=T, size=3)+
  labs(x="Seed circularity",y="T50 germination", title="Vicia")+
  scale_shape_manual(name = "Lifespan",
                     labels = c("Annual","Perennial"),
                     values = c(21, 24))+
  scale_fill_manual(name = "Lifespan",
                    labels = c("Annual", "Perennial"),
                    values = c("firebrick3", "mediumblue")) +
  scale_colour_manual(name = "Lifespan",
                      labels = c("Annual", "Perennial"),
                      values = c("firebrick3", "mediumblue")) +
  theme_classic()
p.vicia.seed.area.x.ppd50

# SPECIES BOXPLOTS ####

# For traits with a significant species effect (So, all of them)

# Lathyrus
seed.wt.p <- ggplot(data=lath.merge, aes(x=reorder(species.short,seed.wt.avg,median, na.rm=T), y=seed.wt.avg, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.wt.p

seed.area.p <- ggplot(data=lath.merge, aes(x=reorder(species.short,seed.area.mn,median, na.rm=T), y=seed.area.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.area.p

seed.circ.p <- ggplot(data=lath.merge, aes(x=reorder(species.short,seed.circ.mn,median, na.rm=T), y=seed.circ.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.circ.p

ppd50.p <- ggplot(data=lath.merge, aes(x=reorder(species.short,ppd50.mn,median, na.rm=T), y=ppd50.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
ppd50.p

viable.prop.p <- ggplot(data=lath.merge, aes(x=reorder(species.short,FINAL.viable.prop,median, na.rm=T), y=FINAL.viable.prop, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
viable.prop.p

height.growth.p <- ggplot(data=lath.merge, aes(x=reorder(species.short,height.perday.mn,median, na.rm=T), y=height.perday.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
height.growth.p

# Phaseolus
seed.wt.p <- ggplot(data=phas.merge, aes(x=reorder(species.short,seed.wt.avg,median, na.rm=T), y=seed.wt.avg, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.wt.p

seed.area.p <- ggplot(data=phas.merge, aes(x=reorder(species.short,seed.area.mn,median, na.rm=T), y=seed.area.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.area.p

seed.circ.p <- ggplot(data=phas.merge, aes(x=reorder(species.short,seed.circ.mn,median, na.rm=T), y=seed.circ.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.circ.p

ppd50.p <- ggplot(data=phas.merge, aes(x=reorder(species.short,ppd50.mn,median, na.rm=T), y=ppd50.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
ppd50.p

germ.prop.p <- ggplot(data=phas.merge, aes(x=reorder(species.short,FINAL.germ.prop,median, na.rm=T), y=FINAL.germ.prop, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
germ.prop.p

height.growth.p <- ggplot(data=phas.merge, aes(x=reorder(species.short,height.perday.mn,median, na.rm=T), y=height.perday.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
height.growth.p

# Vicia
seed.wt.p <- ggplot(data=vicia.merge, aes(x=reorder(species.short,seed.wt.avg,median, na.rm=T), y=seed.wt.avg, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.wt.p

seed.area.p <- ggplot(data=vicia.merge, aes(x=reorder(species.short,seed.area.mn,median, na.rm=T), y=seed.area.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.area.p

seed.circ.p <- ggplot(data=vicia.merge, aes(x=reorder(species.short,seed.circ.mn,median, na.rm=T), y=seed.circ.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
seed.circ.p

ppd50.p <- ggplot(data=vicia.merge, aes(x=reorder(species.short,ppd50.mn,median, na.rm=T), y=ppd50.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
ppd50.p

germ.prop.p <- ggplot(data=vicia.merge, aes(x=reorder(species.short,FINAL.germ.prop,median, na.rm=T), y=FINAL.germ.prop, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
germ.prop.p

height.growth.p <- ggplot(data=vicia.merge, aes(x=reorder(species.short,height.perday.mn,median, na.rm=T), y=height.perday.mn, color=lifespan))+ 
  geom_boxplot()+
  geom_point()
height.growth.p

# reproductive data summaries and concatenation ####

# sum of info from 'packet info' spreadsheet
# use edited one in 2017 folder - fixed inaccuracies in original

repro_effort2017 <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/Original data files/repro_effort2017.xlsx", 
                               sheet = "packet info", na="NA")
View(repro_effort2017)

packet.sum <- repro_effort2017 %>%
  group_by(accession, replicate, individual) %>%
  summarise(    
    date.min = min(date, na.rm=T),
    date.max = max(date, na.rm=T),
    ripe.sum = sum(ripe.pod, na.rm=T),
    dehisc.sum = sum(dehisc.pod, na.rm=T),
    undev.sum = sum(undev.pod, na.rm=T)
  ) %>% ungroup ()
View(packet.sum)

# Concatenate repro counts, weights, and packet info
# By accession, rep, individual

count.2017 <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/Original data files/Reproductive counts & weights.xlsx", 
                         sheet = "counts", na="NA")
View(count.2017)

# merge count & weight data frames
weight.2017 <- read_excel("C:/Users/SterlingH/Desktop/DATA_GOOGLE/2017 STL/DATA/Original data files/Reproductive counts & weights.xlsx", 
                         sheet = "weights", na="NA")
View(weight.2017)

merge1 <- merge(count.2017, weight.2017, by=c("accession","replicate","individual"))
View(merge1)

# merge previous merge with packet info

merge.final <- merge(merge1, packet.sum, by=c("accession","replicate","individual"))
View(merge.final)

# write as excel
write.xlsx(merge.final, file = "C:/Users/SterlingH/Desktop/R_exports/repro.merge.xlsx")

# merge reproductive count / phenology data ####
# combining data for each individual to get the most reproductive information.
# ORIGINALS!!! ####
p.sw <- ggplot(data=seedwt.fig, aes(x=lifespan, y=seed.wt.mn, group=genus, shape=genus, color=genus, fill=genus))
sw.plot <- p.sw + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=seed.wt.mn-seed.wt.sd, ymax=seed.wt.mn+seed.wt.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Average seed mass (mg)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sw.plot

# Seed length
p.sl <- ggplot(data=seed.fig, aes(x=lifespan, y=seed.length.mn, group=genus, shape=genus, color=genus, fill=genus))
sl.plot <- p.sl + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=seed.length.mn-seed.length.sd, ymax=seed.length.mn+seed.length.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Seed length (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sl.plot

# Seed width
p.swd <- ggplot(data=seed.fig, aes(x=lifespan, y=seed.width.mn, group=genus, shape=genus, color=genus, fill=genus))
swd.plot <- p.swd + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=seed.width.mn-seed.width.sd, ymax=seed.width.mn+seed.width.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Seed width (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
swd.plot

# Seed perimeter
p.sp <- ggplot(data=seed.fig, aes(x=lifespan, y=seed.perim.mn, group=genus, shape=genus, color=genus, fill=genus))
sp.plot <- p.sp + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=seed.perim.mn-seed.perim.sd, ymax=seed.perim.mn+seed.perim.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y=bquote("Seed perimeter (mm)"))+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sp.plot

# Seed area
p.sa <- ggplot(data=seed.fig, aes(x=lifespan, y=seed.area.mn, group=genus, shape=genus, color=genus, fill=genus))
sa.plot <- p.sa + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=seed.area.mn-seed.area.sd, ymax=seed.area.mn+seed.area.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y=bquote("Seed area "~(mm^2)))+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sa.plot

# Seed circularity
p.sc <- ggplot(data=seed.fig, aes(x=lifespan, y=seed.circ.mn, group=genus, shape=genus, color=genus, fill=genus))
sc.plot <- p.sc + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=seed.circ.mn-seed.circ.sd, ymax=seed.circ.mn+seed.circ.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Seed circularity")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  # easier to visualize the difference in circularity w/o y starting at 0, since it is a much smaller range of variation
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
sc.plot

# Germination proportion
p.gp <- ggplot(data=germ.fig, aes(x=lifespan, y=germ.prop.mn, group=genus, shape=genus, color=genus, fill=genus))
gp.plot <- p.gp + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=germ.prop.mn-germ.prop.sd, ymax=germ.prop.mn+germ.prop.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Germination proportion")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  # easier to visualize the difference in circularity w/o y starting at 0, since it is a much smaller range of variation
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
gp.plot

# Viable proportion
p.vp <- ggplot(data=germ.fig, aes(x=lifespan, y=viab.prop.mn, group=genus, shape=genus, color=genus, fill=genus))
vp.plot <- p.vp + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=viab.prop.mn-viab.prop.sd, ymax=viab.prop.mn+viab.prop.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Viable proportion")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
vp.plot

# Germination rate (T50)
p.gr <- ggplot(data=germrate.fig, aes(x=lifespan, y=ppd50.fig.mn, group=genus, shape=genus, color=genus, fill=genus))
gr.plot <- p.gr + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=ppd50.fig.mn-ppd50.fig.sd, ymax=ppd50.fig.mn+ppd50.fig.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="T50 (days)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+ # forces y axis to start at 0
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
gr.plot

# Height growth
p.ht <- ggplot(data=growth.fig, aes(x=lifespan, y=height.perday.mn, group=genus, shape=genus, color=genus, fill=genus))
ht.plot <- p.ht + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=height.perday.mn-height.perday.sd, ymax=height.perday.mn+height.perday.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Growth rate (mm/day)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
ht.plot

# Nodes growth
p.nd <- ggplot(data=growth.fig, aes(x=lifespan, y=nodes.perday.mn, group=genus, shape=genus, color=genus, fill=genus))
nd.plot <- p.nd + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=nodes.perday.mn-nodes.perday.sd, ymax=nodes.perday.mn+nodes.perday.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Growth rate (nodes/day)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
nd.plot

# Height 21 days
p.ht.21 <- ggplot(data=growth.fig, aes(x=lifespan, y=height.21.mn, group=genus, shape=genus, color=genus, fill=genus))
ht.21.plot <- p.ht.21 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=height.21.mn-height.21.sd, ymax=height.21.mn+height.21.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Height at 21 days (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
ht.21.plot

# Nodes 21 days
p.nd.21 <- ggplot(data=growth.fig, aes(x=lifespan, y=nodes.21.mn, group=genus, shape=genus, color=genus, fill=genus))
nodes.21.plot <- p.nd.21 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=nodes.21.mn-nodes.21.sd, ymax=nodes.21.mn+nodes.21.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Node number at 21 days")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
nodes.21.plot

# Height 35 days
p.ht.35 <- ggplot(data=growth.fig, aes(x=lifespan, y=height.35.mn, group=genus, shape=genus, color=genus, fill=genus))
ht.35.plot <- p.ht.35 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=height.35.mn-height.35.sd, ymax=height.35.mn+height.35.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Height at 35 days (mm)")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
ht.35.plot

# Nodes 35 days
p.nd.35 <- ggplot(data=growth.fig, aes(x=lifespan, y=nodes.35.mn, group=genus, shape=genus, color=genus, fill=genus))
nodes.35.plot <- p.nd.35 + 
  geom_line(position=position_dodge(.4), size=.75)+
  geom_errorbar(aes(ymin=nodes.35.mn-nodes.35.sd, ymax=nodes.35.mn+nodes.35.sd),
                width=0.2, position = position_dodge(.4), color="black")+
  geom_point(size=3, stroke=1, color="black", position=position_dodge(.4), stat="identity")+
  labs(x="Lifespan",y="Node number at 35 days")+
  scale_x_discrete(labels=c("annual"="Annual", "perennial"="Perennial"))+
  scale_shape_manual(name = "Genus",
                     labels = c("Lathyrus", "Phaseolus", "Vicia"),
                     values = c(21, 22, 24))+
  scale_fill_manual(name = "Genus",
                    labels = c("Lathyrus", "Phaseolus", "Vicia"),
                    values = c("midnightblue", "grey50", "dodgerblue"))+
  scale_colour_manual(name = "Genus",
                      labels = c("Lathyrus", "Phaseolus", "Vicia"),
                      values = c("midnightblue", "grey50", "dodgerblue"))+
  expand_limits(y = 0)+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        panel.background = element_blank())
nodes.35.plot
# Boxplot FIGURES - trait by trait ####

# seed mass
ggplot(data=lath, aes(x=lifespan, y=seed.wt.avg, fill=lifespan))+
  geom_boxplot()+
  labs(title="Lathyrus",x="Lifespan",y="Average single seed mass (mg)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=phas, aes(x=lifespan, y=seed.wt.avg, fill=lifespan))+
  geom_boxplot()+
  labs(title="Phaseolus",x="Lifespan",y="Average single seed mass (mg)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=vic, aes(x=lifespan, y=seed.wt.avg, fill=lifespan))+
  geom_boxplot()+
  labs(title="Vicia",x="Lifespan",y="Average single seed mass (mg)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

p.sw <- ggplot(data=comb, aes(x=lifespan, y=seed.wt.avg, fill=lifespan))
sw.plot <- p.sw + 
  geom_boxplot()+
  labs(x="Lifespan",y="Average single seed mass (mg)")+
  scale_fill_manual(values = c("firebrick","dodgerblue4"))+
  theme_classic()+
  facet_wrap(~genus)
sw.plot

# seed area
ggplot(data=lath, aes(x=lifespan, y=seed.area, fill=lifespan))+
  geom_boxplot()+
  labs(title="Lathyrus",x="Lifespan",y=bquote("Single seed area "~(mm^2)))+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=phas, aes(x=lifespan, y=seed.area, fill=lifespan))+
  geom_boxplot()+
  labs(title="Phaseolus",x="Lifespan",y=bquote("Single seed area "~(mm^2)))+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=vic, aes(x=lifespan, y=seed.area, fill=lifespan))+
  geom_boxplot()+
  labs(title="Vicia",x="Lifespan",y=bquote("Single seed area "~(mm^2)))+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

# seed width
ggplot(data=lath, aes(x=lifespan, y=seed.width, fill=lifespan))+
  geom_boxplot()+
  labs(title="Lathyrus",x="Lifespan",y="Single seed width (mm)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=phas, aes(x=lifespan, y=seed.width, fill=lifespan))+
  geom_boxplot()+
  labs(title="Phaseolus",x="Lifespan",y="Single seed width (mm)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=vic, aes(x=lifespan, y=seed.width, fill=lifespan))+
  geom_boxplot()+
  labs(title="Vicia",x="Lifespan",y="Single seed width (mm)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

# seed length
ggplot(data=lath, aes(x=lifespan, y=seed.length, fill=lifespan))+
  geom_boxplot()+
  labs(title="Lathyrus",x="Lifespan",y="Single seed length (mm)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=phas, aes(x=lifespan, y=seed.length, fill=lifespan))+
  geom_boxplot()+
  labs(title="Phaseolus",x="Lifespan",y="Single seed length (mm)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=vic, aes(x=lifespan, y=seed.length, fill=lifespan))+
  geom_boxplot()+
  labs(title="Vicia",x="Lifespan",y="Single seed length (mm)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

# seed circularity
ggplot(data=lath, aes(x=lifespan, y=seed.circ, fill=lifespan))+
  geom_boxplot()+
  labs(title="Lathyrus",x="Lifespan",y="Single seed circularity")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=phas, aes(x=lifespan, y=seed.circ, fill=lifespan))+
  geom_boxplot()+
  labs(title="Phaseolus",x="Lifespan",y="Single seed circularity")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=vic, aes(x=lifespan, y=seed.circ, fill=lifespan))+
  geom_boxplot()+
  labs(title="Vicia",x="Lifespan",y="Single seed circularity")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

# growth rate (height)

ggplot(data=comb.plant[comb.plant$genus=="Lathyrus",], aes(x=lifespan, y=height.perday.21.35, fill=lifespan))+
  geom_boxplot()+
  labs(title="Lathyrus",x="Lifespan",y="Growth rate (mm / day)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=comb.plant[comb.plant$genus=="Phaseolus",], aes(x=lifespan, y=height.perday.21.35, fill=lifespan))+
  geom_boxplot()+
  labs(title="Phaseolus",x="Lifespan",y="Growth rate (mm / day)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=comb.plant[comb.plant$genus=="Vicia",], aes(x=lifespan, y=height.perday.21.35, fill=lifespan))+
  geom_boxplot()+
  labs(title="Vicia",x="Lifespan",y="Growth rate (mm / day)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

# growth rate (nodes)

ggplot(data=comb.plant[comb.plant$genus=="Lathyrus",], aes(x=lifespan, y=nodes.perday.21.35, fill=lifespan))+
  geom_boxplot()+
  labs(title="Lathyrus",x="Lifespan",y="Growth rate (nodes / day)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=comb.plant[comb.plant$genus=="Phaseolus",], aes(x=lifespan, y=nodes.perday.21.35, fill=lifespan))+
  geom_boxplot()+
  labs(title="Phaseolus",x="Lifespan",y="Growth rate (nodes / day)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

ggplot(data=comb.plant[comb.plant$genus=="Vicia",], aes(x=lifespan, y=nodes.perday.21.35, fill=lifespan))+
  geom_boxplot()+
  labs(title="Vicia",x="Lifespan",y="Growth rate (nodes / day)")+
  scale_fill_manual(values = c("red3","royalblue3"))+
  theme_classic()

# LINEAR models - SEED TRAITS

Anova(lm(data=comb.acc, seed.wt.avg~genus + lifespan + genus*lifespan 
         + genus/lifespan/species.short), 
      type = "II")

Anova(lm(data=comb.acc, seed.area~genus + lifespan + genus*lifespan 
         + genus/lifespan/species.short), 
      type = "II")

Anova(lm(data=comb.acc, seed.length~genus + lifespan + genus*lifespan 
         + genus/lifespan/species.short), 
      type = "II")

Anova(lm(data=comb.acc, seed.width~genus + lifespan + genus*lifespan 
         + genus/lifespan/species.short), 
      type = "II")

Anova(lm(data=comb.acc, seed.circ~genus + lifespan + genus*lifespan 
         + genus/lifespan/species.short), 
      type = "II")


# GENUS-BY-GENUS GERMINATION MODELS ####

# GERMINATION T50

# Subset since need separate models for each genus
lath.comb.germ <- comb.germ.clean[comb.germ.clean$genus=="Lathyrus",]
phas.comb.germ <- comb.germ.clean[comb.germ.clean$genus=="Phaseolus",]
vic.comb.germ <- comb.germ.clean[comb.germ.clean$genus=="Vicia",]

# Lathyrus
# ORIGINAL
model <- lmer(data=lath.comb.germ, ppd50~lifespan + lifespan/species.short
              + (1 | species.short:accession) + (1 | age) + (1 | sandpaper))

# FINAL
# Nonsignificant: age, sandpaper; keeping age, since biologically relevant
lath.gmrt.modreduced <- lmer(data=lath.comb.germ, ppd50~lifespan + lifespan/species.short
                             + (1 | species.short:accession) + (1 | age))

rand(model)
rand(lath.gmrt.modreduced)
anova(lath.gmrt.modreduced)
emmeans(object=lath.gmrt.modreduced, specs = pairwise ~ lifespan) # gives you the p value
lath.gmrt.sum <- emmeans(object=lath.gmrt.modreduced, specs = c("lifespan")) # gives you a workable structure for a dataframe
lath.gmrt.df <- as.data.frame(lath.gmrt.sum)

# Phaseolus
# ORIGINAL
model <- lmer(data=phas.comb.germ, ppd50~lifespan + lifespan/species.short
              + (1 | species.short:accession) + (1 | age) + (1 | sandpaper))

# FINAL
# Nonsignificant: age, sandpaper
phas.gmrt.modreduced <- lmer(data=phas.comb.germ, ppd50~lifespan + lifespan/species.short
                             + (1 | species.short:accession) + (1 | age))

rand(model)
rand(phas.gmrt.modreduced)
anova(phas.gmrt.modreduced)
emmeans(object=phas.gmrt.modreduced, specs = pairwise ~ lifespan) # gives you the p value
phas.gmrt.sum <- emmeans(object=phas.gmrt.modreduced, specs = c("lifespan")) # gives you a workable structure for a dataframe
phas.gmrt.df <- as.data.frame(phas.gmrt.sum)

# Vicia
# ORIGINAL
model <- lmer(data=vic.comb.germ, ppd50~lifespan + lifespan/species.short
              + (1 | species.short:accession) + (1 | age) + (1 | sandpaper))

# FINAL
# Nonsignificant: accession, sandpaper
vic.gmrt.modreduced <- lmer(data=vic.comb.germ, ppd50~lifespan + lifespan/species.short
                            + (1 | species.short:accession) + (1 | age))

rand(model)
rand(vic.gmrt.modreduced)
anova(vic.gmrt.modreduced)
emmeans(object=vic.gmrt.modreduced, specs = pairwise ~ lifespan) # gives you the p value
vic.gmrt.sum <- emmeans(object=vic.gmrt.modreduced, specs = c("lifespan")) # gives you a workable structure for a dataframe
vic.gmrt.df <- as.data.frame(vic.gmrt.sum)

# GERMINATION PROPORTION
# subset
lath.germ.prop <- merge.acc.all[merge.acc.all$genus=="Lathyrus",]
phas.germ.prop <- merge.acc.all[merge.acc.all$genus=="Phaseolus",]
vic.germ.prop <- merge.acc.all[merge.acc.all$genus=="Vicia",]

# Lathyrus
# ORIGINAL
model <- lmer(data=lath.germ.prop, FINAL.germ.prop~lifespan + lifespan/species.short
              + (1 | age) + (1 | sandpaper))

# FINAL 
# Nonsignificant: sandpaper. Can't have accession, since data was collected at the accession level.
lath.gmprop.modreduced <- model <- lmer(data=lath.germ.prop, FINAL.germ.prop~lifespan + lifespan/species.short
                                        + (1 | age))

rand(model)
rand(lath.gmprop.modreduced) 
anova(lath.gmprop.modreduced)
emmeans(object=lath.gmprop.modreduced, specs = pairwise ~ lifespan) # gives you the p value
lath.gmprop.sum <- emmeans(object=lath.gmprop.modreduced, specs = c("lifespan")) # gives you a workable structure for a dataframe
lath.gmprop.df <- as.data.frame(lath.gmprop.sum)

# Phaseolus
# ORIGINAL
model <- lmer(data=phas.germ.prop, FINAL.germ.prop~lifespan + lifespan/species.short
              + (1 | age) + (1 | sandpaper))

# FINAL 
# Nonsignificant: age, sandpaper.
phas.gmprop.modreduced <- model <- lmer(data=phas.germ.prop, FINAL.germ.prop~lifespan + lifespan/species.short
                                        + (1 | age))

rand(model)
rand(phas.gmprop.modreduced) 
anova(phas.gmprop.modreduced)
emmeans(object=phas.gmprop.modreduced, specs = pairwise ~ lifespan) # gives you the p value
phas.gmprop.sum <- emmeans(object=phas.gmprop.modreduced, specs = c("lifespan")) # gives you a workable structure for a dataframe
phas.gmprop.df <- as.data.frame(phas.gmprop.sum)

# Vicia
# ORIGINAL
model <- lmer(data=vic.germ.prop, FINAL.germ.prop~lifespan + lifespan/species.short
              + (1 | age) + (1 | sandpaper))

# FINAL 
# Nonsignificant: age, sandpaper.
vic.gmprop.modreduced <- model <- lmer(data=vic.germ.prop, FINAL.germ.prop~lifespan + lifespan/species.short
                                       + (1 | age))

rand(model)
rand(vic.gmprop.modreduced) 
anova(vic.gmprop.modreduced)
emmeans(object=vic.gmprop.modreduced, specs = pairwise ~ lifespan) # gives you the p value
vic.gmprop.sum <- emmeans(object=vic.gmprop.modreduced, specs = c("lifespan")) # gives you a workable structure for a dataframe
vic.gmprop.df <- as.data.frame(vic.gmprop.sum)

# Old PCA plots ####
# group by genera
ag <- fviz_pca_ind(a.pca,
                   geom.ind = "point", # show points only (nbut not "text")
                   col.ind = merge.acc.all.d$genus, # color by groups
                   palette = c("midnightblue", "grey50", "dodgerblue"), 
                   addEllipses = T, label = "var", mean.point=F,
                   col.var = "black", repel = T,
                   legend.title = "Genus", title="Genus",
                   xlab="PC1 (41.2%)", ylab="PC2 (24.4%)")

ag <- ag + scale_shape_manual(values=c(16,15,17))+ 
  scale_x_continuous(limits=c(-7.5, 9), breaks=seq(-6, 9, 3), expand=c(0, 0))+
  scale_y_continuous(limits=c(-7.5, 6), breaks=seq(-6, 6, 3), expand=c(0, 0))+ 
  theme_classic()
ag

# group by lifespan
al <- fviz_pca_ind(a.pca,
                   geom.ind = "point", # show points only (nbut not "text")
                   col.ind = merge.acc.all.d$lifespan, # color by groups
                   palette = c("firebrick", "dodgerblue3"),
                   addEllipses = T, label = "var", mean.point=F,
                   col.var = "black", repel = T,
                   legend.title="Lifespan", title="Lifespan",
                   xlab="PC1 (41.2%)", ylab="PC2 (24.4%)")

al <- al + scale_shape_manual(values=c(16,15,17))+ 
  scale_x_continuous(limits=c(-7.5, 9), breaks=seq(-6, 9, 3), expand=c(0, 0))+
  scale_y_continuous(limits=c(-7.5, 6), breaks=seq(-6, 6, 3), expand=c(0, 0))+ 
  theme_classic()
al

# Variables plot from factoextra
var.pca <- fviz_pca_var(a.pca, col.var = "cos2",
                        gradient.cols = c("blue", "orange", "red"), 
                        repel = T, title="Variables",
                        xlab="PC1 (41.2%)", ylab="PC2 (24.4%)")
var.pca

# Additional PCs (colored by genera)
# PC1 vs. PC3
pc13 <- fviz_pca_biplot(a.pca, geom.ind = "point", 
                        col.ind = merge.acc.all.d$genus, axes=c(1,3),
                        palette = c("midnightblue", "grey50", "dodgerblue"), 
                        addEllipses = T, label = "var", mean.point=F,
                        col.var = "black", repel = T,
                        legend.title = "Genus",
                        xlab="PC1 (45.0%)", ylab="PC3 (10.6%)")
pc13 <- pc13 + scale_shape_manual(values=c(16,15,17)) + scale_x_continuous(limits=c(-6,10), breaks=c(-6,-3,0,3,6,9)) +
  scale_y_continuous(limits=c(-4,4), breaks=c(-4,-2,0,2,4)) + theme_classic()
pc13

pc14 <- fviz_pca_biplot(a.pca, geom.ind = "point", 
                        col.ind = merge.acc.all.d$genus, axes=c(1,4),
                        palette = c("midnightblue", "grey50", "dodgerblue"), 
                        addEllipses = T, label = "var", mean.point=F,
                        col.var = "black", repel = T,
                        legend.title = "Genus",
                        xlab="PC1 (45.0%)", ylab="PC4 (7.6%)")
pc14 <- pc14 + scale_shape_manual(values=c(16,15,17)) + scale_x_continuous(limits=c(-6,10), breaks=c(-6,-3,0,3,6,9)) +
  scale_y_continuous(limits=c(-4,4), breaks=c(-4,-2,0,2,4)) + theme_classic()
pc14

pc23 <- fviz_pca_biplot(a.pca, geom.ind = "point", 
                        col.ind = merge.acc.all.d$genus, axes=c(2,3),
                        palette = c("midnightblue", "grey50", "dodgerblue"), 
                        addEllipses = T, label = "var", mean.point=F,
                        col.var = "black", repel = T,
                        legend.title = "Genus",
                        xlab="PC2 (24.4%)", ylab="PC3 (10.6%)")
pc23 <- pc23 + scale_shape_manual(values=c(16,15,17)) + scale_x_continuous(limits=c(-6.5,6), breaks=c(-6,-3,0,3,6)) +
  scale_y_continuous(limits=c(-4,4), breaks=c(-4,-2,0,2,4)) + theme_classic()
pc23

pc24 <- fviz_pca_biplot(a.pca, geom.ind = "point", 
                        col.ind = merge.acc.all.d$genus, axes=c(2,4),
                        palette = c("midnightblue", "grey50", "dodgerblue"), 
                        addEllipses = T, label = "var", mean.point=F,
                        col.var = "black", repel = T,
                        legend.title = "Genus",
                        xlab="PC2 (24.4%)", ylab="PC4 (7.6%)")
pc24 <- pc24 + scale_shape_manual(values=c(16,15,17)) + scale_x_continuous(limits=c(-6.5,6), breaks=c(-6,-3,0,3,6)) +
  scale_y_continuous(limits=c(-4,4), breaks=c(-4,-2,0,2,4)) + theme_classic()
pc24

pc34 <- fviz_pca_biplot(a.pca, geom.ind = "point", 
                        col.ind = merge.acc.all.d$genus, axes=c(3,4),
                        palette = c("midnightblue", "grey50", "dodgerblue"), 
                        addEllipses = T, label = "var", mean.point=F,
                        col.var = "black", repel = T,
                        legend.title = "Genus",
                        xlab="PC3 (10.6%)", ylab="PC4 (7.6%)")
pc34 <- pc34 + scale_shape_manual(values=c(16,15,17)) + scale_x_continuous(limits=c(-4,4), breaks=c(-4,-2,0,2,4)) +
  scale_y_continuous(limits=c(-4,4), breaks=c(-4,-2,0,2,4)) + theme_classic()
pc34

lath.seed.pca <- fviz_pca_biplot(ls.pca, geom.ind = "point", pointshape=16,
                                 col.ind = lath.seed$lifespan,
                                 palette = c("firebrick", "dodgerblue3"),
                                 addEllipses = T, label = "var", mean.point=F,
                                 col.var = "black", repel = T,
                                 legend.title = "Lifespan", title="Lathyrus seed PCA",
                                 xlab="PC1 (69.0%)", ylab="PC2 (23.2%)")
lath.seed.pca

# Variables plot
ls.var.pca <- fviz_pca_var(ls.pca, col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                           repel = TRUE,
                           xlab="PC1 (69.0%)", ylab="PC2 (23.2%)")
ls.var.pca

phas.seed.pca <- fviz_pca_biplot(ps.pca, geom.ind = "point", pointshape=15,
                                 col.ind = phas.seed$lifespan,
                                 palette = c("firebrick", "dodgerblue3"),
                                 addEllipses = T, label = "var", mean.point=F,
                                 col.var = "black", repel = T,
                                 legend.title = "Lifespan", title="Phaseolus seed PCA",
                                 xlab="PC1 (67.1%)", ylab="PC2 (23.3%)")
phas.seed.pca

ps.pca.species <- fviz_pca_biplot(ps.pca, geom.ind = "point",
                                  col.ind = phas.seed$species.short,
                                  addEllipses = T, label = "var", mean.point=F,
                                  col.var = "black", repel = T,
                                  legend.title = "Lifespan", title="Phaseolus seed PCA",
                                  xlab="PC1 (67.1%)", ylab="PC2 (23.3%)")
ps.pca.species

# Variables plot
ps.var.pca <- fviz_pca_var(ps.pca, col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                           repel = TRUE, title="Phaseolus seed variables",
                           xlab="PC1 (67.1%)", ylab="PC2 (23.3%%)")
ps.var.pca

vic.seed.pca <-fviz_pca_biplot(vs.pca, geom.ind = "point", pointshape=17,
                               col.ind = vicia.seed$lifespan,
                               palette = c("firebrick", "dodgerblue3"),
                               addEllipses = T, label = "var", mean.point=F,
                               col.var = "black", repel = T,
                               legend.title = "Lifespan", title="Vicia seed PCA",
                               xlab="PC1 (69.9%)", ylab="PC2 (22.4%)")
vic.seed.pca

# Variables plot
vs.var.pca <- fviz_pca_var(ps.pca, col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                           repel = TRUE, title="Vicia seed variables",
                           xlab="PC1 (69.9%)", ylab="PC2 (22.4%%)")
vs.var.pca

lath.veg.pca <- fviz_pca_biplot(lv.pca, geom.ind = "point", pointshape=16,
                                col.ind = lath.veg$lifespan,
                                palette = c("firebrick", "dodgerblue3"),
                                addEllipses = T, label = "var", mean.point=F,
                                col.var = "black", repel = T,
                                legend.title = "Lifespan", title="Lathyrus vegetative PCA",
                                xlab="PC1 (50.9%)", ylab="PC2 (22.4%)")
lath.veg.pca

# Variables plot
lv.var.pca <- fviz_pca_var(lv.pca, col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                           repel = TRUE, title="Lathyrus vegetative variables",
                           xlab="PC1 (38.4%)", ylab="PC2 (33.0%)")
lv.var.pca

phas.veg.pca <- fviz_pca_biplot(pv.pca, geom.ind = "point", pointshape=15,
                                col.ind = phas.veg$lifespan,
                                palette = c("firebrick", "dodgerblue3"),
                                addEllipses = T, label = "var", mean.point=F,
                                col.var = "black", repel = T,
                                legend.title = "Lifespan", title="Phaseolus vegetative PCA",
                                xlab="PC1 (50.1%)", ylab="PC2 (37.2%)")
phas.veg.pca

# Variables plot
pv.var.pca <- fviz_pca_var(pv.pca, col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                           repel = TRUE, title="Phaseolus vegetative variables",
                           xlab="PC1 (50.1%)", ylab="PC2 (37.2%)")
pv.var.pca

vic.veg.pca <- fviz_pca_biplot(vv.pca, geom.ind = "point", pointshape=17,
                               col.ind = vicia.veg$lifespan,
                               palette = c("firebrick", "dodgerblue3"),
                               addEllipses = T, label = "var", mean.point=F,
                               col.var = "black", repel = T,
                               legend.title = "Lifespan", title="Vicia vegetative PCA",
                               xlab="PC1 (54.8%)", ylab="PC2 (30.1%)")
vic.veg.pca

# Variables plot
vv.var.pca <- fviz_pca_var(vv.pca, col.var = "cos2",
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                           repel = TRUE, title="Vicia vegetative variables",
                           xlab="PC1 (54.8%)", ylab="PC2 (30.1%)")
vv.var.pca

# TOTAL DATA

lath.tot <- merge.acc.all %>% filter(genus=="Lathyrus")
lt <- lath.tot %>% select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, nodes.21.mn, nodes.35.mn, nodes.perday.mn)
estim_ncpPCA(lt, scale=T) # estimates # of dimenstions needed in the next imputation step (ncp)
lt1 <- imputePCA(as.data.frame(lt), ncp=3, scale=T) # imputes missing values (NAs)
lt.pca <- prcomp(lt1$completeObs, scale=T)
get_eigenvalue(lt.pca) # Eigenvalues for PCs
fviz_eig(lt.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(lt.pca, col.var = "black")
fviz_pca_biplot(lt.pca, geom.ind = "point", pointshape=16,
                col.ind = lath.tot$lifespan,
                palette = c("firebrick", "dodgerblue3"),
                addEllipses = T, label = "var", mean.point=F,
                col.var = "black", repel = T,
                legend.title = "Lifespan", title="Lathyrus total PCA",
                xlab="PC1 (42.4%)", ylab="PC2 (21.4%)")

phas.tot <- merge.acc.all %>% filter(genus=="Phaseolus")
pt <- phas.tot %>% select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, nodes.21.mn, nodes.35.mn, nodes.perday.mn)
pt.pca <- prcomp(pt, scale=T)
get_eigenvalue(pt.pca) # Eigenvalues for PCs
fviz_eig(pt.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(pt.pca, col.var = "black")
fviz_pca_biplot(pt.pca, geom.ind = "point", pointshape=16,
                col.ind = phas.tot$lifespan,
                palette = c("firebrick", "dodgerblue3"),
                addEllipses = T, label = "var", mean.point=F,
                col.var = "black", repel = T,
                legend.title = "Lifespan", title="Phaseolus total PCA",
                xlab="PC1 (40.2%)", ylab="PC2 (25.7%)")

vicia.tot <- merge.acc.all %>% filter(genus=="Vicia")
vt <- vicia.tot %>% select(seed.wt.avg, seed.length.mn, seed.width.mn, seed.perim.mn, seed.area.mn, seed.circ.mn, seed.round.mn, ppd50.mn, FINAL.germ.prop, height.21.mn, height.35.mn, height.perday.mn, nodes.21.mn, nodes.35.mn, nodes.perday.mn)
estim_ncpPCA(vt, scale=T) # estimates # of dimenstions needed in the next imputation step (ncp)
vt1 <- imputePCA(as.data.frame(vt), ncp=2, scale=T) # imputes missing values (NAs)
vt.pca <- prcomp(vt1$completeObs, scale=T)
get_eigenvalue(vt.pca) # Eigenvalues for PCs
fviz_eig(vt.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(vt.pca, col.var = "black")
fviz_pca_biplot(vt.pca, geom.ind = "point", pointshape=16,
                col.ind = vicia.tot$lifespan,
                palette = c("firebrick", "dodgerblue3"),
                addEllipses = T, label = "var", mean.point=F,
                col.var = "black", repel = T,
                legend.title = "Lifespan", title="Vicia total PCA",
                xlab="PC1 (55.6%)", ylab="PC2 (24.2%)")
# PCA - Seeds ####

# Seed PCA (by individual seeds)
# only the active quantitative variables will be used initially
s <- comb.seed.clean %>% dplyr::select(seed.circ, seed.roundness, seed.length, seed.width, seed.perim, seed.area)
# using prcomp (SVD). PCA() [FactoMineR] is another SVD method but returns slightly different scores
# prcomp is preferred since it gives most results PCA does as well as PC loadings (FactoMineR does not).
# though other things like eigenvalues seem to be the same for both
s.pca <- prcomp(s, scale=T)
get_eigenvalue(s.pca) # Eigenvalues for PCs
fviz_eig(s.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(s.pca, col.var = "black") # correlation circle
fviz_cos2(s.pca, choice = "var", axes = 1:2) # Cos2, quality of representation
fviz_contrib(s.pca, choice = "var", axes = 1, top = 10) # Contributions of variables to PC1
fviz_contrib(s.pca, choice = "var", axes = 2, top = 10) # Contributions of variables to PC2
ind <- get_pca_ind(s.pca) # gives you coordinates, cos2, & contribution values for all individuals
# note that ind$coord is has the same values as s.pca$x
s.pca$rotation # PC loadings from each variable
sqrt(1/nrow(s.pca$rotation)) # cutoff for 'important' loadings (s.pca$rotation refers to the 6 columns of PCs)
# the sum squares of all loadings for an individual PC must sum to one, 
# this is the calculation for what the loadings would be if all variables contributed equally to that PC
# any variable that has a larger loading than this value contributes more than one 
# variable's worth of information and would be regarded as an important contributor to that PC

# average PC scores by accession (seed)
s.accession <- comb.seed.clean$accession
s.pca.bind <- cbind.data.frame(s.accession, s.pca$x) # combine accession names with PC scores for individuals (ensure the dataframe is in the same order)
seed.acc.pca <- s.pca.bind %>%
  group_by(s.accession) %>%
  summarise(
    sPC1 = mean(PC1, na.rm=T),
    sPC2 = mean(PC2, na.rm=T),
    sPC3 = mean(PC3, na.rm=T),
  ) %>% ungroup () # ungroup allows you to still mutate / summarize grouping variables of seed.acc.pca
View(seed.acc.pca)
names(seed.acc.pca)[1] <- "accession" # change colname to accession for merging

# Individuals PCA figure
# group by genera
fviz_pca_ind(s.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = comb.seed.clean$genus, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = T, # Concentration ellipses
             mean.point=F,
             legend.title = "Genera"
)

# group by lifespan
fviz_pca_ind(s.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = comb.seed.clean$lifespan, # color by groups
             addEllipses = T, # Concentration ellipses
             mean.point=F,
             legend.title = "Lifespan"
)

# Biplot PCA: genus
p <- fviz_pca_biplot(s.pca, geom.ind = "point", 
                     col.ind = comb.seed.clean$genus, 
                     palette = c("midnightblue", "grey50", "dodgerblue"), 
                     addEllipses = T, label = "var", mean.point=F,
                     col.var = "black", repel = T,
                     legend.title = "Genus", title="Seed PCA (by individual)",
                     xlab="PC1 (70.0%)", ylab="PC2 (23.7%)")
p + ggplot2::scale_shape_manual(values=c(16,15,17))

# Biplot PCA: lifespan
p <- fviz_pca_biplot(s.pca, geom.ind = "point", pointshape=16,
                     col.ind = comb.seed.clean$lifespan,
                     palette = c("firebrick", "dodgerblue3"),
                     addEllipses = T, label = "var", mean.point=F,
                     col.var = "black", repel = T,
                     legend.title = "Lifespan", title="Seed PCA (by individual)",
                     xlab="PC1 (70.0%)", ylab="PC2 (23.7%)")
# PCA - vegetative growth ####

p.raw <- comb.plant.clean.new %>% dplyr::select(height.21.new, height.35.new, height.perday.21.35, nodes.21.new, nodes.35.new, nodes.perday.21.35)
p <- p.raw[rowSums(is.na(p.raw)) != ncol(p.raw), ] # remove rows where all data is NA
# p has missing values, so an iterative PCA method is used to impute NAs
estim_ncpPCA(p, scale=T, method.cv="Kfold") # estimates # of dimenstions needed in the next imputation step (ncp)
p1 <- imputePCA(as.data.frame(p), ncp=5, scale=T) # new df with imputed values in place of NAs
p.pca <- prcomp(p1$completeObs, scale=T) # prcomp has also be used instead of PCA() (Josse & Husson 2016)
get_eigenvalue(p.pca) # Eigenvalues for PCs
fviz_eig(p.pca, addlabels = TRUE) # Scree plot
fviz_pca_var(p.pca, col.var = "black") # correlation circle
fviz_cos2(p.pca, choice = "var", axes = 1:2) # Cos2, quality of representation
fviz_contrib(p.pca, choice = "var", axes = 1, top = 10) # Contributions of variables to PC1
fviz_contrib(p.pca, choice = "var", axes = 2, top = 10) # Contributions of variables to PC2
p.pca$rotation # PC loadings from each variable
sqrt(1/nrow(p.pca$rotation))
# check imputation
write.xlsx(p1, file = "C:/Users/SterlingH/Desktop/R_exports/veg_impute.xlsx")
write.xlsx(comb.plant.clean.new, file = "C:/Users/SterlingH/Desktop/R_exports/veg_compare.xlsx")

# average PC scores by accession (vegetative)
p.accession <- comb.plant.clean.new$accession
p.pca.bind <- cbind.data.frame(p.accession, p.pca$x) # combine accession names with PC scores for individuals (ensure the dataframe is in the same order)
veg.acc.pca <- p.pca.bind %>%
  group_by(p.accession) %>%
  summarise(
    vPC1 = mean(PC1, na.rm=T),
    vPC2 = mean(PC2, na.rm=T),
    vPC3 = mean(PC3, na.rm=T),
    vPC4 = mean(PC4, na.rm=T),
  ) %>% ungroup () # ungroup allows you to still mutate / summarize grouping variables of seed.acc.pca
View(veg.acc.pca)
names(veg.acc.pca)[1] <- "accession" # change colname to accession for merging

# group by genera
p <- fviz_pca_ind(p.pca,
                  geom.ind = "point", # show points only (but not "text")
                  col.ind = comb.plant.clean.new$genus, # color by groups
                  palette = c("midnightblue", "grey50", "dodgerblue"),
                  addEllipses = T, # Concentration ellipses
                  mean.point=F,
                  legend.title = "Genera")

p + ggplot2::scale_shape_manual(values=c(16,15,17))

# group by lifespan
fviz_pca_ind(p.pca,
             geom.ind = "point", pointshape=16, # show points only (nbut not "text")
             col.ind = comb.plant.clean.new$lifespan, # color by groups
             palette = c("firebrick", "dodgerblue3"),
             addEllipses = T, # Concentration ellipses
             mean.point=F,
             legend.title = "Lifespan")

# Biplot PCA: genus
p <- fviz_pca_biplot(p.pca, geom.ind = "point", 
                     col.ind = comb.plant.clean.new$genus, 
                     palette = c("midnightblue", "grey50", "dodgerblue"), 
                     addEllipses = T, ellipse.level=0.95, label = "var", mean.point=F,
                     col.var = "black", repel = T,
                     legend.title = "Genus", title="Vegetative PCA (by individual)",
                     xlab="PC1 (63.3%)", ylab="PC2 (20.5%)")
p + ggplot2::scale_shape_manual(values=c(16,15,17))

# Biplot PCA: lifespan
fviz_pca_biplot(p.pca, geom.ind = "point", pointshape=16,
                col.ind = comb.plant.clean.new$lifespan,
                palette = c("firebrick", "dodgerblue3"),
                addEllipses = T, label = "var", mean.point=F,
                col.var = "black", repel = T,
                legend.title = "Lifespan", title="Vegetative PCA (by individual)",
                xlab="PC1 (63.3%)", ylab="PC2 (20.5%)")
