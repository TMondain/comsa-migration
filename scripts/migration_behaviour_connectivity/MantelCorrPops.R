rm(list=ls())

library(tidyverse)
# library()
library("raster")
library("sp")
library("rgdal")
library("geosphere")
library("rgeos")
library("ecodist")
library("readr")

#load functions
source("scripts/custom_functions.R")

# load migration summary
gm <- read_csv("data/movement_data/GLS_mig_R.csv")
head(gm)

# label the individual populations from Scotland (not used in paper)
gm2 <- gm %>% mutate(loc2 = ifelse(loc == "Sedbergh", loc,
                                   ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                                          ifelse(loc == "Scotland" & indiv == "Bird1", "Scotland Spey.",
                                                 ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Spey.",loc))))) 

# Calculate distances between birds for England and Scotland
gm_MigCon <- subset(gm2, loc != "Senegal") %>% 
  group_by(Location = loc) %>% 
  summarise(# number tagged
    Number_tagged = length(brd_lon),
    # Mean breeding longitude & lat
    Tagging_long = mean(brd_lon),
    Tagging_lat = mean(brd_lat),  
    # Max distance between breeding sites
    Tagging_max_dist = (max(lower(distm(cbind(brd_lon, brd_lat), fun = distHaversine))))/1000,
    # Mean distance between breeding sites
    Tagging_median_dist = (meanDist(brd_lon, brd_lat))/1000, 
    # Mean winter longitude & lat
    Mean_final_long = mean(win_lon), 
    Mean_final_lat = mean(win_lat),
    # Max distance between winter sites
    Final_max_dist = (max(lower(distm(cbind(win_lon, win_lat), fun = distHaversine))))/1000,
    # Mean distance between winter sites
    Final_median_dist = (meanDist(win_lon, win_lat))/1000, 
    # Distance between mean breding and winter site
    Migration_dist = (distHaversine(cbind(mean(brd_lon), mean(brd_lat)), cbind(mean(win_lon), mean(win_lat))))/1000) 

# Calculate distances between birds for Senegal only
gm_MigCon_sen <- subset(gm2, loc2 == "Senegal") %>% 
  group_by(Location = loc2) %>% 
  summarise(# number tagged
    Number_tagged = length(brd_lon),
    # Mean wintering longitude & lat
    Tagging_long = mean(win_lon),
    Tagging_lat = mean(win_lat),  
    # Max distance between wintering sites
    Tagging_max_dist = (max(lower(distm(cbind(win_lon, win_lat), fun = distHaversine))))/1000,
    # Mean distance between wintering sites
    Tagging_median_dist = (meanDist(win_lon, win_lat))/1000, 
    # Mean winter longitude & lat
    Mean_final_long = mean(brd_lon), 
    Mean_final_lat = mean(brd_lat),
   # Max distance between breeding sites
    Final_max_dist = (max(lower(distm(cbind(brd_lon, brd_lat), fun = distHaversine))))/1000,
    # Mean distance between breeding sites
    Final_median_dist = (meanDist(brd_lon, brd_lat))/1000, 
    # Distance between mean wintering and breeding site
    Migration_dist = (distHaversine(cbind(mean(win_lon), mean(win_lat)), cbind(mean(brd_lon), mean(brd_lat))))/1000) 

mc <- rbind(gm_MigCon,gm_MigCon_sen)
mc

write.csv(mc, file = "outputs/Migratory_mixing_pops.csv")


### mantel correlation test
# do birds from same location mix at end destination?
data.frame(Tagging_Location = "Breeding", 
      subset(gm2, loc != "Senegal") %>% summarise('mantelR' = mantelCor(brd_lon, brd_lat, win_lon, win_lat),
                                                  'mantelP' = mantelCorP(brd_lon, brd_lat, win_lon, win_lat)))

# mantel test for Senegal out of interest (not used in paper)
gm2s <- subset(gm2, loc == "Senegal")

mantel(lower(distm(cbind(gm2s$brd_lon, gm2s$brd_lat), fun = distHaversine)) ~
         lower(distm(cbind(gm2s$win_lon, gm2s$win_lat), fun = distHaversine)))
