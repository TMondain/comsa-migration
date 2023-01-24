rm(list = ls())

library(mgcv)
library(tidyverse)
library(rgdal)
library(raster)
library(lubridate)
library(adehabitatHR)

source('scripts/custom_functions.R')


#####################################################
#####          PLOTTING MIGRATION MAPS          #####
#####################################################

# load world map
world.shp <- readOGR("data/worldmap.geojson", verbose = F)

# crop to palearctic
pal.shp <- world.shp %>% crop(., extent(-80, 155, -40, 90)) %>% fortify


## load migration summary
mg <- read_csv("data/movement_data/GLS_mig_R.csv")
head(mg)

# tidy dates
g <- mg %>% mutate(brd_dep = ymd_hms(brd_dep),
                   brd_dep_jd = yday(brd_dep),
                   brd_dep_yr = year(brd_dep),
                   win_arr = ymd_hms(win_arr),
                   win_arr_jd = yday(win_arr),
                   win_arr_yr = year(win_arr),
                   win_dep = ymd_hms(win_dep),
                   win_dep_jd = yday(win_dep),
                   win_dep_yr = year(win_dep),
                   brd_arr = ymd_hms(brd_arr),
                   brd_arr_jd = yday(brd_arr),
                   brd_arr_yr = year(brd_arr))


## load Migration schedule
ms <- read_csv("data/movement_data/Schedule_AllIndivs_MeanWint.csv")
head(ms)

# set senegal winter site
ms$mig[ms$loc == "Senegal" & ms$lat < 25] <- "Winter"

# make sure spring/autumn labels are consistent
ms <- ms %>% mutate(mig = ifelse(mig == "Autumn", "autumn", 
                                 ifelse(mig == "Spring", "spring", 
                                        ifelse(mig == "Winter", "winter", mig))))


## load positions of all individuals
mp <- read_csv("data/movement_data/Positions_AllIndivs.csv")
head(mp)

# autumn/spring consistency and remove missing points
mp <- mp %>% mutate(mig = gsub("Autumn", replacement = "autumn", x = mig),
                    mig = gsub("Spring", replacement = "spring", x = mig)) %>% 
  na.omit 

# create initial dataset
mpt2 <- mp %>% 
  # filter((jd > 165 & jd < 240) | (jd > 75 & jd < 145)) %>% ## Remove clearly erroneous locations from wintering grounds in the tracks
  mutate(loc2 = ifelse(loc == "Sedbergh", loc,
                                    ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                                           ifelse(loc == "Scotland" & indiv == "Bird1", "Scotland Spey.",
                                                  ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Spey.",loc)))),
                      mig2 = ifelse(is_wint == "NotWinter", mig, "winter")) %>% 
  filter((jd < 70 | jd > 90) & (jd < 254 | jd > 274) & site > 0) %>% # remove 20 days either side of equinoxes
  dplyr::select(lon, lat, mig2, loc) %>% # subset data
  na.omit #%>% 
# mutate(lat = ma((lat), n = 4), 
#        lon = ma((lon), n = 4)) 

# if want to try and plot senegal in one migration
mpt2$loc[grepl(mpt2$loc, pattern = "Senegal")] <- "Senegal"

# When birds reach certain thresholds assume they have reached breeding locations
mpt2$lat[mpt2$loc == "Sedbergh" & mpt2$lat > 50] <- NA
mpt2$lat[grepl(pattern = "Scotland", mpt2$loc) & mpt2$lat > 52] <- NA
mpt2$lat[mpt2$loc == "Senegal" & mpt2$lat > 57] <- NA

# some weird stopover sites listed as in wintering grounds
mpt2$mig2[mpt2$lat < 20] <- "winter"

# remove breeding sites
mpt2 <- na.omit(mpt2)


#####################################################
#####          KERNEL DENSITY ANALYSES          #####
#####################################################

# Run kernel density analyses for each migration period separately
# 1. spring
# 2. autumn
# 3. winter

#########################
####    1. Spring

k_spr <- subset(mpt2, mig2 == "spring")

# create a grid extent for the kd function to look for densities within
x <- seq(min(k_spr$lon) - 20, max(k_spr$lon) + 20, by=0.5) # resolution is the pixel size you desire 
y <- seq(min(k_spr$lat) - 20, max(k_spr$lat) + 20, by=0.5)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

# create spatial data frame
spr <- SpatialPointsDataFrame(coordinates((cbind(k_spr$lon, k_spr$lat))),
                              data = k_spr, proj4string=CRS("+init=epsg:3395"))

# run KD analysis
kd_spr <- kernelUD(spr[,4], grid = xy,
                   extent = 1, h = "href", hlim = c(0.1, 7))

# Get 75 KD region
c95_spr <- cbind(fortify(getverticeshr(kd_spr, percent = 75)), mig = "spring")

# plot for funsies
ggplot() + 
  geom_polygon(data = c95_spr, aes(x = long, y = lat, group = group, fill = id), alpha = 0.7) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1)  +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(-5, 72)) +
  theme_bw()


#########################
####    2. Autumn

k_aut <- subset(mpt2, mig2 == "autumn")

# create a grid extent for the kd function to look for densities within
x <- seq(min(k_aut$lon) - 20, max(k_aut$lon) + 20, by=0.5) # resolution is the pixel size you desire 
y <- seq(min(k_aut$lat) - 20, max(k_aut$lat) + 20, by=0.5)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

# create spatial data frame
aut <- SpatialPointsDataFrame(coordinates((cbind(k_aut$lon, k_aut$lat))),
                              data = k_aut, proj4string=CRS("+init=epsg:3395"))

# run KD analysis
kd_aut <- kernelUD(aut[,4], grid = xy,
                   extent = 1, h = "href", hlim = c(0.1, 7))

# Get 75 KD region
c95_aut <- cbind(fortify(getverticeshr(kd_aut, percent = 75)), mig = "autumn")

# plot for funsies
ggplot() + 
  geom_polygon(data = c95_aut, aes(x = long, y = lat, group = group, fill = id), alpha = 0.7) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1)  +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(-5, 72)) +
  theme_bw()


#########################
####    3. Winter

k_win <- subset(mpt2, mig2 == "winter" & lat < 35)

# create a grid extent for the kd function to look for densities within
x <- seq(min(k_win$lon) - 20, max(k_win$lon) + 20, by=0.5) # resolution is the pixel size you desire 
y <- seq(min(k_win$lat) - 20, max(k_win$lat) + 20, by=0.5)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

# create spatial data frame
win <- SpatialPointsDataFrame(coordinates((cbind(k_win$lon, k_win$lat))),
                              data = k_win, proj4string=CRS("+init=epsg:3395"))

# run KD analysis
kd_win <- kernelUD(win[,4], grid = xy,
                   extent = 1, h = "href", hlim = c(0.1, 7))

# Get 75 KD region
c95_win <- cbind(fortify(getverticeshr(kd_win, percent = 75)), mig = "winter")

# plot for funsies
ggplot() + 
  geom_polygon(data = c95_win, aes(x = long, y = lat, group = group, fill = id), alpha = 0.7) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1)  +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(-5, 72)) +
  theme_bw()


## Combine them all
kd_all <- rbind(c95_spr, c95_aut, c95_win)
head(kd_all)



###########################################################
#####          PLOTTING MIGRATION MAPS, FIGURE        #####
###########################################################

## set breeding sites data frame for plotting
brd_sites <- data.frame(loc = c("Sedbergh", "Scotland", "Senegal", "Senegal", "Senegal", "Senegal"), 
                        lat = c(54.3, 57.95, c(g$brd_lat[g$loc == "Senegal"])),
                        lon = c(-2.55, -3.9, c(g$brd_lon[g$loc == "Senegal"]))) 


## setting location for scottish birds by breeding location separately
g <- g %>% mutate(loc2 = ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                                ifelse(loc == "Scotland" & (indiv == "Bird1" | indiv == "Bird2"), "Scotland Spey.", 
                                       ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Suth.",loc))))

# manipulate dataset for plotting tracks
pl_c <- mp %>% 
  filter((jd > 165 & jd < 240) | (jd > 75 & jd < 145)) %>% ## Remove clearly erroneous locations from wintering grounds in the tracks
  na.omit %>% 
  mutate(lat = ma((lat), n = 4),
         lon = ma((lon), n = 4)) %>% 
  mutate(loc2 = ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                       ifelse(loc == "Scotland" & (indiv == "Bird1" | indiv == "Bird2"), "Scotland Spey.", 
                              ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Suth.",loc)))) %>% 
  filter((jd < 70 | jd > 90) & (jd < 254 | jd > 274) & site > 0) ## filter out movements and equinox


# When birds reach certain thresholds assume they have reached breeding locations
pl_c$lat[pl_c$loc == "Sedbergh" & pl_c$lat > 52] <- 54.4
pl_c$lon[pl_c$loc == "Sedbergh" & pl_c$lat > 52] <- -2.55
pl_c$lat[pl_c$loc == "Scotland Spey." & pl_c$lat > 56] <- 58.5
pl_c$lon[pl_c$loc == "Scotland Spey." & pl_c$lat > 56] <- -4.3
pl_c$lat[pl_c$loc == "Scotland Suth." & pl_c$lat > 55] <- 57.4
pl_c$lon[pl_c$loc == "Scotland Suth." & pl_c$lat > 55] <- -3.5

comb_mig_pl <- ggplot() + 
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.6) +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(-10, 72)) +
  geom_polygon(data = kd_all, aes(x = long, y = lat, group = group, fill = id), alpha = 0.5) + 
  # geom_path(data = subset(pl_c, lat>10), aes(x = lon, y = lat, group = indiv, colour = loc), alpha = 0.6, size = 0.4) +
  geom_point(data = brd_sites, aes(x = lon, y = lat, colour = loc, shape = "Breeding"), size = 3.5) +
  geom_point(data = g, aes(x = win_lon, y = win_lat, colour = loc, shape = "Wintering"), size = 3.5) +
  facet_grid(~mig, labeller=labeller(mig=c(spring='Spring',autumn='Autumn',winter='Winter'))) + 
  theme_bw() + 
  theme(text = element_text(size = 15), panel.spacing = unit(1.5, "lines")) + 
  guides(fill=guide_legend(title="Tagging location"),
         shape=guide_legend(title = "Site"),
         colour = 'none') +
  xlab("") + ylab("") +
  scale_fill_discrete(labels = c('Scotland', 'Cumbria', 'Senegal'))

comb_mig_pl

# save
ggsave(comb_mig_pl, file = "outputs/movements/Fig_1_Combined_kd_mig_pl_75.tiff", device = "tiff", width = 12, height = 8, dpi = 600)




################################################################
####     Find intersection between different populations    ####
################################################################

### spatial points data frames

ps <- data.frame(stage = "Spring", kerneloverlaphr(kd_spr, percent = 75, method = "HR", conditional = T))
pa <- data.frame(stage = "Autumn", kerneloverlaphr(kd_aut, percent = 75, method = "HR", conditional = T))
pw <- data.frame(stage = "Winter", kerneloverlaphr(kd_win, percent = 75, method = "HR", conditional = T))

p <- rbind(pa,ps,pw)

p <- cbind(location = rep(c("Scotland", "Sedbergh", "Senegal"), 3), p)

rownames(p) <- NULL
p



