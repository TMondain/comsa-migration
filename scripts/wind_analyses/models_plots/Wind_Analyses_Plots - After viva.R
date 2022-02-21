rm(list = ls())

##############################################
#####     Wind analyses and plotting     #####
##############################################

library(tidyverse)
library(gridExtra)
library(rWind)
library(rgdal)
library(raster)
library(lme4)
library(cowplot)
library(lubridate)
library(patchwork)

source('scripts/custom_functions.R')

#####     load rasters for plotting     #####
load("data/wind_analyses/wind_cost_rasters/w_aut_ras_final")
load("data/wind_analyses/wind_cost_rasters/w_spr_ras_final")


#####     load simulated bird tracks for plotting     #####
load("data/wind_analyses/simulated_bird_tracks/simulated_birds_100_finalV2")


#####     load simulated bird costs     #####
load("data/wind_analyses/flight_costs/cost_out_sim_100inds_finalV2")
s_out <- do.call("rbind", r_sim_out)
head(s_out)


#####     load real bird costs     #####
cst_mig <- read.csv("data/wind_analyses/flight_costs/cost_mig_individuals_final.csv")[,-1]


#####     load geolocation error simulation     #####
glc_error <- read.csv("data/wind_analyses/flight_costs/cost_mig_individuals_simulated_geolocation_error_final.csv")
glc_err <- glc_error[,c(2:3, 5:7)]
glc_err$type <- "real_gls_err"

colnames(glc_err) <- c("cost", "indiv", "loc", "mig", "c_ind", "type")
head(glc_err)

#####     create combined real and simulated bird df     #####
head(s_out)

s_p <- pivot_longer(s_out, cols = c("cost_aut_ind", "cost_spr_ind"), names_to = "mig", values_to = "cost_ind")

sp_t <- s_p[,c(2,4:7)] %>% mutate(mig = ifelse(mig == "cost_aut_ind", "autumn", "spring"),
                                  type = "sim")
colnames(sp_t) <- c("cost", "indiv", "loc", "mig", "c_ind", "type")
cst_mig$type <- "real"

# combined dataset
com_d <- rbind(sp_t, cst_mig, glc_err)
head(com_d)



#####     Get migration positions real birds     #####
# positions of all individuals
mp <- read_csv("data/movement_data/Positions_AllIndivs.csv")

mp <- mp %>% mutate(mig = gsub(pattern = "Spring", replacement = "spring", x = mig),
                    mig = gsub(pattern = "Autumn", replacement = "autumn", x = mig))


# tidy up positions dataset
# remove equinox
mp <- mp %>% 
  na.omit %>% 
  group_by(indiv) %>% 
  mutate(lat = c(smooth(lat, twiceit = T)),
         lon = c(smooth(lon, twiceit = T))) %>% 
  ungroup()



###############################################################
#####     Plot 1 : Boxplot of real vs simulated costs     #####
###############################################################

head(com_d)

c_ind_pl <- ggplot(data = subset(com_d, type != "real_gls_err"), aes(x = loc, y = c_ind, fill = type)) +
  geom_boxplot() +
  facet_wrap(~mig, labeller=labeller(mig=c(spring='Spring',autumn='Autumn'))) + 
  guides(fill=guide_legend(title="Bird type"),
         colour = 'none') +
  xlab("Tagging location") + ylab("Cost index") + 
  scale_fill_manual(labels = c("Observed", "Simulated"), values = c("#F8766D", "#619CFF")) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        legend.title = element_blank(),
        strip.background = element_rect(fill="lightgrey")) +
  scale_x_discrete(labels = c("Cumbria", "Senegal", "Scotland"))
c_ind_pl

# # save
# ggsave(c_ind_pl, file = "outputs/wind_analyses/Cost_index_BoxplotV2.tiff",
#        device = "tiff", dpi = 600, height = 7, width = 8)


c_ind_pl <- ggplot(data = com_d, aes(x = loc, y = c_ind, fill = type)) +
  geom_boxplot() +
  facet_grid(~mig, labeller=labeller(mig=c(spring='Spring',autumn='Autumn'))) + 
  guides(fill=guide_legend(title="Bird type"),
         colour = 'none') +
  xlab("Tagging location") + ylab("Cost index") + 
  scale_fill_manual(labels = c("Observed", "Geolocation error", "Simulated"), values = c("#F8766D", "#00BA38", "#619CFF")) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        legend.title = element_blank(),
        strip.background = element_rect(fill="lightgrey")) +
  scale_x_discrete(labels = c("Cumbria", "Senegal", "Scotland"))
c_ind_pl

# # save
# ggsave(c_ind_pl, file = "outputs/wind_analyses/Cost_index_BoxplotV2_geoloc_error.tiff",
#        device = "tiff", dpi = 600, height = 7, width = 8)



#####     Haven't done any plots of raw costs yet

###############################################
#####     Plot 2 : Raw costs boxplots     #####
###############################################

head(s_out)
head(cst_mig)
head(glc_error)

ge <- glc_error[, c(3,5,6,2)]
ge$type <- 'Geoloc_err'
colnames(ge) <- c( "indiv", "loc", "mig", "cost", "type")
head(ge)

rs <- pivot_longer(s_out, cols = c(cost_aut, cost_spr), names_to = "mig", values_to = "cost")[,-c(1:3)]
head(rs)

rs <- rs %>% mutate(mig = ifelse(mig == "cost_aut", "autumn", "spring"),
                    type = "Simulated")

colnames(rs) <- c( "indiv", "loc", "mig", "cost", "type")

cst_mig$type <- "Real"

raw_comb <- rbind(rs, cst_mig[,c(2:4, 1, 6)], ge)

raw_c_pl <- ggplot(data = raw_comb, aes(x = loc, y = cost, fill = type)) + geom_boxplot() +
  facet_grid(~mig, labeller=labeller(mig=c(spring='Spring',autumn='Autumn'))) + 
  guides(fill=guide_legend(title="Bird type"),
         colour = 'none') +
  xlab("Tagging location") + ylab("Cost") + 
  # scale_fill_manual(labels = c("Real", "Simulated"), values = c("#F8766D", "#619CFF")) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = c("Cumbria", "Senegal", "Scotland"))
raw_c_pl

ggsave(raw_c_pl, file = "outputs/wind_analyses/Cost_raw_BoxplotV2.tiff", device = "tiff", dpi = 600, height = 8, width = 8)



#####################################################
#####     Plot 3 : Example track and raster     #####
#####################################################

# load world map
world.shp <- readOGR("data/worldmap.geojson", verbose = F)
pal.shp <- world.shp %>% crop(., extent(-25, 10, -5, 60)) %>% fortify


head(mp)
mp <- mp %>% filter((jd > 165 & jd < 240) | (jd > 75 & jd < 145)) %>% 
  subset((jd <250 | jd > 280) & (jd < 65 | jd > 95)) %>%
  group_by(indiv) %>% mutate(lat = ma(lat, n=2),
                             lon = ma(lon, n=2))

ggplot(mp, aes(x = lon, y = lat, colour = mig, group = interaction(indiv, mig))) + geom_path() +
  facet_wrap(~indiv)

unique(mp$indiv)

ras <- as(w_aut_ras[[7]], "SpatialPixelsDataFrame")
ra <- as.data.frame(ras)


ras_s <- as(w_spr_ras[[7]], "SpatialPixelsDataFrame")
rs <- as.data.frame(ras_s)

dk_spr <- ggplot() + geom_tile(data = rs, aes(x= x, y = y, fill = speed)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.2) +
  coord_map("mercator", xlim = c(-30, 8), ylim = c(0, 60)) +
  geom_path(data = subset(mp, indiv == "DK" & mig == "spring"), aes(x = lon, y = lat, colour = mig, group = mig)) +
  xlab("") + ylab("") +
  theme_bw()+
  guides(colour = FALSE) +
  theme(text = element_text(size = 15)) + 
  scale_fill_continuous(name = "Wind speed")

dk_aut <- ggplot() + geom_tile(data = ra, aes(x= x, y = y, fill = speed)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.2) +
  coord_map("mercator", xlim = c(-30, 8), ylim = c(0, 60)) +
  geom_path(data = subset(mp, indiv == "DK" & mig == "autumn"), aes(x = lon, y = lat, colour = mig, group = mig)) +
  xlab("") + ylab("") +
  theme_bw() +
  guides(colour = FALSE) +
  theme(text = element_text(size = 15)) + 
  scale_fill_continuous(name = "Wind speed")

ex_pl_ras <- dk_aut|dk_spr + plot_layout(guides = 'collect')
ex_pl_ras

getwd()

ggsave(ex_pl_ras, file = "outputs/Example_Plot_raster_DK2.tiff", device = "tiff", dpi = 600, width = 8, height = 5)


rs_d <- rs %>% mutate(lat = mround(y, 3),
                      lon = mround(x, 3)) %>% 
  group_by(lat,lon) %>% 
  summarise(dir = median(direction),
            spd = median(speed))

head(rs_d)  

w_spring <- ggplot() + geom_tile(data = rs, aes(x= x, y = y, fill = speed)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.2) +
  coord_map("mercator", xlim = c(-25, 8), ylim = c(0, 60))+
  geom_spoke(data = rs_d, aes(x = lon, y = lat, angle = dir, radius = scales::rescale(spd, c(0.5, 2))), colour = "white", 
             arrow = arrow(length = unit(.05, 'inches'))) +
  geom_path(data = subset(mp, indiv == "DK" & mig == "spring"), aes(x = lon, y = lat, colour = mig, group = mig)) + 
  xlab("") + ylab("") +
  theme_bw() +
  guides(colour = FALSE) +
  theme(text = element_text(size = 15)) + 
  scale_fill_continuous(name = "Wind speed")


ra_d <- ra %>% mutate(lat = mround(y, 3),
                      lon = mround(x, 3)) %>% 
  group_by(lat,lon) %>% 
  summarise(dir = median(direction),
            spd = median(speed)) %>% 
  mutate(spd = ifelse(spd>12, 12, spd))

head(ra_d)  

ra$speed <- ifelse(ra$speed>12, 12, ra$speed)

w_aut <- ggplot() + geom_tile(data = ra, aes(x= x, y = y, fill = speed)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.2) +
  coord_map("mercator", xlim = c(-25, 8), ylim = c(0, 60)) +
  geom_spoke(data = ra_d, aes(x = lon, y = lat, angle = dir, radius = scales::rescale(spd, c(0.5, 2))), colour = "white", 
             arrow = arrow(length = unit(.05, 'inches'))) +
  geom_path(data = subset(mp, indiv == "DK" & mig == "autumn"), aes(x = lon, y = lat, colour = mig, group = mig)) + 
  xlab("") + ylab("") +
  theme_bw() +
  guides(colour = FALSE) +
  theme(text = element_text(size = 15),
        legend.position = 'none') + 
  scale_fill_continuous(name = "Wind speed")


# ex_pl_ras_wind <- grid.arrange(w_aut, w_spring, ncol = 2)

ex_pl_ras_wind <- w_aut|w_spring + plot_layout(guides = 'collect')
ex_pl_ras_wind

# ggsave(ex_pl_ras_wind, file = "outputs/Example_Plot_raster_and_Wind_DKV2.tiff", device = "tiff", dpi = 600, width = 8, height = 5)




#####################################################
#####      Plot 4 : Simulated bird tracks       #####
#####################################################


head(sim_t)
pl_t <- sim_t %>% mutate(loc = as.character(loc),
                         loc2 = ifelse(loc == "Sedbergh", "Sedbergh",
                                       ifelse(loc == "Scotland" & (r_id != "Bird1" & r_id != "Bird2"), "Scotland Suth.", 
                                              ifelse(loc == "Scotland" & r_id == "Bird1", "Scotland Spey.",
                                                     ifelse(loc == "Scotland" & r_id == "Bird2", "Scotland Spey.", loc))))) %>% 
  group_by(indiv, mig, loc, r_id) %>% 
  mutate(lon = smooth(lon, twiceit = T),
         lat = smooth(lat, twiceit = T),
         lon_seq = seq(1, length(lon)))

plt_sum <- pl_t %>% 
  group_by(r_id, mig, loc, lon_seq) %>% 
  summarise(lon = mean(lon),
            lat = mean(lat))


ap <- ggplot(data = pl_t, aes(x = lon, y = lat, group = interaction(indiv, mig, loc, r_id), colour = loc2)) + 
  geom_path() +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.2) +
  coord_map("mercator", xlim = c(-30, 35), ylim = c(0, 75)) +
  labs(colour = "Tagging location") +
  ylab("") + xlab("") +
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  scale_colour_discrete(labels = c("Scotland Spey.", "Scotland Suth.", 
                                   "Cumbria", "Senegal"))


sp <- ggplot(data = plt_sum, aes(x = lon, y = lat, group = interaction(mig, loc, r_id), colour = loc2)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.2) +
  coord_map("mercator", xlim = c(-30, 35), ylim = c(0, 75)) + 
  geom_path() +
  labs(colour = "Tagging location") +
  ylab("") + xlab("") +
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  scale_colour_discrete(labels = c("Scotland Spey.", "Scotland Suth.", 
                                   "Cumbria", "Senegal"))

cp <- ap + sp + plot_layout(guides = 'collect')

# cp <- plot_grid(ap, sp, ncol = 2)#, labels = c("a", "b"))
cp

ggsave(cp, file = "outputs/Comb_Sim_Sum_pl2.tiff", device = "tiff", width = 10, height = 5)
