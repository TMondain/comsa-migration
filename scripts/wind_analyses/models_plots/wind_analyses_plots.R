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
# library(cowplot)
library(lubridate)
library(patchwork)

library(rnaturalearth)
library(sf)

source('scripts/custom_functions.R')

#####     load world map     #####
# world.shp <- readOGR("data/worldmap.geojson", verbose = F)
# world_coordinates <- map_data("world")
world.shp <- ne_countries(scale = "medium", returnclass = "sf")
plot(st_geometry(world.shp))

#####     load rasters for plotting     #####
w_aut_ras <- readRDS("data/wind_analyses/wind_cost_rasters/wind_raster_autumn.rds")
w_spr_ras <- readRDS("data/wind_analyses/wind_cost_rasters/wind_raster_spring.rds")

#####     load simulated bird tracks for plotting     #####
simulated_bird_tracks <- readRDS("data/wind_analyses/simulated_birds_n100.rds")


#####     load migration costs     #####
com_d <- read.csv('data/wind_analyses/combined_real_simulated_costs.csv') %>% 
  mutate(loc = ifelse(loc == 'Sedbergh', 'Cumbria', loc))


#####     load pressure data     #####
pressure_df <- read.csv('data/wind_analyses/flight_costs/pressure_levels_per_relocation.csv')



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
  mutate(lat = c(smooth(lat, twiceit = TRUE)),
         lon = c(smooth(lon, twiceit = TRUE))) %>% 
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
        strip.background = element_rect(fill="lightgrey"))
c_ind_pl

# # save
# ggsave(c_ind_pl, file = "outputs/wind_analyses/Cost_index_comparison.tiff",
#        device = "tiff", dpi = 600, height = 7, width = 8)


c_ind_gls_pl <- ggplot(data = com_d, aes(x = loc, y = c_ind, fill = type)) +
  geom_boxplot() +
  facet_grid(~mig, labeller=labeller(mig=c(spring='Spring',autumn='Autumn'))) + 
  guides(fill=guide_legend(title="Bird type"),
         colour = 'none') +
  xlab("Tagging location") + ylab("Cost index") + 
  scale_fill_manual(labels = c("Observed", "Geolocation error", "Simulated"), values = c("#F8766D", "#00BA38", "#619CFF")) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        legend.title = element_blank(),
        strip.background = element_rect(fill="lightgrey"))
c_ind_gls_pl

# # save
# ggsave(c_ind_gls_pl, file = "outputs/wind_analyses/cost_index_geoloc_error_comparison.tiff",
#        device = "tiff", dpi = 600, height = 7, width = 8)

# flip the plot around to get observed birds' autumn and spring migration
# next to each other
ggplot(subset(com_d), aes(x = loc, y = c_ind, fill = mig)) +
  geom_boxplot() + 
  guides(fill=guide_legend(title="Migration"),
         colour = 'none') +
  xlab("Tagging location") + ylab("Cost index") +
  facet_wrap(~type,
             labeller=labeller(type=c(real='Observed', 
                                      real_gls_err='Geolocation error', 
                                      sim = 'Simulated'))) + 
  theme_classic()


###############################################
#####     Plot 2 : Raw costs boxplots     #####
###############################################

## Somewhat meaningless because of different numbers of relocations 

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

raw_c_pl <- ggplot(data = com_d, aes(x = loc, y = cost, fill = type)) + geom_boxplot() +
  facet_grid(~mig, labeller=labeller(mig=c(spring='Spring',autumn='Autumn'))) + 
  guides(fill=guide_legend(title="Bird type"),
         colour = 'none') +
  xlab("Tagging location") + ylab("Cost") + 
  # scale_fill_manual(labels = c("Real", "Simulated"), values = c("#F8766D", "#619CFF")) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  scale_x_discrete(labels = c("Cumbria", "Senegal", "Scotland"))
raw_c_pl

# ggsave(raw_c_pl, file = "outputs/wind_analyses/Cost_raw_BoxplotV2.tiff", device = "tiff", dpi = 600, height = 8, width = 8)



#####################################################
#####     Plot 3 : Example track and raster     #####
#####################################################

pal.shp <- sf::st_make_valid(world.shp) %>% 
  st_crop(xmin = -25, xmax = 10, 
          ymin = -5, ymax = 60)# %>% fortify


head(mp)
mp <- mp %>% filter((jd > 165 & jd < 240) | (jd > 75 & jd < 145)) %>% 
  subset((jd <250 | jd > 280) & (jd < 65 | jd > 95)) %>%
  group_by(indiv) %>% mutate(lat = ma(lat, n=2),
                             lon = ma(lon, n=2))

ggplot(mp, aes(x = lon, y = lat, colour = mig, group = interaction(indiv, mig))) + 
  geom_path() +
  facet_wrap(~indiv)

unique(mp$indiv)

ras <- as(w_aut_ras[[7]][[1]], "SpatialPixelsDataFrame")
ra <- as.data.frame(ras) %>% 
  mutate(mig = 'autumn')


ras_s <- as(w_spr_ras[[7]][[1]], "SpatialPixelsDataFrame")
rs <- as.data.frame(ras_s) %>% 
  mutate(mig = 'spring')

dk_spr <- ggplot() + geom_tile(data = rs, aes(x= x, y = y, fill = speed)) +
  geom_sf(data = pal.shp,
          fill = NA, colour = "black", linewidth = 0.4) +
  coord_map("mercator", xlim = c(-25, 8), ylim = c(0, 60)) +
  geom_path(data = subset(mp, indiv == "DK" & mig == "spring"), 
            aes(x = lon, y = lat, colour = mig, group = mig)) +
  xlab("") + ylab("") +
  theme_bw()+
  guides(colour = 'none') +
  theme(text = element_text(size = 15)) + 
  scale_fill_continuous(name = "Wind speed")

dk_aut <- ggplot() + geom_tile(data = ra, aes(x= x, y = y, fill = speed)) +
  geom_sf(data = pal.shp,
          fill = NA, colour = "black", size = 0.4) +
  coord_map("mercator", xlim = c(-25, 8), ylim = c(0, 60)) +
  geom_path(data = subset(mp, indiv == "DK" & mig == "autumn"), aes(x = lon, y = lat, colour = mig, group = mig)) +
  xlab("") + ylab("") +
  theme_bw() +
  guides(colour = FALSE) +
  theme(text = element_text(size = 15)) + 
  scale_fill_continuous(name = "Wind speed")

ex_pl_ras <- dk_aut|dk_spr + plot_layout(guides = 'collect')
ex_pl_ras


# ggsave(ex_pl_ras, file = "outputs/Example_Plot_raster_DK2.tiff", device = "tiff", dpi = 600, width = 8, height = 5)

dk_comb <- rbind(ra, rs)

dk_examp <- ggplot() + 
  geom_tile(data = dk_comb, aes(x= x, y = y, fill = speed)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.4) +
  coord_map("mercator", xlim = c(-25, 8), ylim = c(0, 60)) +
  geom_path(data = subset(mp, indiv == "DK"), aes(x = lon, y = lat, colour = mig, group = mig)) +
  xlab("") + ylab("") +
  theme_bw() +
  scale_colour_manual(values = c('#F8766D','#F8766D')) +
  guides(colour = 'none') +
  theme(text = element_text(size = 15)) + 
  scale_fill_continuous(name = "Wind speed") +
  facet_wrap(~mig, ncol = 2,
             labeller = as_labeller(c(autumn = 'Autumn', spring = 'Spring')))
dk_examp

# ggsave(dk_examp, file = "outputs/Example_Plot_raster_DK2_patch.tiff", device = "tiff", dpi = 600, width = 8, height = 5)



rs_d <- rs %>% mutate(lat = mround(y, 3),
                      lon = mround(x, 3)) %>% 
  group_by(lat,lon) %>% 
  summarise(dir = median(direction),
            spd = median(speed))

head(rs_d)  

w_spring <- ggplot() + geom_tile(data = rs, aes(x= x, y = y, fill = speed)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.4) +
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
               fill = NA, colour = "black", size = 0.4) +
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

# resize palearctic flyway
pal.shp <- world.shp %>% crop(., extent(-80, 155, -40, 90)) %>% fortify


head(simulated_bird_tracks)
pl_t <- simulated_bird_tracks %>% 
  mutate(loc = as.character(loc),
         loc2 = ifelse(loc == "Sedbergh", "Sedbergh",
                       ifelse(loc == "Scotland" & (r_id != "Bird1" & r_id != "Bird2"), "Scotland Suth.", 
                              ifelse(loc == "Scotland" & r_id == "Bird1", "Scotland Spey.",
                                     ifelse(loc == "Scotland" & r_id == "Bird2", "Scotland Spey.", loc))))) %>% 
  group_by(indiv, mig, loc, r_id) %>% 
  mutate(lon = c(smooth(lon, twiceit = T)),
         lat = c(smooth(lat, twiceit = T)),
         lon_seq = seq(1, length(lon)))

plt_sum <- pl_t %>% 
  group_by(r_id, mig, loc, lon_seq) %>% 
  summarise(lon = mean(lon),
            lat = mean(lat))


ap <- ggplot(data = pl_t, aes(x = lon, y = lat, group = interaction(indiv, mig, loc, r_id), colour = loc)) + 
  geom_path() +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.4) +
  coord_map("mercator", xlim = c(-30, 35), ylim = c(0, 75)) +
  labs(colour = "Tagging location") +
  ylab("") + xlab("") +
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  scale_colour_discrete(labels = c("Scotland", 
                                   "Cumbria", "Senegal"))


sp <- ggplot(data = plt_sum, aes(x = lon, y = lat, group = interaction(mig, loc, r_id), colour = loc)) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.4) +
  coord_map("mercator", xlim = c(-30, 35), ylim = c(0, 75)) + 
  geom_path() +
  labs(colour = "Tagging location") +
  ylab("") + xlab("") +
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  scale_colour_discrete(labels = c("Scotland", 
                                   "Cumbria", "Senegal"))

cp <- ap + sp + plot_layout(guides = 'collect')
cp

ggsave(cp, file = "outputs/Comb_Sim_Sum_pl2.tiff", device = "tiff", width = 10, height = 5)
ggsave(sp, file = "outputs/Comb_Sim_average_tracks.tiff", device = "tiff", width = 6, height = 5)




########################################################
#####      Plot 5 : Optimal migration height       #####
########################################################

head(pressure_df)
head(mp) 

# get Europe
euna <- sf::st_make_valid(world.shp) %>% 
  st_crop(xmin = -25, xmax = 40, 
          ymin = -5, ymax = 70)

mpl <- ggplot() +
  geom_sf(data = euna, fill = 'lightgrey') +
  # geom_path(data = mp, 
  #           aes(x = lon, y = lat, colour = mig, group = interaction(indiv, mig)),
  #           alpha = 0.75) +
  ylim(-5,70) +
  theme_bw()
mpl

# altitude by migration with points
ggplot(subset(pressure_df), 
       aes(lat, pressure_index, colour = migration)) +
  geom_point() +
  geom_smooth() +
  ylim(0.5,4.2) +
  theme_classic() +
  facet_wrap(~migration)

# altitude by migration, smooth only
alt_p <- ggplot(pressure_df, 
                aes(lat, pressure_index, colour = migration)) +
  # geom_point() +
  geom_smooth() +
  xlab('Latitude') +
  ylab('Pressure index') +
  ylim(0, 4.2) +
  xlim(-5, 70) +
  theme_bw() +
  coord_flip()
alt_p

mpl + alt_p +
  plot_layout(ncol = 2, widths = c(2,1))


library(gamm4)

m1 <- gamm4(pressure_index ~ s(lat, by = factor(migration)),
            random = ~(1|indiv),
            data = pressure_df)
summary(m1$gam)
m1

ggplot(pressure_df, aes(x = factor(pressure_index), y = lat, fill = migration)) +
  geom_boxplot()

library(viridis)
ggplot(subset(pressure_df, cost >0), aes(lon, lat, colour = log(cost))) +
  geom_point() +
  scale_colour_viridis()
