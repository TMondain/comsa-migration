rm(list = ls())

library(mgcv)
library(tidyverse)
library(rgdal)
library(raster)
library(lubridate)
library(adehabitatHR)


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



mp_t <- mp %>% filter((jd > 165 & jd < 240) | (jd > 75 & jd < 145))

ma <- function(arr, n=15){
  res = arr
  for(i in n:length(arr)){
    res[i] = mean(arr[(i-n+1):i])
  }
  res
}



mp_t <- subset(mp_t, lat > 5) %>% group_by(indiv) %>% mutate(lat = ma((lat), n = 4),
                                                             lon = ma((lon), n = 4))
mp_t$lat[mp_t$loc == "Sedbergh" & mp_t$lat > 52] <- 54.4
mp_t$lon[mp_t$loc == "Sedbergh" & mp_t$lat > 52] <- -2.55
mp_t$lat[mp_t$loc == "Scotland" & mp_t$lat > 57] <- 58
mp_t$lon[mp_t$loc == "Scotland" & mp_t$lat > 57] <- -4



mp_t$lat[mp_t$loc == "Scotland" & mp_t$lat > 20 & mp_t$lon < -15] <- NA

ms$lat[ms$indiv == "KT" & ms$lat > 58.6] <- NA

head(g)

ggplot(data = mp_t, aes(x = lon, y = lat, group = indiv, colour = loc)) + geom_path(alpha = 0.1, size = 1.5) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1) +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(0, 72)) +
  geom_point(data = subset(ms, mig == "spring"),
             aes(x = lon, y = lat, fill = loc, size = stopover, pch = loc), alpha = 0.5) +
  geom_point(data = subset(ms, mig == "autumn"),
             aes(x = lon, y = lat, fill = loc, size = stopover, pch = loc), alpha = 0.5) +
  geom_point(data = g, aes(x = win_lon, y = win_lat, fill = loc), pch = 23, size = 4) +
  geom_point(data = g, aes(x = brd_lon, y = brd_lat, fill = loc, pch = loc), size = 4) +
  facet_wrap(~mig) +
  theme_bw()

library(gganimate)
library(plotly)

eg_pl <- ggplot()  +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.3) +
  coord_map("mercator", xlim = c(-30, 5), ylim = c(0, 60)) + 
  geom_path(data = subset(mp_t, indiv == 'CA')[-c(1:35,160:166),], aes(x = lon, y = lat, group = mig, colour = mig), 
            alpha = 0.8, size = 1, lineend = 'round')+
  theme_bw() +
  theme(legend.position = 'none') +
  xlab("") + ylab("") +
  transition_reveal(jd)


# animate(eg_pl, height = 800, width =400)
# anim_save("Outputs/CAmigration.gif")

# ggsave(eg_pl, file = "Outputs/CAmigration.jpeg", device = 'jpeg',
#        dpi = 300, height = 5, width = 3)


ggplot(data = mp_t, aes(x = lon, y = lat, group = indiv, colour = loc)) + geom_path(alpha = 0.1, size = 1.5) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1) +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(0, 72)) +
  geom_point(data = subset(ms, mig == "spring"),
             aes(x = lon, y = lat, fill = loc, size = stopover, pch = "Stopovers"), alpha = 0.5) +
  geom_point(data = subset(ms, mig == "autumn"),
             aes(x = lon, y = lat, fill = loc, size = stopover, pch = "Stopovers"), alpha = 0.5) +
  geom_point(data = g, aes(x = win_lon, y = win_lat, fill = loc, pch = "Wintering"), size = 4) +
  geom_point(data = g, aes(x = brd_lon, y = brd_lat, fill = loc, pch = "Breeding"), size = 4) +
  facet_wrap(~mig) +
  theme_bw()


# mp_t$lat[mp_t$lat < 15] <- NA

ggplot() + geom_path(data = subset(mp_t, lat>15), aes(x = lon, y = lat, group = indiv, colour = loc), alpha = 0.7, size = 0.2) + 
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1) +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(0, 72)) +
  geom_point(data = subset(ms, mig == "spring"), 
             aes(x = lon, y = lat, fill = loc, size = stopover)) +
  geom_point(data = subset(ms, mig == "autumn"), 
             aes(x = lon, y = lat, fill = loc, size = stopover)) +
  geom_density_2d(data = subset(mp_t, lat < 35), aes(x = lon, y = lat),
                  bins = 4, binwidth = 10000, colour = "red", size = 1, n = 500, h = 20) +#, pch = 17, size = 4) +
  geom_point(data = g, aes(x = brd_lon, y = brd_lat, fill = loc, pch = loc), pch = 18, size = 4) +
  facet_wrap(~mig) +
  xlab("") + ylab("") +
  theme_bw()


ggplot(data = g, aes(x = win_lon, y = win_lat)) + geom_density_2d(bins = 3)



head(mp_t)

mp_t2 <- mp %>% mutate(loc2 = ifelse(loc == "Sedbergh", loc,
                                     ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                                            ifelse(loc == "Scotland" & indiv == "Bird1", "Scotland Spey.",
                                                   ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Spey.",
                                                          ifelse(indiv == "KL", "Senegal_1",
                                                                 ifelse(indiv == "KP", "Senegal_2",
                                                                        ifelse(indiv == "KT", "Senegal_3",
                                                                               ifelse(indiv == "KW", "Senegal_4", "")))))))),
                       mig2 = ifelse(is_wint == "NotWinter", mig, "winter")) %>% 
  filter((jd < 70 | jd > 90) & (jd < 254 | jd > 274) & site > 0)

mp_t2 <- mp %>% mutate(loc2 = ifelse(loc == "Sedbergh", loc,
                                     ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                                            ifelse(loc == "Scotland" & indiv == "Bird1", "Scotland Spey.",
                                                   ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Spey.",loc)))),
                       mig2 = ifelse(is_wint == "NotWinter", mig, "winter")) %>% 
  filter((jd < 70 | jd > 90) & (jd < 254 | jd > 274) & site > 0)


mpt2 <- mp_t2 %>% ungroup %>%  dplyr::select(lon, lat, mig2, loc2) %>% na.omit %>% 
  mutate(lat = ma((lat), n = 4),
         lon = ma((lon), n = 4)) 

# if want to try and plot senegal in one migration
mpt2$loc2[grepl(mpt2$loc2, pattern = "Senegal")] <- "Senegal"

# if want to remove breeding locations
mpt2$lat[mpt2$loc2 == "Sedbergh" & mpt2$lat > 50] <- NA

mpt2$lat[grepl(pattern = "Scotland", mpt2$loc2) & mpt2$lat > 52] <- NA

mpt2$lat[mpt2$loc2 == "Senegal" & mpt2$lat > 57] <- NA

# some weird stopover sites listed as in wintering grounds
mpt2$mig2[mpt2$lat < 20] <- "winter"


mpt2 <- na.omit(mpt2)

spdf <- SpatialPointsDataFrame(coordinates((cbind(mpt2$lon, mpt2$lat))),
                               data = mpt2)

kd <- kernelUD(spdf[,4],
               extent = 1, h = "href", hlim = c(0.1, 7))

c95 <- (fortify(getverticeshr(kd, percent = 75)))
head(c95)

ggplot(data = c95, aes(x = long, y = lat, group = group, fill = id)) + 
  geom_polygon(alpha = 0.2) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1) +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(-5, 72)) 


##################################
##     Split the migrations     ##
##################################

k_spr <- subset(mpt2, mig2 == "spring")

# create a grid extent for the kd function to look for densities within
x <- seq(min(k_spr$lon) - 20, max(k_spr$lon) + 20, by=0.5) # resolution is the pixel size you desire 
y <- seq(min(k_spr$lat) - 20, max(k_spr$lat) + 20, by=0.5)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

spr <- SpatialPointsDataFrame(coordinates((cbind(k_spr$lon, k_spr$lat))),
                              data = k_spr, proj4string=CRS("+init=epsg:3395"))
kd_spr <- kernelUD(spr[,4], grid = xy,
                   extent = 1, h = "href", hlim = c(0.1, 7))

c95_spr <- cbind(fortify(getverticeshr(kd_spr, percent = 75)), mig = "spring")

ggplot() + 
  geom_polygon(data = c95_spr, aes(x = long, y = lat, group = group, fill = id), alpha = 0.2) +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.1)  +
  coord_map("mercator", xlim = c(-30, 30), ylim = c(-5, 72))



k_aut <- subset(mpt2, mig2 == "autumn")


# create a grid extent for the kd function to look for densities within
x <- seq(min(k_aut$lon) - 20, max(k_aut$lon) + 20, by=0.5) # resolution is the pixel size you desire 
y <- seq(min(k_aut$lat) - 20, max(k_aut$lat) + 20, by=0.5)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

aut <- SpatialPointsDataFrame(coordinates((cbind(k_aut$lon, k_aut$lat))),
                              data = k_aut, proj4string=CRS("+init=epsg:3395"))
kd_aut <- kernelUD(aut[,4], grid = xy,
                   extent = 1, h = "href", hlim = c(0.1, 7))

c95_aut <- cbind(fortify(getverticeshr(kd_aut, percent = 75)), mig = "autumn")

ggplot() + 
  geom_polygon(data = c95_aut, aes(x = long, y = lat, group = group, fill = id), alpha = 0.2)




k_win <- subset(mpt2, mig2 == "winter" & lat < 35)

# create a grid extent for the kd function to look for densities within
x <- seq(min(k_win$lon) - 20, max(k_win$lon) + 20, by=0.5) # resolution is the pixel size you desire 
y <- seq(min(k_win$lat) - 20, max(k_win$lat) + 20, by=0.5)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

win <- SpatialPointsDataFrame(coordinates((cbind(k_win$lon, k_win$lat))),
                              data = k_win, proj4string=CRS("+init=epsg:3395"))
# win<-spTransform(win,CRS("+proj=utm +south +zone=55 +datum=WGS84"))
# win 
kd_win <- kernelUD(win[,4], grid = xy,
                   extent = 1, h = "href", hlim = c(0.1, 7))

c95_win <- cbind(fortify(getverticeshr(kd_win, percent = 75)), mig = "winter")

ggplot() + 
  geom_polygon(data = c95_win, aes(x = long, y = lat, group = group, fill = id), alpha = 0.2)


image(kd_aut)
image(kd_spr)
image(kd_win)


kd_all <- rbind(c95_spr, c95_aut, c95_win)
head(kd_all)

head(g)
brd_sites <- data.frame(loc = c("Sedbergh", "Scotland Spey.", "Scotland Suth.", "Senegal", "Senegal", "Senegal", "Senegal"), 
                        lat = c(54.3, 58.5, 57.4, c(g$brd_lat[g$loc == "Senegal"])),
                        lon = c(-2.55, -4.3, -3.5, c(g$brd_lon[g$loc == "Senegal"])))

g <- g %>% mutate(loc2 = ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                                ifelse(loc == "Scotland" & (indiv == "Bird1" | indiv == "Bird2"), "Scotland Spey.", 
                                       ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Suth.",loc))))


getwd()

pl_c <- mp %>% filter((jd > 165 & jd < 240) | (jd > 75 & jd < 145)) %>% na.omit %>% 
  mutate(lat = ma((lat), n = 4),
         lon = ma((lon), n = 4)) %>% 
  mutate(loc = ifelse(loc == "Scotland" & (indiv != "Bird1" & indiv != "Bird2"), "Scotland Suth.", 
                      ifelse(loc == "Scotland" & (indiv == "Bird1" | indiv == "Bird2"), "Scotland Spey.", 
                             ifelse(loc == "Scotland" & indiv == "Bird2", "Scotland Suth.",loc)))) %>% 
  filter((jd < 70 | jd > 90) & (jd < 254 | jd > 274) & site > 0)


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
  geom_point(data = g, aes(x = win_lon, y = win_lat, colour = loc2, shape = "Wintering"), size = 3.5) +
  facet_grid(~mig, labeller=labeller(mig=c(spring='Spring',autumn='Autumn',winter='Winter'))) + 
  theme_bw() + 
  theme(text = element_text(size = 15), panel.spacing = unit(1.5, "lines")) + 
  guides(fill=guide_legend(title="Tagging location"),
         shape=guide_legend(title = "Site"),
         colour = FALSE) +
  xlab("") + ylab("") +
  scale_fill_discrete(labels = c("Scotland Spey.", "Scotland Suth.", 
                                 "Cumbria", "Senegal"))
comb_mig_pl

getwd()
ggsave(comb_mig_pl, file = "Outputs/Combined_kd_mig_pl_75_darker.tiff", device = "tiff", width = 12, height = 8, dpi = 600)

###   Try and animated plot
###   with ggplotly 
eg_pl2 <- ggplot()  +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.3) +
  coord_map("mercator", xlim = c(-30, 5), ylim = c(0, 60)) + 
  geom_point(data = subset(pl_c, mig =='autumn') , aes(x = lon, y = lat, group = interaction(mig, indiv), colour = indiv, frame = jd), 
             alpha = 0.8, size = 1)+
  theme_bw() +
  theme(legend.position = 'none') +
  xlab("") + ylab("") 

fig <- ggplotly(eg_pl2, layer = 'geom_path')

fig

###   with gganimate 

eg_pl3 <- ggplot()  +
  geom_polygon(data = pal.shp, aes(x = long, y = lat, group = group),
               fill = NA, colour = "black", size = 0.3) +
  coord_map("mercator", xlim = c(-30, 25), ylim = c(0, 70)) + 
  geom_path(data = pl_c, aes(x = lon, y = lat, group = interaction(mig, indiv), colour = loc), 
            alpha = 0.3, size = 1)+
  theme_bw() +
  theme(legend.position = 'none') +
  xlab("") + ylab("") +
  transition_reveal(jd) +
  facet_wrap(~mig)

eg_pl3

animate(eg_pl3, height = 900, width =1600)
anim_save("Outputs/AllmigrationByLoc.gif")


################################################################
####     Find intersection between different populations    ####
################################################################

library(raster)
### spatial points data frames

ps <- data.frame(stage = "Spring", kerneloverlaphr(kd_spr, percent = 75, method = "HR", conditional = T))
pa <- data.frame(stage = "Autumn", kerneloverlaphr(kd_aut, percent = 75, method = "HR", conditional = T))
pw <- data.frame(stage = "Winter", kerneloverlaphr(kd_win, percent = 75, method = "HR", conditional = T))

p <- rbind(pa,ps,pw)

p <- cbind(location = rep(c("Scotland Spey.", "Scotland Suth.", "Sedbergh", "Senegal"), 3), p)

rownames(p) <- NULL
p

### to be determined whether this is useful
# it's dependent on the projection which I can't seem to get right...
kernel.area(kd_spr, percent = 95, unin = "km", unout = "km2")
kernel.area(kd_aut, percent = 95, unin = "km", unout = "km2")
kernel.area(kd_win, percent = 95, unin = "km", unout = "km2")

