rm(list = ls())


## go through each one in turn to remove superfluous packages!
library(tidyverse)
library(lubridate) # dates
library(gtools) # for mixedsort() function
library(sp) # for spatialpoints()

library(RNCEP)
library(rWind)
library(terra)

source('scripts/custom_functions.R')
source("scripts/wind_analyses/calculate_costs/FunctionsSimulateSpringAut_Simple_after_viva.R")
# source("scripts/wind_analyses/calculate_costs/FixedWind_dl_2.r")



###################################################################
#### calculate cost of individuals and paired simulated tracks ####
###################################################################


##### load in individuals

# Data frame has breeding locations 
GLS_mig <- read_csv("data/movement_data/GLS_mig_R.csv")

# wintering grounds averaged for birds that showed weird mid-winter movements
ms <- read_csv("data/movement_data/Schedule_AllIndivs_MeanWint.csv")

ms$mig[ms$loc == "Senegal" & ms$lat < 25] <- "Winter"

# positions of all individuals
mp <- read_csv("data/movement_data/Positions_AllIndivs.csv")

# tidying
mp <- mp %>% mutate(mig = gsub("Autumn", replacement = "autumn", x = mig),
                    mig = gsub("Spring", replacement = "spring", x = mig)) %>% 
  # subset((jd <256 | jd > 276) & (jd < 69 | jd > 89)) %>% # subset equinox
  na.omit %>% 
  group_by(indiv, mig) %>% 
  mutate(lat = c(smooth(lat, twiceit = T)),
         lon = c(smooth(lon, twiceit = T)))

# make sure bird 7 (Scotland) all set to autumn because failed.
mp$mig[mp$indiv=='Bird7'] <- "autumn"

# define the altitudes to download
pressure_levels <- c(1000, 925, 850, 700) # sea level, 779, 1502 and 3130 m a.s.l



######################################################
####   Step 1: wind rasters for each individual   ####
######################################################

########## get weather data for each population separately
## store as different objects

# names of birds
inds <- unique(mp$indiv)

#output rasters
w_aut_ras <- list()
w_spr_ras <- list()


#output flow dispersion files
fd_aut_out <- list()
fd_spr_out <- list()

# define the grid to download
lat_range <- c(-10, 72.5)      #latitude range
lon_range <- c(-30, 40)     #longitude range

for(w in 1:length(inds)) {
  print(paste(inds[w], w, sep = " "))
  
  # get dates for wind sequence
  d_m <- subset(ms, indiv == inds[w])
  head(d_m)
  
  #### autumn
  # get breeding departure
  dep_br <- d_m[1,4]
  
  # get winter arrival
  arr_w <- na.omit(d_m$Arrival[d_m$mig=="Winter"])[1]
  
  #### spring
  # winter departure
  dep_w <- na.omit(d_m$Departure[d_m$mig=="Winter"])
  dep_w <- dep_w[length(dep_w)]
  
  # Spring arrival
  arr_br <- d_m[dim(d_m)[1],3]
  
  
  # for senegal and scotland use different sequence because of dates formatted differently
  if(unique(d_m$loc) == "Senegal" | unique(d_m$loc) == "Scotland") {
    dt_a <- seq(ymd_hms(dep_br), ymd_hms(arr_w)+days(5), by = "1 day") # autumn
    
    if((length(dep_w)>0)) {
      dt_s <- seq(ymd_hms(dep_w)-days(5), ymd_hms(arr_br), by = "1 day") # spring
      
    } 
    
  }
  
  if(unique(d_m$loc) == "Sedbergh") {
    # create date sequence for wind download
    dt_a <- seq(dmy_hms(dep_br), dmy_hms(arr_w)+days(5), by = "1 day") # autumn
    dt_s <- seq(dmy_hms(dep_w)-days(5), dmy_hms(arr_br), by = "1 day") # spring
    
    # edits for certain individuals after manual checking
    if(unique(d_m$indiv) == "ET") {
      dt_s <- seq(dmy_hms(dep_w)-days(35), dmy_hms(arr_br), by = "1 day") # spring
    }
    
    if(unique(d_m$indiv) == "HE") {
      dt_s <- seq(dmy_hms(dep_w)-days(30), dmy_hms(arr_br), by = "1 day") # spring
    }
  }
  
  
  # if bird didn't return back to the breeding grounds - tag failed!
  # need to include some dummy data to keep code running
  if(length(dep_w)==0){
    print(paste("no spring migration - dummy sequence", inds[w]))
    
    if(unique(d_m$loc) == "Senegal" | unique(d_m$loc) == "Scotland") {
      dt_s <- seq(ymd_hms("2014-04-21 00:16:30"), ymd_hms("2014-04-22 00:16:30"), by = "1 day")
    }
    
    if(unique(d_m$loc) == "Sedbergh") {
      # create date sequence for wind download
      dt_s <- seq(ymd_hms("2014-04-21 00:16:30"), ymd_hms("2014-04-22 00:16:30"), by = "1 day")
    }
  }
  
  # new code to download wind for multiple altitudes
  
  #### workflow
  
  #define the time sequences, dates, months, years
  aut_dates <- get_dates(dt_a)
  spr_dates <- get_dates(dt_s)
  
  # get average wind for autumn
  rast_mean_aut <- lapply(pressure_levels, FUN = function(x) {
    
    wnd_dl <- download_wind(pressure_level = x,
                            date_sequence = aut_dates$date_seq,
                            months_minmax = aut_dates$month_range,
                            years_minmax = aut_dates$year_range,
                            lat_minmax = lat_range,
                            lon_minmax = lon_range,
                            aggreg_function = "mean",
                            return_raster = TRUE)
    
    # # plot to check
    # plot(wnd_dl)
    
    # stack the rasters in right order
    wind_terr <- c(wnd_dl$wind_dir,wnd_dl$wind_speed)
    
    # disaggregate to 0.5 resolution 2.5 to 0.5 = 5
    wnd_fine <- disagg(wind_terr, fact = 5, method = 'bilinear')
    
    # convert to raster
    wnd <- raster::stack(wnd_fine)
    names(wnd) <- c("direction", "speed") # rename for use with rWind
    
    # back to rWind workflow to get transition probabilities between all cells
    flow_disp <- flow.dispersion(wnd, 
                                 type="active",
                                 output="transitionLayer")
    flow_disp <- gdistance::geoCorrection(flow_disp, type="r", multpl=FALSE, scl=TRUE)
    
    return(list(wind_rast = wnd,
                flow_dispersion = flow_disp))
    
  })
  
  # get average wind for spring
  rast_mean_spr <- lapply(pressure_levels, FUN = function(x) {
    
    wnd_dl <- download_wind(pressure_level = x,
                            date_sequence = spr_dates$date_seq,
                            months_minmax = spr_dates$month_range,
                            years_minmax = spr_dates$year_range,
                            lat_minmax = lat_range,
                            lon_minmax = lon_range,
                            aggreg_function = "mean",
                            return_raster = TRUE)
    
    # # plot to check
    # plot(wnd_dl)
    
    # stack the rasters in right order
    wind_terr <- c(wnd_dl$wind_dir,wnd_dl$wind_speed)
    
    # disaggregate to 0.5 resolution 2.5 to 0.5 = 5
    wnd_fine <- disagg(wind_terr, fact = 5, method = 'bilinear')
    
    # convert to raster
    wnd <- raster::stack(wnd_fine)
    names(wnd) <- c("direction", "speed") # rename for use with rWind
    
    # back to rWind workflow to get transition probabilities between all cells
    flow_disp <- flow.dispersion(wnd, 
                                 type="active",
                                 output="transitionLayer")
    flow_disp <- gdistance::geoCorrection(flow_disp, type="r", multpl=FALSE, scl=TRUE)
    
    return(list(wind_rast = wnd, 
                flow_dispersion = flow_disp))
    
  })
  
  names(rast_mean_aut) <- pressure_levels
  names(rast_mean_spr) <- pressure_levels
  
  # save wind raster
  w_aut_ras[[w]] <- lapply(1:length(rast_mean_aut), function(x) rast_mean_aut[[x]]$wind_rast)
  w_spr_ras[[w]] <- lapply(1:length(rast_mean_spr), function(x) rast_mean_spr[[x]]$wind_rast)
  
  names(w_aut_ras[[w]]) <- pressure_levels
  names(w_spr_ras[[w]]) <- pressure_levels
  
  # save transition layer
  fd_aut_out[[w]] <- lapply(1:length(rast_mean_aut), function(x) rast_mean_aut[[x]]$flow_dispersion)
  fd_spr_out[[w]] <- lapply(1:length(rast_mean_spr), function(x) rast_mean_spr[[x]]$flow_dispersion)
  
  names(fd_aut_out[[w]]) <- pressure_levels
  names(fd_spr_out[[w]]) <- pressure_levels
  
}

# name the raster layers names of the individuals
names(w_aut_ras) <- inds
names(w_spr_ras) <- inds

names(fd_aut_out) <- inds
names(fd_spr_out) <- inds

## output rasters
dir.create("data/wind_analyses/wind_cost_rasters/", recursive = TRUE)

# saveRDS(w_aut_ras, file = "data/wind_analyses/wind_cost_rasters/wind_raster_autumn.rds")
# saveRDS(w_spr_ras, file = "data/wind_analyses/wind_cost_rasters/wind_raster_spring.rds")
w_aut_ras <- readRDS("data/wind_analyses/wind_cost_rasters/wind_raster_autumn.rds")
w_spr_ras <- readRDS("data/wind_analyses/wind_cost_rasters/wind_raster_spring.rds")


## output flow dispersion files
# saveRDS(fd_aut_out, file = "data/wind_analyses/wind_cost_rasters/flow_dispersion_autumn.rds")
# saveRDS(fd_spr_out, file = "data/wind_analyses/wind_cost_rasters/flow_dispersion_spring.rds")
fd_aut_out <- readRDS("data/wind_analyses/wind_cost_rasters/flow_dispersion_autumn.rds")
fd_spr_out <- readRDS("data/wind_analyses/wind_cost_rasters/flow_dispersion_spring.rds")

## In the calculation of costs, need to ignore birds in spring that didn't return to 
## the breeding grounds
## Also, their wind data / rasters are not accurate. 


#############################################################
####     Step 2: Simulate n birds for each real bird     ####
#############################################################


# Simulated birds travel between the breeding and wintering ground
# of each real bird

inds <- unique(ms$indiv)

out_si <- list()

for(s in 1:length(inds)) {
  print(inds[s])
  
  d_m <- subset(ms, indiv == inds[s])
  # head(d_m)
  
  mig <- c("autumn", "spring")
  
  for(m in mig) {
    
    n_i <- gsub(pattern = "_", replacement = "", x = inds[s])
    
    
    if(m == "autumn"){
      
      # start loc
      id_s <- subset(GLS_mig, indiv == n_i)
      start <- c(id_s$brd_lon, id_s$brd_lat)
      
      # end loc
      id_e <- subset(GLS_mig, indiv == n_i)
      end <- c(id_e$win_lon, id_e$win_lat)
      
      if(is.na(end[1]) | is.na(end[2])) {
        print(paste("using predefined wintering ground", inds[s]))
        end <- c(-16, 10)
      }
      
      s_aut <- simp_sim(start = start, end = end, print_out = F, n = 100, move = "south")
      
    }
    
    if(m == "spring"){
      
      # start loc
      id_s <- subset(GLS_mig, indiv == n_i)
      start <- c(id_s$win_lon, id_s$win_lat)
      
      if(is.na(start[1]) | is.na(start[2])) {
        print(paste("using predefined wintering ground", inds[s]))
        start <- c(-16, 10)
      }
      
      # end loc
      id_e <- subset(GLS_mig, indiv == n_i)
      end <- c(id_e$brd_lon, id_e$brd_lat)
      
      s_spr <- simp_sim(start = start, end = end, print_out = F, n = 100, move = "north")
      
    }
    
  }
  
  
  s_aut <- data.frame(s_aut, mig = "autumn", loc = unique(d_m$loc), r_id = inds[s])
  s_spr <- data.frame(s_spr, mig = "spring", loc = unique(d_m$loc), r_id = inds[s])
  
  sim_df <- rbind(s_aut, s_spr)
  out_si[[s]] <- sim_df
  
}


sim_t <- do.call("rbind", out_si)
head(sim_t)

# ## output simulated bird tracks
# saveRDS(sim_t, file = "data/wind_analyses/simulated_birds_n100.rds")
sim_t <- readRDS("data/wind_analyses/simulated_birds_n100.rds")


##########################################################
####     Step 3: Get cost of real bird migrations     ####
##########################################################

# this section gets the cost of the real bird migrations and the
# cost associated with geolocation error

# wintering grounds averaged for birds that showed weird mid-winter movements
ms <- read_csv("data/movement_data/Schedule_AllIndivs_RawWint.csv")

# positions of all individuals
mp <- read_csv("data/movement_data/Positions_AllIndivs.csv")

# check the individuals
unique(mp$indiv)

# clean up some minor errors in the datasets
mp <- mp %>% mutate(mig = gsub(pattern = "Spring", replacement = "spring", x = mig),
                    mig = gsub(pattern = "Autumn", replacement = "autumn", x = mig))

mp$mig[mp$indiv=='Bird7'] <- "autumn"

# name the senegalese wintering grounds
ms$mig[ms$loc == "Senegal" & ms$lat < 25] <- "Winter" 

# tidy up positions dataset
# remove equinox
mp <- mp %>% 
  na.omit %>% 
  group_by(indiv) %>% 
  mutate(lat = c(smooth(lat, twiceit = T)),
         lon = c(smooth(lon, twiceit = T)))


# # This code makes a sequence for each senegalese bird from the breeding grounds
# # to the start of the track and vice versa
# # because the tracks don't go all the way to the breeding grounds
# # doesn't have any effect on final results. 
#
# head(GLS_mig)
# gm <- subset(GLS_mig, tag_loc == "sen")
# 
# ## autumn
# ## KL
# 
# ## breeding site
# brd <- gm[gm$indiv == "KL",][,6:7]
# 
# # start of autumn track
# au <- head(mp[mp$indiv == "KL",], 1)[,6:7]
# 
# a_lat <- seq(brd$breed_lat, au$lat, by = -0.5)
# a_lon <- seq(brd$breed_lon, au$lon, length = length(a_lat))
# 
# aut_tr_KL <- cbind(lon = a_lon, lat = a_lat)
# 
# rbind(aut_tr_KL,crds)
# 
# # end of spring track
# sp <- tail(mp[mp$indiv == "KL",], 1)[,6:7]
# 
# s_lat <- seq(sp$lat, brd$breed_lat, by = 0.5)
# s_lon <- seq(sp$lon, brd$breed_lon, length = length(s_lat))
# 
# spr_tr_KL <- cbind(lon = s_lon, lat = s_lat)
# 
# rbind(crds, spr_tr_KL)


all_out <- list()
all_out_sim <- list()
pressure_level_out <- list()

i_r <- unique(ms$indiv)
mig_pos <- list()

for(i in 1:length(i_r)){
  
  print(i_r[i])
  
  ############################################
  ####    Cost of real bird migrations    ####
  ############################################
  
  # subset main position df by indiv
  ind_pos <- subset(mp, indiv == i_r[i]) %>% mutate(tm = ymd(tm))
  
  # subset main schedule df by indiv
  ind_sched <- suppressWarnings(subset(ms, indiv == i_r[i]) %>% mutate(Arrival = dmy_hms(Arrival),
                                                                       Departure = dmy_hms(Departure)))
  
  # senegal and scotland have different date format 
  if(unique(ind_sched$loc == "Senegal") | unique(ind_sched$loc) == "Scotland") {
    ind_sched <- subset(ms, indiv == i_r[i]) %>% mutate(Arrival = ymd_hms(Arrival),
                                                        Departure = ymd_hms(Departure))
  }
  
  # extract the wintering grounds of the individual
  ad_win <- na.omit(ind_sched[ind_sched$mig == "Winter",])
  ad_win <- ad_win[3:4]
  
  # bird2 and bird 7 are different, have NAs
  if(i_r[i] == "Bird2" | i_r[i] == "Bird7") {
    ad_win <- (ind_sched[ind_sched$mig == "Winter",])
    ad_win <- ad_win[3:4]
  }
  
  # get only the positions not in winter
  # add a few extra days to avoid missing the start
  ind_pos_nwint <- subset(ind_pos, ymd(tm) < (ymd_hms(ad_win$Arrival[1]) + days(5)) | 
                            ymd(tm) > (ymd_hms(ad_win$Departure[length(ad_win$Departure)]) - days(5)))
  
  # Because senegalese birds tagged in senegal have to change it from | to &
  if(unique(ind_pos$loc) == "Senegal"){
    ind_pos_nwint <- subset(ind_pos, ymd(tm) < (ymd_hms(ad_win$Arrival[1]) + days(5)) & 
                              ymd(tm) > (ymd_hms(ad_win$Departure[length(ad_win$Departure)]) - days(5)))
  }
  
  
  if(i_r[i] == "AC") {
    ind_pos_nwint <- subset(ind_pos, ymd(tm) < (ymd_hms(ad_win$Arrival[1]) + days(5)) | 
                              ymd(tm) > (ymd_hms(ad_win$Departure[length(ad_win$Departure)]) - days(30)))
  }
  
  if(i_r[i] == "ET") {
    ind_pos_nwint <- subset(ind_pos, ymd(tm) < (ymd_hms(ad_win$Arrival[1]) + days(5)) | 
                              ymd(tm) > (ymd_hms(ad_win$Departure[length(ad_win$Departure)]) - days(35)))
  }
  
  if(i_r[i] == "HE") {
    ind_pos_nwint <- subset(ind_pos, ymd(tm) < (ymd_hms(ad_win$Arrival[1]) + days(5)) | 
                              ymd(tm) > (ymd_hms(ad_win$Departure[length(ad_win$Departure)]) - days(30)))
  }
  
  mig_pos[[i]] <- ind_pos_nwint 
  
  # get the corresponding raster for each real individual
  # check that it's the right individual
  if(i_r[i]!=names(fd_aut_out)[i]) stop('!! raster order does not match individual order')
  
  fd_aut <- fd_aut_out[[i]]
  fd_spr <- fd_spr_out[[i]]
  
  
  cost_out <- list()
  cost_out_sim <- list()
  
  # pressure levlels data frame to store everything relevant to real birds
  pressure_level_df <- data.frame()
  
  # loop through migrations
  migs <- unique(ind_pos_nwint$mig)
  
  for(m in 1:length(migs)) {
    print(migs[m])
    
    # get the real bird
    i1 <- subset(ind_pos_nwint, mig == migs[m])
    crds <- na.omit(i1[,6:7])
    # str(crds)
    
    
    # for real bird tracks
    cst_ind_aut <- list()
    cst_ind_spr <- list()
    
    
    ######
    # get the cost between coord x and x - 1 
    for(x in 2:dim(crds)[1]) {
      
      
      if(crds[x,2] > 70){
        print(paste("latitude north of raster extext", i_r[i], migs[m], mround(crds[x,2],1)))
        
        next
      }
      
      if(crds[x-1,2] > 70){
        print(paste("latitude north of raster extext", i_r[i], migs[m], mround(crds[x-1,2],1)))
        
        next
      }
      
      
      if(crds[x,2] < (-10)) {
        print(paste("latitude south of raster extext", i_r[i], migs[m], mround(crds[x,2],1)))
        
        next
      }
      
      if(crds[x-1,2] < (-10)) {
        print(paste("latitude south of raster extext", i_r[i], migs[m], mround(crds[x-1,2],1)))
        
        next
      }
      
      if(crds[x,1]<=(-30)) {
        print(paste("longitude west of raster extent", i_r[i], migs[m], mround(crds[x,1],1)))
        
        next
      }
      
      if(crds[x-1,1]<=(-30)) {
        print(paste("longitude west of raster extent", i_r[i], migs[m], mround(crds[x-1,1],1)))
        
        next
      }
      
      ## calculate cost of relocations  
      if(migs[m] == "autumn") {
        
        ### real birds
        # calculate dispersal cost for autumn relocation
        # loop through all the different altitudes and pick the minimum cost
        dispersal_costs_aut <- 
          sapply(1:length(fd_aut), FUN = function(i){ 
            gdistance::costDistance(fd_aut[[i]], 
                                    fromCoords = SpatialPoints(crds[x-1,]), 
                                    toCoords = SpatialPoints(crds[x,]))
          })
        
        # get the relocation with the smallest cost
        cd_aut <- min(dispersal_costs_aut)
        
        pressure_level_df <- rbind(pressure_level_df,
                                   data.frame(indiv = i_r[i],
                                              migration = migs[m],
                                              crd_pos = x,
                                              lon = crds[x,1],
                                              lat = crds[x,2],
                                              pressure_index = which.min(dispersal_costs_aut),
                                              pressure_level = pressure_levels[which.min(dispersal_costs_aut)],
                                              cost = cd_aut,
                                              all_costs = paste(round(dispersal_costs_aut, 1), collapse = '_')))
        
        cst_ind_aut[[x]] <- cd_aut
        
      }
      
      if(migs[m] == "spring") {
        
        dispersal_costs_spr <- 
          sapply(1:length(fd_spr), FUN = function(i){ 
            gdistance::costDistance(fd_spr[[i]], 
                                    fromCoords = SpatialPoints(crds[x-1,]), 
                                    toCoords = SpatialPoints(crds[x,]))
          })
        
        # get the relocation with the smallest cost
        cd_spr <- min(dispersal_costs_spr)
        
        pressure_level_df <- rbind(pressure_level_df,
                                   data.frame(indiv = i_r[i],
                                              migration = migs[m],
                                              crd_pos = x,
                                              lon = crds[x,1],
                                              lat = crds[x,2],
                                              pressure_index = which.min(dispersal_costs_spr),
                                              pressure_level = pressure_levels[which.min(dispersal_costs_spr)],
                                              cost = cd_spr,
                                              all_costs = paste(round(dispersal_costs_spr, 1), collapse = '_')))
        
        cst_ind_spr[[x]] <- cd_spr
      }
      
      
    }
    
    if(migs[m] == "autumn") {
      t_aut <- do.call("c",cst_ind_aut)
      t_aut <- t_aut[is.finite(t_aut)]
      df <- data.frame(cost = sum(na.omit(t_aut)), indiv = i_r[i], 
                       loc = unique(i1$loc), mig = migs[m], c_ind = sum(na.omit(t_aut))/(dim(crds)[1]))
      
    }
    
    if(migs[m] == "spring") {
      t_spr <- do.call("c",cst_ind_spr)
      t_spr <- t_spr[is.finite(t_spr)]
      df <- data.frame(cost = sum(na.omit(t_spr)), indiv = i_r[i], 
                       loc = unique(i1$loc), mig = migs[m], c_ind = sum(na.omit(t_spr))/(dim(crds)[1]))
      
    }
    
    ## output of real bird summed costs across all relocations for one migration
    cost_out[[m]] <- df
    
    #################################
    ####    Geolocation error    ####
    #################################
    
    ##############
    ####    Generate simulated 'real' birds dataframe
    
    # generate 100 individuals for each real bird
    # tracks are the same length as the real bird tracks
    # with lon and lat sampled from a normal distribution to 
    # simulate geolocation error.
    # Error is ~ 250km for lat and ~100 ish for longitude
    # 1 degree =~111km, so sample with 1 stdev for lat and 0.5 stdevs for longitude
    # This is generated for each real bird and each migration
    # then bind the true coordinates (from the real bird) to the simulated data
    
    
    # subset the outputted pressure data from above to the right migration 
    # and use the stored coordinates to generate the simulated tracks.
    # need to do this because some of the true coordinates fall outside
    # of the area of interest
    pressure_mig_subset <- pressure_level_df[pressure_level_df$migration == migs[m],]
    
    # simulate tracks
    lon_sim <- data.frame(replicate(100,rnorm(length(pressure_mig_subset$lon), 
                                              pressure_mig_subset$lon, 
                                              0.5))) %>% 
      rownames_to_column()
    
    lat_sim <- data.frame(replicate(100,rnorm(length(pressure_mig_subset$lat), 
                                              pressure_mig_subset$lat, 
                                              1))) %>% 
      rownames_to_column()
    
    # pivot to 'long' format
    lons_t <- lon_sim %>% pivot_longer(cols = 2:101, values_to = "lon") %>% arrange(name)
    lats_t <- lat_sim %>% pivot_longer(cols = 2:101, values_to = "lat") %>% arrange(name)
    
    # create the simulated tracks dataframe and smooth
    locs_sim <- cbind(lons_t, 
                      lats_t[,3],
                      true_lon = pressure_mig_subset$lon,
                      true_lat = pressure_mig_subset$lat,
                      alt_ind = c(pressure_mig_subset$pressure_index),
                      alt = c(pressure_mig_subset$pressure_level)) %>% 
      group_by(name) %>% 
      mutate(lon = c(smooth(lon, twiceit = TRUE)),
             lat = c(smooth(lat, twiceit = TRUE)))
    head(locs_sim)
    
    # sort by simulated bird name
    sim_inds <- unique(locs_sim$name) %>% mixedsort()
    
    
    ######################################################
    ####    Geolocation error start of calculation    ####
    ######################################################
    
    # output list of 
    cost_individual_sim <- list()
    
    for(s in 1:length(sim_inds)){
      print(paste('geolocation error bird', sim_inds[s]))
      
      cst_track_aut_sim <- list()
      cst_track_spr_sim <- list()
      
      crds_sim <- subset(locs_sim, name == sim_inds[s])
      
      # get the cost between coord x and x - 1 
      for(x_s in 2:dim(crds_sim)[1]) {
        
        
        if(crds_sim[x_s,4] > 70){
          print(paste("latitude north of raster extext", sim_inds[s], migs[m], mround(crds_sim[x_s,4],1)))
          
          next
        }
        
        if(crds_sim[x_s-1,4] > 70){
          print(paste("latitude north of raster extext", sim_inds[s], migs[m], mround(crds_sim[x_s-1,4],1)))
          
          next
        }
        
        
        if(crds_sim[x_s,4] < (-10)) {
          print(paste("latitude south of raster extext", sim_inds[s], migs[m], mround(crds_sim[x_s,4],1)))
          
          next
        }
        
        if(crds_sim[x_s-1,4] < (-10)) {
          print(paste("latitude south of raster extext", sim_inds[s], migs[m], mround(crds_sim[x_s-1,4],1)))
          
          next
        }
        
        if(crds_sim[x_s,3]<=(-30)) {
          print(paste("longitude west of raster extent", sim_inds[s], migs[m], mround(crds_sim[x_s,3],1)))
          
          next
        }
        
        if(crds_sim[x_s-1,3]<=(-30)) {
          print(paste("longitude west of raster extent", sim_inds[i], migs[m], mround(crds_sim[x_s-1,3],1)))
          
          next
        }
        
        
        #################################################
        ##   Cost of the geolocation error simulation  ##
        #################################################
        if(migs[m] == "autumn") {
          
          # use the altitude from the corresponding real bird
          cd_aut_sim <- gdistance::costDistance(fd_aut[[crds_sim$alt_ind[x_s]]], 
                                                SpatialPoints(crds_sim[x_s-1,3:4]), 
                                                SpatialPoints(crds_sim[x_s,3:4]))
          
          cst_track_aut_sim[[x_s]] <- cd_aut_sim
          
        }
        
        if(migs[m] == "spring") {
          
          # use the altitude from the corresponding real bird
          cd_spr_sim <- gdistance::costDistance(fd_spr[[crds_sim$alt_ind[x_s]]], 
                                                SpatialPoints(crds_sim[x_s-1,3:4]), 
                                                SpatialPoints(crds_sim[x_s,3:4]))
          
          cst_track_spr_sim[[x_s]] <- cd_spr_sim
          
          
        }
        
        
      }
      
      if(migs[m] == "autumn") {
        t_aut <- do.call("c",cst_track_aut_sim)
        t_aut <- t_aut[is.finite(t_aut)]
        
        # save the DFs for each individual within each migration
        cost_individual_sim[[s]] <- data.frame(cost = sum(na.omit(t_aut)), 
                                               r_in = i_r[i], 
                                               indiv = sim_inds[s], 
                                               loc = unique(i1$loc), 
                                               mig = migs[m], 
                                               c_ind = sum(na.omit(t_aut))/(dim(crds)[1]))
        
      }
      
      if(migs[m] == "spring") {
        t_spr <- do.call("c",cst_track_spr_sim)
        t_spr <- t_spr[is.finite(t_spr)]
        
        # save the DFs for each individual within each migration
        cost_individual_sim[[s]] <- data.frame(cost = sum(na.omit(t_spr)), 
                                               r_in = i_r[i], 
                                               indiv = sim_inds[s], 
                                               loc = unique(i1$loc), 
                                               mig = migs[m], 
                                               c_ind = sum(na.omit(t_spr))/(dim(crds)[1]))
        
      }
      
    }
    
    # store the migration of the simulated bird
    cost_out_sim[[m]] <- do.call("rbind", cost_individual_sim)
    
  }
  
  all_out[[i]] <- do.call("rbind", cost_out) # real birds
  all_out_sim[[i]] <- do.call("rbind", cost_out_sim) # simulated birds
  
  pressure_level_out[[i]] <- pressure_level_df 
  
}

# Cost of the tracks of real birds
cst_mig <- do.call("rbind", all_out)
head(cst_mig)

## saving and checking files

# write.csv(cst_mig, file = "data/wind_analyses/flight_costs/cost_mig_individuals.csv")
cst_mig <- read.csv("data/wind_analyses/flight_costs/cost_mig_individuals.csv")[,-1]


# cost of the tracks accounting for the geolocation error
cst_mig_sim_rl <- do.call("rbind", all_out_sim)

# write.csv(cst_mig_sim_rl, file = "data/wind_analyses/flight_costs/cost_mig_individuals_simulated_geolocation_error.csv")
cst_mig_sim_rl <- read.csv("data/wind_analyses/flight_costs/cost_mig_individuals_simulated_geolocation_error.csv")
head(cst_mig_sim_rl)


# Pressure df real birds
pressure_df <- do.call("rbind", pressure_level_out)

# write.csv(pressure_df, file = "data/wind_analyses/flight_costs/pressure_levels_per_relocation.csv")
press_per_loc <- read.csv("data/wind_analyses/flight_costs/pressure_levels_per_relocation.csv")
head(press_per_loc)


###########################################################
####     Step 4: Get cost of simulated bird tracks     ####
###########################################################

## This code works, but takes a LONG time (>12hr). I instead used the NERC Jasmin
## super computer to run this. This is contained in the two scripts
## simulated_bird_costs_function.R, and
## submit_simulated_bird_costs_function.R

# each track's cost is simulated relative to the conditions
# the real bird experienced during it migration
head(sim_t)
str(sim_t)
# smooth the tracks of simulated birds to be like
# real birds
sim_t <- sim_t %>% group_by(indiv, mig, r_id) %>% 
  mutate(lon = c(smooth(lon, twiceit = TRUE)),
         lat = c(smooth(lat, twiceit = TRUE)))

r_birds <- unique(sim_t$r_id)

r_sim_out <- list()

# start with getting the rasters from real birds
for(r in 1:length(r_birds)){
  print(r_birds[r])
  
  # get the corresponding raster for each real individual
  fd_aut <- fd_aut_out[[r]]
  fd_spr <- fd_spr_out[[r]]
  
  # get the simulated tracks corresponding to each real individual 
  sim_t_r <- subset(sim_t, r_id == r_birds[r])
  
  cost_out <- list()
  
  # for each simulated individual in the dataset
  for(i in 1:length(unique(sim_t_r$indiv))) {
    print(i)
    
    # get simulated individual 
    sim_i <- subset(sim_t_r, indiv == i)
    
    head(sim_i)
    
    cst_ind_aut <- list()
    cst_ind_spr <- list()
    
    # to get the cost for each migration
    migs <- c("autumn", "spring")
    
    out_m_aut <- list()
    out_m_spr <- list()
    
    for(m in migs){
      print(m)
      sim_i_m <- subset(sim_i, mig == m)
      
      for(x in 2:length(sim_i_m$lat)) {
        
        # print(x)
        
        if(unique(sim_i_m$mig) == "autumn") {
          crds <- sim_i_m[,1:2]
          cd_aut <- gdistance::costDistance(fd_aut[[sample(1:length(fd_aut), 1)]], 
                                            SpatialPoints(crds[x-1,]), 
                                            SpatialPoints(crds[x,]))
          cst_ind_aut[[x]] <- cd_aut
        }
        
        
        if(unique(sim_i_m$mig) == "spring") {
          crds <- na.omit(sim_i_m[,1:2])
          cd_spr <- costDistance(fd_spr[[sample(1:length(fd_spr), 1)]], 
                                 SpatialPoints(crds[x-1,]), 
                                 SpatialPoints(crds[x,]))
          cst_ind_spr[[x]] <- cd_spr
        }
        
        
      }
      
      t_aut <- do.call("c",cst_ind_aut)
      
      t_spr <- do.call("c",cst_ind_spr)
      
      out_m_aut[[m]] <- t_aut
      out_m_spr[[m]] <- t_spr
      
    }
    
    ca <- do.call("c", out_m_aut)
    cs <- do.call("c", out_m_spr)
    
    df <- data.frame(cost_aut = sum(na.omit(ca[is.finite(ca)])), 
                     cost_aut_ind = sum(na.omit(ca[is.finite(ca)]))/length(sim_i_m$lon), 
                     cost_spr = sum(na.omit(cs[is.finite(cs)])), 
                     cost_spr_ind = sum(na.omit(cs[is.finite(cs)]))/length(sim_i_m$lon),
                     sim_indiv = i, r_id = r_birds[r], loc = unique(sim_i_m$loc))
    
    cost_out[[i]] <- df
    
  }
  
  
  c_o <- do.call("rbind", cost_out)
  
  r_sim_out[[r]] <- c_o
  
}

sim_output <- do.call(rbind, r_sim_out)

# write.csv(sim_output, file = "data/wind_analyses/flight_costs/combined_simulation_costs.csv")
s_out <- read.csv("data/wind_analyses/flight_costs/combined_simulation_costs.csv")
head(s_out)

