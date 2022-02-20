
#############################################################
####  Functions Autumn and spring simulating tracks ####
############################################################

##################################################
###       Autumn and Spring random tracks      ###
##################################################

simp_sim <- function(start, end, n, by = 0.5, print_out = T, move = c("north", "south")) {
  
  mround <- function(x,base){
    base*round(x/base)
  }
  
  lon_out <- list()
  
  reps <- list()
  
  if(move == "south") {
    lat <- seq(from = start[2], to = end[2], by = -by)
    lon_shift <- seq(abs(start[1]),abs(end[1]), length = length(lat))
    
    if(start[1]>0){
      lon_shift <- seq((-start[1]),(-end[1]), length = length(lat))
    }
    
  }
  
  if(move == "north") {
    lat <- seq(from = start[2], to = end[2], by = by)
    lon_shift <- seq((start[1]),(end[1]), length = length(lat))
    
  }
  
  
  for(j in 1:n) {
    
    if(print_out == T){
      print(j)
    }
    
    for(i in 1:length(lat)){
      if(move == "south") {
        lon_val <- rnorm(1, mean = -lon_shift[i], sd = 4)
        lon_val
        
        if(lat[i]<end[2]+4 & lat[i]>end[2]+2) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = 3.5)
        }
        
        if(lat[i]<end[2]+2.5 & lat[i]>end[2]+0.5) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = 3)
        }
        
        if(lat[i]<end[2]+1) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = 0.01)
        }
        
        while (abs(lon_val) > abs(lon_shift[i])+12) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = 1)
        }
        
        lon_out[[i]] <- mround(lon_val, 0.5) 
      }
      
      if(move == "north") {
        lon_val <- rnorm(1, mean = lon_shift[i], sd = 4)
        lon_val
        
        if(lat[i]>(end[2]-4) & lat[i]<end[2]-2) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = 3.5)
        }
        
        if(lat[i]>end[2]-2.5 & lat[i]<end[2]-0.5) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = 3)
        }
        
        if(lat[i]>end[2]-0.5) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = 0.01)
        }
        
        while (abs(lon_val) > abs(lon_shift[i])+12) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = 1)
        }
        
        lon_out[[i]] <- mround(lon_val, 0.5) #- lon_shift[i]
      }
      
    }
    
    l <- do.call("rbind", lon_out)
    
    l <- rbind(start[1], l)
    
    reps[[j]] <- l
    
  }
  
  tracks_aut <- do.call("rbind", reps)
  lat <- c(start[2], lat)
  pl_aut <- data.frame(lon = tracks_aut, lat = rep(lat, n), indiv = rep(1:n, each=length(lat)))
  pl_aut
}




library(tidyverse)
library(lubridate)
library(rgdal)
library(raster)
library(gridExtra)
library(rWind)
library(rworldmap)  
library(readr)
library(gdistance)
library(shape)  


mround <- function(x,base){
  base*round(x/base)
}


setwd("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/WindAnalysis/SimpleSimSprAutFunctionRealComp")
setwd("C:/Users/mondainm/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/WindAnalysis/SimpleSimSprAutFunctionRealComp")
source("FunctionsSimulateSpringAut_Simple.R")
source("FixedWind_dl_2.r")

setwd("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/WindAnalysis/SimpleSimSprAutFunctionRealComp/Final")


###################################################################
#### calculate cost of individuals and paired simulated tracks ####
###################################################################


##### load in individuals

# Data frame has breeding locations 
GLS_mig <- read_csv("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/GLS_Summary_Movements_Stopovers/GLS_mig_R.csv")

# wintering grounds averaged for birds that showed weird mid-winter movements
ms <- read_csv("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/Raw_Tracks_Mig_Schedules/All_Combined/Schedule_AllIndivs_MeanWint.csv")

ms$mig[ms$loc == "Senegal" & ms$lat < 25] <- "Winter"

# positions of all individuals
mp <- read_csv("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/Raw_Tracks_Mig_Schedules/All_Combined/Positions_AllIndivs.csv")

mp <- mp %>% mutate(mig = gsub("Autumn", replacement = "autumn", x = mig),
                    mig = gsub("Spring", replacement = "spring", x = mig)) %>% 
  subset((jd <256 | jd > 276) & (jd < 69 | jd > 89)) %>% 
  na.omit %>% 
  group_by(indiv, mig) %>% mutate(lat = smooth(lat, twiceit = T),
                                  lon = smooth(lon, twiceit = T))

mp$mig[mp$indiv=='Bird7'] <- "autumn"

######################################################
####   Step 1: wind rasters for each individual   ####
######################################################


########## get weather data for each population separately
## store as different objects

inds <- unique(mp$indiv)


#output rasters
w_aut_ras <- list()
w_spr_ras <- list()


#output flow dispersion files
fd_aut_out <- list()
fd_spr_out <- list()


for(w in 1:length(inds)) {
  print(paste(inds[w], w, sep = " "))
  
  # get dates for wind sequence
  d_m <- subset(ms, indiv == inds[w])
  head(d_m)
  
  #### autumn
  dep_br <- d_m[1,4]
  
  
  arr_w <- na.omit(d_m$Arrival[d_m$mig=="Winter"])[1]
  
  #### spring
  dep_w <- na.omit(d_m$Departure[d_m$mig=="Winter"])
  dep_w <- dep_w[length(dep_w)]
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
  
  
  # download wind
  ww_aut <- w_t(dt_a, -30, 40, -10, 72.5)
  ww_spr <- w_t(dt_s, -30, 40, -10, 72.5)
  
  # get the mean wind
  w_mean_aut <- wind.mean(ww_aut)
  w_mean_spr <- wind.mean(ww_spr)
  
  # convert to a raster
  r_mean_aut <- wind2raster(w_mean_aut)
  r_mean_spr <- wind2raster(w_mean_spr)
  
  # save raster
  w_aut_ras[[w]] <- r_mean_aut
  w_spr_ras[[w]] <- r_mean_spr
  
  
  fd_aut <- flow.dispersion(r_mean_aut, type="active",
                            output="transitionLayer")
  fd_aut <- geoCorrection(fd_aut, type="r", multpl=FALSE, scl=TRUE)
  
  fd_aut_out[[w]] <- fd_aut
  
  
  fd_spr <- flow.dispersion(r_mean_spr, type="active",
                            output="transitionLayer")
  fd_spr <- geoCorrection(fd_spr, type="r", multpl=FALSE, scl=TRUE)
  
  fd_spr_out[[w]] <- fd_spr
  
  
}

getwd()
#output rasters
# save(w_aut_ras, file = "w_aut_ras_final")
# save(w_spr_ras, file = "w_spr_ras_final")
load("w_aut_ras_final")
load("w_spr_ras_final")


#output flow dispersion files
# save(fd_aut_out, file = "fd_aut_out_final")
# save(fd_spr_out, file = "fd_spr_out_final")
load("fd_aut_out_final")
load("fd_spr_out_final")

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
      
      s_aut <- simp_sim(start = start, end = end, by = dim(mp[mp$indiv == inds[s] & mp$mig == m,])[1], print_out = F, n = 100, move = "south")
      
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
      
      if(dim(mp[mp$indiv == inds[s] & mp$mig == m,])[1]==0){
        next
      }
      
      s_spr <- simp_sim(start = start, end = end, by = dim(mp[mp$indiv == inds[s] & mp$mig == m,])[1], print_out = F, n = 100, move = "north")
      
    }
    
  }
  
  
  s_aut <- data.frame(s_aut, mig = "autumn", loc = unique(d_m$loc), r_id = inds[s])
  s_spr <- data.frame(s_spr, mig = "spring", loc = unique(d_m$loc), r_id = inds[s])
  
  sim_df <- rbind(s_aut, s_spr)
  out_si[[s]] <- sim_df
  
}


sim_t <- do.call("rbind", out_si)
head(sim_t)

# save(sim_t, file = "simulated_birds_100_finalV2")
load("simulated_birds_100_finalV2")


ggplot(data = subset(sim_t, r_id == "KW" & indiv == 1), aes(x = lon, y = lat, group = r_id, colour=mig)) +
  geom_path() 
ggplot(data = sim_t, aes(x = lon, y = lat, group = interaction(indiv, mig), colour=loc)) +
  geom_path() + facet_wrap(~mig)



###########################################################
####     Step 3: Get cost of simulated bird tracks     ####
###########################################################

# each track's cost is simulated relative to the conditions
# the real bird experienced during it migration
head(sim_t)

# smooth the tracks of simulated birds to be like
# real birds
sim_t <- sim_t %>% group_by(indiv, mig, r_id) %>% 
  mutate(lon = smooth(lon, twiceit = T),
         lat = smooth(lat, twiceit = T))

r_birds <- unique(mp$indiv)

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
          cd_aut <- costDistance(fd_aut, SpatialPoints(crds[x-1,]), SpatialPoints(crds[x,]))
          cst_ind_aut[[x]] <- cd_aut
        }
        
        
        if(unique(sim_i_m$mig) == "spring") {
          crds <- na.omit(sim_i_m[,1:2])
          cd_spr <- costDistance(fd_aut, SpatialPoints(crds[x-1,]), SpatialPoints(crds[x,]))
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

# saved at number 18. 
# save(r_sim_out, file = "cost_out_sim_100inds_finalV2")
load("cost_out_sim_100inds_finalV2")


s_out <- do.call("rbind", r_sim_out)

head(s_out)

s<-s_out
colnames(s) <- c("cost_aut_rl", "cost_aut_ind", "cost_spr_rl", "cost_spr_ind", "sim_indiv", "r_id",      "loc")

pivot_longer(s, cols = c("cost_aut_rl", "cost_spr_rl","cost_aut_ind", "cost_spr_ind"), 
             names_to = c("mig", "type"), names_sep = "cost_", values_to = c("cost_rl", "cost_ind"))


# get raw cost into long format 
s_raw_c <- pivot_longer(s_out, cols = c("cost_aut", "cost_spr"), names_to = "mig", values_to = "cost")


# get cost index into long format
s_p <- pivot_longer(s_out, cols = c("cost_aut_ind", "cost_spr_ind"), names_to = "mig", values_to = "cost_ind")

ggplot(data = s_p, aes(cost_ind, fill = mig)) + geom_histogram() +
  facet_wrap(~mig)
ggplot(data = s_p, aes(x = r_id, y = cost_ind, colour = mig)) + geom_boxplot()
ggplot(data = s_p, aes(x = loc, y = cost_ind, colour = mig)) + geom_boxplot()



##########################################################
####     Step 4: Get cost of real bird migrations     ####
##########################################################

# wintering grounds averaged for birds that showed weird mid-winter movements
ms <- read_csv("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/Raw_Tracks_Mig_Schedules/All_Combined/Schedule_AllIndivs_RawWint.csv")

# positions of all individuals
mp <- read_csv("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/Raw_Tracks_Mig_Schedules/All_Combined/Positions_AllIndivs.csv")

# only sedbergh birds
# mp <- subset(mp, loc == "Sedbergh")
mp <- mp %>% mutate(mig = gsub(pattern = "Spring", replacement = "spring", x = mig),
                    mig = gsub(pattern = "Autumn", replacement = "autumn", x = mig))
unique(mp$indiv)

# name the senegalese wintering grounds
ms$mig[ms$loc == "Senegal" & ms$lat < 25] <- "Winter" 

# tidy up positions dataset
# remove equinox
mp <- mp %>% 
  # subset((jd <256 | jd > 276) & (jd < 69 | jd > 89)) %>% # removed because I account for equinoxes lower down in the for loop based on each individual's mig sched  
  na.omit %>% 
  group_by(indiv) %>% mutate(lat = smooth(lat, twiceit = T),
                             lon = smooth(lon, twiceit = T))

mp$mig[mp$indiv=='Bird7'] <- "autumn"


# # need to make a sequence for each senegalese bird from the breeding grounds
# # to the start of the track and vice versa
# # because the tracks don't go all the way to the breeding grounds
# # Might not matter because am dividing by the number of cells the birds cross.... 
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

i_r <- unique(ms$indiv)
mig_pos <- list()

for(i in 1:length(i_r)){
  
  print(i_r[i])
  
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
  
  # bird2 is different, has NAs
  if(i_r[i] == "Bird2") {
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
  fd_aut <- fd_aut_out[[i]]
  fd_spr <- fd_spr_out[[i]]
  
  
  cost_out <- list()
  
  ##########################
  ##  For each migration  ##
  ##########################
  
  migs <- unique(ind_pos_nwint$mig)
  
  for(m in 1:length(migs)) {
    print(migs[m])
    
    i1 <- subset(ind_pos_nwint, mig == migs[m])
    crds <- na.omit(i1[,6:7])
    # str(crds)
    
    cst_ind_aut <- list()
    cst_ind_spr <- list()
    
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
      
      ##############
      ##   COST   ##
      ##############
      if(migs[m] == "autumn") {
        cd_aut <- costDistance(fd_aut, SpatialPoints(crds[x-1,]), SpatialPoints(crds[x,]))
        
        cst_ind_aut[[x]] <- cd_aut
      }
      
      if(migs[m] == "spring") {
        cd_spr <- costDistance(fd_spr, SpatialPoints(crds[x-1,]), SpatialPoints(crds[x,]))
        
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
    
    
    cost_out[[m]] <- df
    
  }
  
  all_out[[i]] <- do.call("rbind", cost_out)
  
}

cst_mig <- do.call("rbind", all_out)
head(cst_mig)
getwd()
# write.csv(cst_mig, file = "cost_mig_individuals_final.csv")
cst_mig <- read.csv("cost_mig_individuals_final.csv")[,-1]


ggplot(data = cst_mig, aes(x = c_ind, fill = loc)) + geom_histogram() +
  facet_wrap(~mig)


ggplot() + geom_boxplot(data = cst_mig, aes(x = loc, y = c_ind, colour = mig)) +
  geom_boxplot(data = s_p, aes(x = loc, y = cost_ind, colour = mig))

summary(lm(c_ind ~ mig, data = cst_mig))


head(cst_mig)
head(s_p)

sp_t <- s_p[,c(2,4:7)] %>% mutate(mig = ifelse(mig == "cost_aut_ind", "autumn", "spring"),
                                  type = "sim")
colnames(sp_t) <- c("cost", "indiv", "loc", "mig", "c_ind", "type")

cst_mig$type <- "real"

com_d <- rbind(sp_t, cst_mig)


##### Remove birds that don't have return trips
com_d$c_ind[com_d$mig == "spring" & com_d$indiv == "AD"] <- NA
com_d$c_ind[com_d$mig == "spring" & com_d$indiv == "Bird2"] <- NA
com_d$c_ind[com_d$mig == "spring" & com_d$indiv == "Bird7"] <- NA


ggplot(data = com_d, aes(x = loc, y = c_ind, colour = type)) + geom_boxplot() +
  facet_wrap(~mig)


ggplot(data = subset(com_d, type == "sim"), aes(c_ind, fill = loc)) + geom_histogram() +
  geom_vline(data = subset(com_d, type == "real"), aes(xintercept = c_ind, colour = loc)) +
  facet_wrap(mig~loc)

ggplot(data = subset(com_d, type == "sim"), aes(c_ind, fill = "simulated")) + geom_density() +
  geom_density(data = subset(com_d, type == "real"), aes(c_ind, fill = "real")) +
  facet_wrap(mig~loc)


ggplot(data = com_d, aes(x = loc, y = c_ind, fill = type)) + geom_boxplot() +
  facet_wrap(~mig)


## check that the tracks look okay
m_pos_nw <- do.call("rbind", mig_pos)
head(m_pos_nw)

ggplot(data = subset(m_pos_nw, loc == "Sedbergh"), aes(x = lon, y = lat, group = indiv, colour = mig)) + geom_path() +
  geom_point(aes(x = -2.5, y = 54.5), size = 4) +
  geom_point(aes(x = -16, y = 10), size = 4) +
  facet_wrap(mig~indiv) 
ggplot(data = subset(m_pos_nw, loc == "Scotland"), aes(x = lon, y = lat, group = indiv, colour = mig)) + geom_path() +
  geom_point(aes(x = -2.5, y = 54.5), size = 4) +
  geom_point(aes(x = -16, y = 10), size = 4) +
  facet_wrap(~indiv) 
ggplot(data = subset(m_pos_nw, loc == "Senegal"), aes(x = lon, y = lat, group = indiv, colour = mig)) + geom_path() +
  geom_point(aes(x = -2.5, y = 54.5), size = 4) +
  geom_point(aes(x = -16, y = 10), size = 4) +
  facet_wrap(~indiv) 


ggplot(data = subset(mp, loc == "Sedbergh"), aes(x = lon, y = lat, group = indiv, colour = mig)) + geom_path() +
  geom_point(aes(x = -2.5, y = 54.5), size = 4) +
  geom_point(aes(x = -16, y = 10), size = 4) +
  facet_wrap(mig~indiv) 

ggplot(data = cst_mig, aes(x = indiv, y = cost, colour = mig)) + geom_point()


