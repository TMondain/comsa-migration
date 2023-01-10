
###########################################
####    Functions for wind analysis    ####
###########################################

####################################################
####       Autumn and Spring random tracks      ####
###################################################

simp_sim <- function(start, end, n, print_out = T, move = c("north", "south"), SD = 4) {
  
  mround <- function(x,base){
    base*round(x/base)
  }
  
  lon_out <- list()
  
  reps <- list()
  
  if(move == "south") {
    lat <- seq(from = start[2], to = end[2], by = -0.5)
    lon_shift <- seq(abs(start[1]),abs(end[1]), length = length(lat))
    
    if(start[1]>0){
      lon_shift <- seq((-start[1]),(-end[1]), length = length(lat))
    }
    
  }
  
  if(move == "north") {
    lat <- seq(from = start[2], to = end[2], by = 0.5)
    lon_shift <- seq((start[1]),(end[1]), length = length(lat))
    
  }
  
  
  for(j in 1:n) {
    
    if(print_out == T){
      print(j)
    }
    
    for(i in 1:length(lat)){
      if(move == "south") {
        lon_val <- rnorm(1, mean = -lon_shift[i], sd = SD)
        lon_val
        
        if(lat[i]<end[2]+4 & lat[i]>end[2]+2) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = SD-(SD/2))
        }
        
        if(lat[i]<end[2]+2.5 & lat[i]>end[2]+0.5) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = SD-(SD/3))
        }
        
        if(lat[i]<end[2]+1) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = SD-(SD/4))
        }
        
        while (abs(lon_val) > abs(lon_shift[i])+12) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = 1)
        }
        
        lon_out[[i]] <- mround(lon_val, 0.5) 
      }
      
      if(move == "north") {
        lon_val <- rnorm(1, mean = lon_shift[i], sd = SD)
        lon_val
        
        if(lat[i]>(end[2]-4) & lat[i]<end[2]-2) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = SD-(SD/2))
        }
        
        if(lat[i]>end[2]-2.5 & lat[i]<end[2]-0.5) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = SD-(SD/3))
        }
        
        if(lat[i]>end[2]-0.5) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = SD-(SD/4))
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


################################
####     wind functions     ####
################################

# convert u and v to speed and direction in 0-360 degrees
calc_wind_speed_dir <- function(u, v, wind_180 = FALSE) {
  # Calculate the wind speed
  wind_speed <- sqrt(u^2 + v^2)
  
  # Calculate the wind direction (in degrees)
  wind_dir <- atan2(u, v) * 180 / pi
  
  # Make sure direction is on the scale 0-360
  wind_dir <- ifelse(wind_dir<0, wind_dir + 360, wind_dir)
  
  return(data.frame(wind_speed, wind_dir))
}


# download wind data for specicfic time periods and pressure levels
download_wind <- function(pressure_level, # what level to download data for
                          date_sequence,
                          months_minmax,
                          years_minmax,
                          lat_minmax,
                          lon_minmax,
                          aggreg_function = "mean",
                          return_raster = TRUE) {
  
  require(RNCEP)
  require(terra)
  require(reshape2)
  
  print("!! Downloading v-wind data")
  
  # download vwnd data
  data_vwind <- NCEP.gather(variable = "vwnd",    #name of the variable
                            level = pressure_level, # pressure level 850hPa
                            months.minmax = months_minmax,
                            years.minmax = years_minmax,
                            lat.southnorth = lat_minmax,
                            lon.westeast = lon_minmax,
                            return.units = TRUE,
                            reanalysis2 = TRUE)
  
  # keep only the dates when the birds were actually migrating
  dates <- gsub('-', '_', as.Date(date_seq)) # convert to hms and add underscores
  
  # search for the date names in the names of the 3rd dimension
  data_vwind <- data_vwind[,, grep(paste(dates, collapse = '|'), 
                                   dimnames(data_vwind)[[3]])]
  
  # aggregate - average across all days
  agg_vwind <- NCEP.aggregate(data_vwind, 
                              YEARS = FALSE, 
                              MONTHS = FALSE,
                              DAYS = FALSE,
                              HOURS = FALSE,
                              fxn = aggreg_function)
  
  # pivot to long format
  agg_vwind_lng <- melt(agg_vwind[,,1])
  colnames(agg_vwind_lng) <- c("lat", "lon", "vwind")
  
  # reorder columns
  agg_vwind_lng <- agg_vwind_lng[,c(2,1,3)]
  
  print("!! Downloading u-wind data")
  
  # download uwnd data
  data_uwind <- NCEP.gather(variable = "uwnd",    #name of the variable
                            level = pressure_level, # pressure level 850hPa
                            months.minmax = months_minmax,
                            years.minmax = years_minmax,
                            lat.southnorth = lat_minmax,
                            lon.westeast = lon_minmax,
                            return.units = TRUE,
                            reanalysis2 = TRUE)
  
  # keep only the dates when the birds were actually migrating
  data_uwind <- data_uwind[,, grep(paste(dates, collapse = '|'), 
                                   dimnames(data_uwind)[[3]])]
  
  # aggregate - average across all days
  agg_uwind <- NCEP.aggregate(data_uwind, 
                              YEARS = FALSE, 
                              MONTHS = FALSE,
                              DAYS = FALSE,
                              HOURS = FALSE,
                              fxn = aggreg_function)
  
  # pivot to long format
  agg_uwind_lng <- melt(agg_uwind[,,1])
  colnames(agg_uwind_lng) <- c("lat", "lon", "uwind")
  
  # reorder columns
  agg_uwind_lng <- agg_uwind_lng[,c(2,1,3)]
  
  # bind the two datasets together
  wind_df <- cbind(agg_vwind_lng, uwind = agg_uwind_lng$uwind)
  head(wind_df)
  
  print("!! calculating wind speed and direction")
  dirsp <- calc_wind_speed_dir(u = wind_df$uwind, 
                               v = wind_df$vwind)
  
  # bind speed and direction to wind
  wind_df <- cbind(wind_df, dirsp)
  
  if(return_raster) {
    wnd_rst <- rast(wind_df, type = 'xyz')
    return(wnd_rst)
  } else {
    return(wind_df)
  }
}



