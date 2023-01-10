
#load the packages
library(RNCEP)
library(lubridate) #date and time manipulation
library(tidyverse) #data manipulation and visualization
library(RColorBrewer) #color schemes
library(sf) #to import a spatial object and to work with geom_sf in ggplot2
library(reshape2)
library(viridis)
library(terra)

library(rWind)
library(gdistance)

# Define the function
calc_wind_speed_dir <- function(u, v, wind_180 = FALSE) {
  # Calculate the wind speed
  wind_speed <- sqrt(u^2 + v^2)
  
  # Calculate the wind direction (in degrees)
  wind_dir <- atan2(u, v) * 180 / pi
  
  # Make sure direction is on the scale 0-360
  wind_dir <- ifelse(wind_dir<0, wind_dir + 360, wind_dir)
  
  return(data.frame(wind_speed, wind_dir))
}


download_wind <- function(pressure_level, # what level to download data for
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



#### workflow

#define the necessary arguments
month_range <- c(1,12)     #period of months
year_range <- c(2016,2016) #period of years

lat_range <- c(-10, 72.5)      #latitude range
lon_range <- c(-30, 40)     #longitude range

#### function not getting any values past 40 lat!!!!

pressure_levels <- c(1000, 925, 850, 700) # sea level, 779, 1502 and 3130 m a.s.l


wind_flow <- lapply(pressure_levels, FUN = function(x) {
  
  wnd_dl <- download_wind(pressure_level = x,
                          months_minmax = month_range,
                          years_minmax = year_range,
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

crds <- data.frame(lon = c(-9.1, -7.3, -8, -3, 5), lat = c(22.1, 32.1, 41, 59, 67))

names(wind_flow) <- pressure_levels
names(wind_flow)

plot(wind_flow[[1]]$wind_rast$direction)
points(crds)


### calculate min cost
pressure_level_df <- data.frame(crd_pos = 1,
                                lon = crds[1,1],
                                lat = crds[1,2],
                                pressure_index = NA,
                                pressure_level = NA,
                                cost = NA,
                                all_costs = NA)

for(x in 2:nrow(crds)){
  
  dispersal_costs <- 
    sapply(1:length(wind_flow), FUN = function(i){ 
      gdistance::costDistance(wind_flow[[i]]$flow_dispersion, 
                              fromCoords = SpatialPoints(crds[x-1,]), 
                              toCoords = SpatialPoints(crds[x,]))
    })
  
  cost <- min(dispersal_costs)
  
  pressure_level_df <- rbind(pressure_level_df, 
                             data.frame(crd_pos = x,
                                        lon = crds[x,1],
                                        lat = crds[x,2],
                                        pressure_index = which.min(dispersal_costs),
                                        pressure_level = pressure_levels[which.min(dispersal_costs)],
                                        cost = cost,
                                        all_costs = paste(round(dispersal_costs, 1), collapse = '_')))
                             
}

pressure_level_df





##### old code

data <- NCEP.gather(variable = "vwnd",    #name of the variable
                    level = 850, # pressure level 850hPa
                    months.minmax = month_range,
                    years.minmax = year_range,
                    lat.southnorth = lat_range,
                    lon.westeast = lon_range,
                    return.units = TRUE,
                    reanalysis2 = TRUE)

av_mnth_vwwind <- NCEP.aggregate(data, 
                                 YEARS = FALSE, 
                                 MONTHS = FALSE,
                                 DAYS = FALSE,
                                 HOURS = FALSE,
                                 fxn = "mean")

image(av_mnth_vwwind[,,1],asp=1)


av_mnth_vwwind_lng <- melt(av_mnth_vwwind[,,1])
colnames(av_mnth_vwwind_lng) <- c("lat", "lon", "vwind")

datauwnd <- NCEP.gather("uwnd",    #name of the variable
                        850, # pressure level 850hPa
                        months.minmax = month_range,
                        years.minmax = year_range,
                        lat.southnorth = lat_range,
                        lon.westeast = lon_range,
                        return.units = TRUE,
                        reanalysis2 = TRUE)

av_mnth_uwwind <- NCEP.aggregate(datauwnd, 
                                 YEARS = FALSE, 
                                 MONTHS = FALSE,
                                 DAYS = FALSE,
                                 HOURS = FALSE,
                                 fxn = "mean")


av_mnth_uwwind_lng <- melt(av_mnth_uwwind[,,1])
colnames(av_mnth_uwwind_lng) <- c("lat", "lon", "uwind")

# they're identical
identical(av_mnth_vwwind_lng$lat, av_mnth_uwwind_lng$lat)
identical(av_mnth_vwwind_lng$lon, av_mnth_uwwind_lng$lon)

wind_df <- cbind(av_mnth_vwwind_lng, uwind = av_mnth_uwwind_lng$uwind)
head(wind_df)

dirsp <- calc_wind_speed_dir(u = wind_df$uwind, 
                             v = wind_df$vwind)

wind_df <- cbind(wind_df, dirsp)

ggplot(wind_df, aes(x=lon, y = lat, fill = wind_speed)) +
  geom_tile()

wnd_rst <- rast(wind_df, type = 'xyz')
wnd_rst

plot(wnd_rst)

wnd <- raster::stack(c(wnd_rst$wind_dir,wnd_rst$wind_speed))
names(wnd) <- c("direction", "speed")

# rwind
fd_aut <- flow.dispersion(wnd, 
                          type="active",
                          output="transitionLayer")
fd_aut <- gdistance::geoCorrection(fd_aut, type="r", multpl=FALSE, scl=TRUE)

### so could basically get the NCEP.aggregate data for each elevation
### then go through each point in turn and query the relevant layers

NCEP.interp()


# remotes::install_github('ErikKusch/KrigR')
library(KrigR)


