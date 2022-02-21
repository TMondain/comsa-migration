# rm(list=ls())
# # supposed to look like:
# https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ncep_global.csv?ugrd10m[(2019-12-10T12:00:00Z):1:(2019-12-10T12:00:00Z)][(-90.0):1:(90.0)][(0.0):1:(359.5)],vgrd10m[(2019-12-10T12:00:00Z):1:(2019-12-10T12:00:00Z)][(-90.0):1:(90.0)][(0.0):1:(359.5)]
# 
# 
# # test
# https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ncep_global.csv?ugrd10m[(2013-06-21T12:00:00Z):1:(2013-06-21T12:00:00Z)][(-90.0):1:(90.0)][(0.0):1:(359.5)],vgrd10m[(2013-06-21T12:00:00Z):1:(2013-06-21T12:00:00Z)][(-90.0):1:(90.0)][(0.0):1:(359.5)]
# 
# 
# # my one
# https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ncep_global.csv?ugrd10m[(2013-06-21T00:00:00Z):1:(2013-06-21T00:00:00Z)][(-10.0):1:(72.5)][(0.0):1:(40)], vgrd10m[(2013-06-21T00:00:00Z):1:(2013-06-21T00:00:00Z)][(-10):1:(72.5)][(0.0):1:(40)]


# # 
# dt_a <- "2013-06-21 00:16:30 UTC"
# # 
# t <- w_t(dt_a, -30, 40, -10, 72.5)
# # # 
# time = dt_a
# lon1 = -30
# lon2 = 40
# lat1 = -10
# lat2 = 72.5
# type = "read-data"

# load package rCAT
library(rCAT)

## first wind.fit_int
wind.fit_int <- function (tmpx) {
  tmpx[,3] <- tmpx[,3] %% 360
  tmpx[tmpx[,3]>=180,3] <- tmpx[tmpx[,3]>=180,3] - 360
  
  ###### DIRECTION
  direction <- atan2(tmpx[,4], tmpx[,5])
  direction <- rad2deg(direction)
  direction[direction < 0] <- 360 + direction[direction < 0]
  
  ###### SPEED
  speed <- sqrt( (tmpx[,4] * tmpx[,4]) + (tmpx[,5] * tmpx[,5]))
  
  ######
  names(tmpx)<- c("time", "lat","lon", "ugrd10m", "vgrd10m")
  res <- cbind(tmpx, dir=direction, speed=speed)
  res <- res[with(res, order(-lat)), ]
  res[,1] <- ymd_hms(res[,1], truncated = 3)
  return(res)
}

w_t <- function(time, lon1, lon2, lat1, lat2, type = "read-data",
                 trace = 1) {

  # load package rCAT
  library(rCAT)
  
  ## first wind.fit_int
  wind.fit_int <- function (tmpx) {
    tmpx[,3] <- tmpx[,3] %% 360
    tmpx[tmpx[,3]>=180,3] <- tmpx[tmpx[,3]>=180,3] - 360
    
    ###### DIRECTION
    direction <- atan2(tmpx[,4], tmpx[,5])
    direction <- rad2deg(direction)
    direction[direction < 0] <- 360 + direction[direction < 0]
    
    ###### SPEED
    speed <- sqrt( (tmpx[,4] * tmpx[,4]) + (tmpx[,5] * tmpx[,5]))
    
    ######
    names(tmpx)<- c("time", "lat","lon", "ugrd10m", "vgrd10m")
    res <- cbind(tmpx, dir=direction, speed=speed)
    res <- res[with(res, order(-lat)), ]
    res[,1] <- ymd_hms(res[,1], truncated = 3)
    return(res)
  }
  
  type <- match.arg(type, c("read-data", "csv"))
  dt <- as_datetime(time)
  resultados <- vector("list", length(dt))
  names(resultados) <- dt
  for (id in seq_along(dt)) {
    yyyy_c <- year(dt[id])
    mm_c <- sprintf("%02d", month(dt[id]))
    dd_c <- sprintf("%02d", day(dt[id]))
    tt_c <- sprintf("%02d", hour(dt[id]))
    testDate <- paste(yyyy_c, "-", mm_c, "-", 
                      dd_c, sep = "")
    if (trace) 
      print(paste(ymd_h(paste(yyyy_c, mm_c, dd_c, tt_c, 
                              sep = "-")), "downloading...", sep = " "))
    tryCatch({
      as.Date(testDate)
      if (lon1 < 0) {
        lon1 <- 360 - (abs(lon1))
      }
      if (lon2 < 0) {
        lon2 <- 360 - (abs(lon2))
      }
      if (lon1 > 180 && lon2 < 180) {
        
        url_west <- paste("https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ncep_global.csv?ugrd10m[(",
                          yyyy_c, "-", mm_c, "-", dd_c, "T",tt_c, ":00:00Z):1:(",
                          yyyy_c, "-", mm_c, "-", dd_c, "T", 
                          tt_c, ":00:00Z)][(", lat1, "):1:(", 
                          lat2, ")][(", lon1, "):1:(359.5)],vgrd10m[(", 
                          yyyy_c, "-", mm_c, "-", dd_c, "T", 
                          tt_c, ":00:00Z):1:(",
                          yyyy_c, "-", mm_c, "-", dd_c, "T", 
                          tt_c, ":00:00Z)][(", lat1, "):1:(", lat2, ")][(", lon1, "):1:(359.5)]", 
                          sep = "")
        
        
        url_east <- paste("https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ncep_global.csv?ugrd10m[(", 
                          yyyy_c, "-", mm_c, "-", dd_c, "T", 
                          tt_c, ":00:00Z):1:(",
                          yyyy_c, "-", mm_c, "-", dd_c, "T", 
                          tt_c, ":00:00Z)][(", lat1, "):1:(", lat2, ")][(0.0):1:(", lon2, ")],vgrd10m[(", 
                          yyyy_c, "-", mm_c, "-", dd_c, "T", 
                          tt_c, ":00:00Z):1:(",
                          yyyy_c, "-", mm_c, "-", dd_c, "T", 
                          tt_c, ":00:00Z)][(", lat1, "):1:(", lat2, ")][(0.0):1:(", lon2,")]", 
                          sep = "")
        
        
        tmp <- rbind(read.csv(url_west, header = FALSE, 
                              skip = 2, stringsAsFactors = FALSE), read.csv(url_east, 
                                                                            header = FALSE, skip = 2, stringsAsFactors = FALSE))
        tmp <- wind.fit_int(tmp)
        if (type == "csv") {
          tmp <- wind.fit_int(tmp)
          fname <- paste("wind_", yyyy_c, "_", 
                         mm_c, "_", dd_c, "_", tt_c, ".csv", 
                         sep = "")
          write.table(tmp, fname, sep = ",", row.names = FALSE, 
                      col.names = TRUE, quote = FALSE)
        }
        else {
          resultados[[id]] <- tmp[, 4:5]
        }
      }
      else {
        url_dir <- paste("https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ncep_global.csv?ugrd10m[(", 
                         yyyy_c, "-", mm_c, "-", dd_c, "T", 
                         tt_c, ":00:00Z):1:(",
                         yyyy_c, "-", mm_c, "-", dd_c, "T", 
                         tt_c, ":00:00Z)][(", lat1, "):1:(", lat2, ")][(",
                         lon1, "):1:(", lon2, ")],vgrd10m[(", 
                         yyyy_c, "-", mm_c, "-", dd_c, "T", 
                         tt_c, ":00:00Z):1:(",
                         yyyy_c, "-", mm_c, "-", dd_c, "T", 
                         tt_c, ":00:00Z)][(", lat1, "):1:(", lat2, ")][(", lon1, "):1:(", lon2,")]", 
                         sep = "")
        
        tmp <- read.csv(url_dir, header = FALSE, skip = 2, 
                        colClasses = c("POSIXct", "double", 
                                       "double", "double", "double"))
        tmp <- wind.fit_int(tmp)
        if (type == "csv") {
          tmp <- wind.fit_int(tmp)
          fname <- paste("wind_", yyyy_c, "_", 
                         mm_c, "_", dd_c, "_", tt_c, ".csv", 
                         sep = "")
          write.table(tmp, fname, sep = ",", row.names = FALSE, 
                      col.names = TRUE, quote = FALSE)
        }
        else {
          resultados[[id]] <- tmp[, 4:5]
        }
      }
    }, error = function(e) {
      cat("ERROR: database not found. Please, check server\n                          connection, date or geographical ranges \n")
    }, warning = function(w) {
      cat("ERROR: database not found. Please, check server\n                            connection, date or geographical ranges  \n")
    })
  }
  if (type == "csv") 
    return(NULL)
  attr(resultados, "lat_lon") <- tmp[, 2:3]
  class(resultados) <- c("rWind_series", "list")
  return(resultados)
}
