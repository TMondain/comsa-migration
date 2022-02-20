
#############################################################
####    Functions Autumn and spring simulating tracks    ####
#############################################################

##################################################
###       Autumn and Spring random tracks      ###
##################################################

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
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = SD - (SD/2))
        }
        
        if(lat[i]<end[2]+2.5 & lat[i]>end[2]+0.5) {
          lon_val <- rnorm(1, mean = -lon_shift[i], sd = SD - (SD/3))
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
        lon_val <- rnorm(1, mean = lon_shift[i], sd = SD)
        lon_val
        
        if(lat[i]>(end[2]-4) & lat[i]<end[2]-2) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = SD-(SD/2))
        }
        
        if(lat[i]>end[2]-2.5 & lat[i]<end[2]-0.5) {
          lon_val <- rnorm(1, mean = lon_shift[i], sd = SD-(SD/3))
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




