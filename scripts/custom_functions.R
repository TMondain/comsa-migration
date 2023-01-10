## Random functions adapted from Finch 2017 supp. matt.

# CUSTOM FUNCTIONS FROM FINCH 2017 ---------------------------------------------------------
# Mean inter-individual distance
medianDist <- function(lon, lat){
  median(lower(distm(cbind(lon, lat), fun = distHaversine)))
}

# Mantel correlation score
# r value
mantelCor <- function(lon1, lat1, lon2, lat2){
  out <- mantel(lower(distm(cbind(lon1, lat1), fun = distHaversine)) ~
                  lower(distm(cbind(lon2, lat2), fun = distHaversine)))
  return(out[c(1)])
}

# p value
mantelCorP <- function(lon1, lat1, lon2, lat2){
  out <- mantel(lower(distm(cbind(lon1, lat1), fun = distHaversine)) ~
                  lower(distm(cbind(lon2, lat2), fun = distHaversine)))
  return(out[c(4)])
}

# Center (without scaling)
center <- function(x){scale(x, scale = F, center = T)}

# moving average function
ma <- function(arr, n=15){
  res = arr
  for(i in n:length(arr)){
    res[i] = mean(arr[(i-n+1):i])
  }
  res
}


mround <- function(x,base){
  base*round(x/base)
}

# get day, month and year sequences

get_dates <- function(date_sequence) {
  
  date_seq <- date_sequence
  month_range <- c(min(month(date_sequence), na.rm = TRUE),
                   max(month(date_sequence), na.rm = TRUE)) #period of months
  year_range <- c(min(year(date_sequence), na.rm = TRUE),
                  max(year(date_sequence), na.rm = TRUE)) #period of years
  
  return(list(date_seq = date_seq,
              month_range = month_range,
              year_range = year_range))
  
}
