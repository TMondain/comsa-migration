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
