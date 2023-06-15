

## Check wind conditions for each bird's migrations

library(raster)
library(tidyverse)


# saveRDS(w_spr_ras, file = "data/wind_analyses/wind_cost_rasters/wind_raster_spring.rds")
w_aut_ras <- readRDS("data/wind_analyses/wind_cost_rasters/wind_raster_autumn.rds")
w_spr_ras <- readRDS("data/wind_analyses/wind_cost_rasters/wind_raster_spring.rds")

# separate autumn direction and speed
dir_aut <- do.call(rbind, lapply(1:length(w_aut_ras), function(x) {
  do.call(rbind, lapply(1:length(w_aut_ras[[x]]), function(i) {
    data.frame(nms = paste(names(w_aut_ras)[x], names(w_spr_ras[[x]])[i], names(w_aut_ras[[x]][[i]])[1]), 
               dr = values(w_aut_ras[[x]][[i]][[1]]))
  }))
}))

spd_aut <- do.call(rbind, lapply(1:length(w_aut_ras), function(x) {
  do.call(rbind, lapply(1:length(w_aut_ras[[x]]), function(i) {
    data.frame(nms = paste(names(w_aut_ras)[x], names(w_spr_ras[[x]])[i], names(w_aut_ras[[x]][[i]])[2]), 
               spd = values(w_aut_ras[[x]][[i]][[2]]))
  }))
}))

# separate spring direction and speed
dir_spr <- do.call(rbind, lapply(1:length(w_spr_ras), function(x) {
  do.call(rbind, lapply(1:length(w_spr_ras[[x]]), function(i) {
    data.frame(nms = paste(names(w_spr_ras)[x], names(w_spr_ras[[x]])[i], names(w_spr_ras[[x]][[i]])[1]), 
               dr = values(w_spr_ras[[x]][[i]][[1]]))
  }))
}))

spd_spr <- do.call(rbind, lapply(1:length(w_spr_ras), function(x) {
  do.call(rbind, lapply(1:length(w_spr_ras[[x]]), function(i) {
    data.frame(nms = paste(names(w_spr_ras)[x], names(w_spr_ras[[x]])[i], names(w_spr_ras[[x]][[i]])[2]), 
               spd = values(w_spr_ras[[x]][[i]][[2]]))
  }))
}))


#### Speed

# autumn 
ggplot(spd_aut, aes(x = nms, y = spd)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# spring 
ggplot(spd_spr, aes(x = nms, y = spd)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### direction

# autumn 
ggplot(dir_aut, aes(x = nms, y = dr)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# spring 
ggplot(dir_spr, aes(x = nms, y = dr)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
