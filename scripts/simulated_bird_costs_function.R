## function to run simulated birds on lotus


###########################################################
####     Step 3: Get cost of simulated bird tracks     ####
###########################################################

calculate_sim_costs <- function(real_indiv, indiv) {
  
  library(rWind)
  library(sp)
  library(gdistance)
  library(tidyverse)
  
  # load simulated birds
  sim_t <- readRDS('../../data/simulated_birds_n100.rds')
  
  # load rasters
  fd_aut_out <- readRDS('../../data/wind_cost_rasters/flow_dispersion_autumn.rds')
  fd_spr_out <- readRDS('../../data/wind_cost_rasters/flow_dispersion_spring.rds')
  
  # smooth the tracks of simulated birds to be like
  # real birds
  sim_t <- sim_t %>% group_by(indiv, mig, r_id) %>% 
    mutate(lon = c(smooth(lon, twiceit = TRUE)),
           lat = c(smooth(lat, twiceit = TRUE)))
  
  print(real_indiv)
  
  # get the corresponding raster for each real individual
  fd_aut <- fd_aut_out[[real_indiv]]
  fd_spr <- fd_spr_out[[real_indiv]]
  
  # get the simulated tracks corresponding to each real individual 
  sim_t_r <- subset(sim_t, r_id == real_indiv)
  
  cost_out <- list()
  
  print(indiv)
  
  # get simulated individual 
  sim_i <- subset(sim_t_r, indiv == indiv)
  
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
        cd_spr <- gdistance::costDistance(fd_spr[[sample(1:length(fd_spr), 1)]], 
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
                   sim_indiv = indiv, 
                   real_indiv = real_indiv, 
                   loc = unique(sim_i_m$loc))
  
  # write the csv out
  write.csv(df, 
            paste0('../../outputs/', real_indiv, '_', indiv, '_random_simulated_costs.csv'))
  
}
