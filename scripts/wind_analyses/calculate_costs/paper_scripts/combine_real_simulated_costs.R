
library(tidyverse)

### combine the simulated and real bird costs

#####     load simulated bird costs     #####
simulated_costs <- read.csv("data/wind_analyses/flight_costs/combined_simulation_costs.csv")

#####     load real bird costs     #####
real_birds_costs <- read.csv("data/wind_analyses/flight_costs/cost_mig_individuals.csv")[,-1]
real_birds_costs$type <- 'real'

head(real_birds_costs)

#####     load geolocation error simulation     #####
glc_err <- read.csv("data/wind_analyses/flight_costs/cost_mig_individuals_simulated_geolocation_error.csv")[,-1]
glc_err$indiv <- NULL
glc_err$type <- "real_gls_err"

colnames(glc_err) <- c("cost", "indiv", "loc", "mig", "c_ind", "type")
head(glc_err)


## edit the simulated birds dataframe for binding
sim_costindex <- simulated_costs[,-1] %>% 
  dplyr::select(real_indiv, sim_indiv, loc, cost_aut_ind, cost_spr_ind) %>%
  pivot_longer(cols = c("cost_aut_ind", "cost_spr_ind"), 
               names_to = "mig", 
               values_to = "cost_ind")

sim_cost <- simulated_costs[,-1] %>% 
  dplyr::select(real_indiv, sim_indiv, loc, cost_aut, cost_spr) %>%
  pivot_longer(cols = c("cost_aut", "cost_spr"), 
               names_to = "mig", 
               values_to = "cost")

simulated_costs_df <- cbind(cost = sim_cost$cost, 
                            sim_costindex) %>% 
  mutate(mig = ifelse(mig == "cost_aut_ind", "autumn", "spring"),
         type = "sim")

# remove the simulated individual name
simulated_costs_df$sim_indiv <- NULL
colnames(simulated_costs_df) <- c("cost", "indiv", "loc", "mig", "c_ind", "type")

head(simulated_costs_df)


#####     create combined real and simulated bird df     #####


com_d <- rbind(simulated_costs_df, real_birds_costs, glc_err)
head(com_d)

write.csv(com_d, 
          file = "data/wind_analyses/combined_real_simulated_costs.csv")
