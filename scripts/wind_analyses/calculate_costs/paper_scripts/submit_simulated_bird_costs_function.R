## submit the simulated costs of comsa migration function
library(rslurm)
library(gtools)
library(tidyverse)

source('../scripts/simulated_bird_costs_function.R')

sim_t <- readRDS('../data/simulated_birds_n100.rds')

pars <- expand.grid(real_indiv = unique(sim_t$r_id),
                    indiv = unique(sim_t$indiv))

pars <- pars[mixedorder(pars$real_indiv),]


#### slurm apply call
comsa_mig_slurm <- slurm_apply(calculate_sim_costs,
                               params = pars,
                               jobname = 'simulated_comsa_costs',
                               nodes = nrow(pars),
                               cpus_per_node = 1,
                               slurm_options = list(partition = 'short-serial-4hr',
                                                    time = '01:59:59',
                                                    mem = 10000,
                                                    output = "sim_sdm_%a.out",
                                                    error = "sim_sdm_%a.err",
                                                    account = "short4hr"),
                               sh_template = "jasmin_submit_sh.txt",
                               submit = T)
pars$BatchID <- comsa_mig_slurm$jobid
pars$JobID <- 0:(nrow(pars)-1)#slurm job ID
write.csv(pars, paste0('_rslurm_simulated_comsa_costs/pars.csv'))#to match to error files


## for testing
# setwd('tester/')
# 
# real_indiv = pars$real_indiv[1]
# indiv = pars$indiv[1]

## For after for combining
fls <- list.files('outputs/', full.names = TRUE)

combined_outputs <- do.call(rbind, lapply(fls, read.csv))

write.csv(combined_outputs, file = 'outputs/combined_simulation_costs.csv')