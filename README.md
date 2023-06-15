
# Analysis of migratory connectivity and wind assistance in the Common Sandpiper

This repository contains the code to run the analyses in the paper "Investigating connectivity and seasonal differences in wind assistance in the migration of Common Sandpipers" - Mondain-Monval _et al._ 2023 (LINK:XXXX). In this, we fitted Common Sandpipers with geolocators in Cumbria and Senegal and tracked them during their migration. We also used the data from Summers _et al._ 2019 who tracked birds from two sites in Scotland (data DOI:10.5441/001/1.9cj0t65k). This paper has two main sections, first we investigated migratory connectivity and overlap between breeding populations using kernel density analyses, and then investigated the influence of wind on migratory birds. We are also the first to report the breeding locations of Common Sandpipers tagged in West Africa.

Data is available from figshare: xxxxx

## Migratory connectivity

This code consists of three scripts, which can be run in any order. A brief summary of each:

1. `PlottingMigrationMaps.R` - Carries out the kernel density analyses for each of the populations, and plots the results
2. `Distances_between_indivs_MantelCorr.R` - Used for the summary table of distances between individuals and the Mantel test
3. `CompareTwoMethodsScot.R` - Compares the outputs of the two different methods of geolocator analysis for the Scottish tagged birds

## Wind analysis

Contains several scripts required to run the wind analysis portion of the analyses. All scripts are coded to run on a PC except for the `submit_simulated_bird_costs_function.R` and `simulated_bird_costs_function.R` which are designed to run on a computer cluster to save time. The scripts must be run in the following order:

1. `calculate_wind_costs.R` - Main script that: downloads wind data across all altitudes, simulates random migratory tracks, simulates geolocator error tracks, and then calculates costs of simulated migrations, real bird migrations and costs accounting for geolocation error.
2. `combine_real_simulated_costs.R` - Combines all the files created in `calculate_wind_costs.R`
3. `wind_analyses_models.R` - Runs all models and outputs tables
4. `wind_analyses_plots.R` - Runs all the plots

There is also the script `check_wind_conditions.R`, which does not need to be run in order (except for after the download). This plots the direction and speed of wind at each altitude during migration window of each individual tracked with a geolocator. 
