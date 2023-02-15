
# Analysis of migratory connectivity and wind assistance in the Common Sandpiper

This repository contains the code to run the analyses in the paper "Investigating connectivity and seasonal differences in wind assistance in the migration of Common Sandpipers" - Mondain-Monval _et al._ 2023 (LINK:XXXX). In this, we fitted Common Sandpipers with geolocators in Cumbria and Senegal and tracked them during their migration. We also used the data from Summers _et al._ 2019 who tracked birds from two sites in Scotland (data DOI:10.5441/001/1.9cj0t65k). This paper has two main sections, first we investigated migratory connectivity and overlap between breeding populations using kernel density analyses, and then investigated the influence of wind on migratory birds. We are also the first to report the breeding locations of Common Sandpipers tagged in West Africa.

Data is available from: xxxxx

## Migratory connectivity

This code consists of three scripts, which can be run in any order. A brief summary of each:

1. PlottingMigrationMaps.R - Carries out the kernel density analyses for each of the populations, and plots the results
2. Distances_between_indivs_MantelCorr.R - Used for the summary table of distances between individuals and the Mantel test
3. CompareTwoMethodsScot.R - Compares the outputs of the two different methods of geolocator analysis for the Scottish tagged birds

## Wind analysis

Contains several scripts required to run the wind analysis portion of the analyses. These must be run in the following order:

1. 
