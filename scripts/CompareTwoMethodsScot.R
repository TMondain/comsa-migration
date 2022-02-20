rm(list = ls())

## correlations between scotland FlightR and Geolight

library(tidyverse) 
library(lubridate)

# load Scottish deployment and retrieval data (from Summers et al. 2019 paper)
dp <- read.csv("data/DeployRetrieveData_Scotland.csv")
head(dp)

# load my calculation of time estimates etc
ms <- read_csv("data/movement_data/GLS_mig_R2.csv")
ms <- ms[ms$loc == "Scotland",]
head(ms)

# get my calculation of breeding departure, winter arrival + depature, and breeding arrival
dp <- dp %>% 
  mutate(ms_brdep = date(ymd_hms(ms$brd_dep[match(animal.id, ms$indiv)])),
         ms_wiarr = date(ymd_hms(ms$win_arr[match(animal.id, ms$indiv)])),
         ms_widep = date(ymd_hms(ms$win_dep[match(animal.id, ms$indiv)])),
         ms_brarr = date(ymd_hms(ms$brd_arr[match(animal.id, ms$indiv)])))

# convert to DOY
dp <- dp %>% mutate(bd_r = yday(dmy(last_bred_day)),
                    wa_r = yday(dmy(first_afr_day)),
                    wd_r = yday(dmy(last_afr_day)),
                    ba_r = yday(dmy(first_day_scot_spring)),
                    
                    bd = yday(ms_brdep),
                    wa = yday(ms_wiarr),
                    wd = yday(ms_widep),
                    ba = yday(ms_brarr))


## check correlations between dates calculated from the two different methods
# breeding departure
cor.test(y = dp$bd, x = dp$bd_r)

# winter arrival
cor.test(y = dp$wa, x = dp$wa_r)

# winter departure
cor.test(y = dp$wd, x = dp$wd_r, use =  "complete.obs")

# breeding arrival
cor.test(y = dp$ba, x = dp$ba_r, use =  "complete.obs")


