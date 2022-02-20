rm(list = ls())

## correlations between scotland FlightR and Geolight

library(tidyverse) 
library(lubridate)

dp <- read.csv("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Secondary data/Summers et al Data/DeployRetrieveData.csv")
head(dp)

ms <- read_csv("C:/Users/tmond/OneDrive - Lancaster University/PhD Work/Analyses/GLS analyses/GLS_Summary_Movements_Stopovers/GLS_mig_R2.csv")
ms <- ms[ms$loc == "Scotland",]
head(ms)

dp <- dp %>% 
  mutate(ms_brdep = date(ymd_hms(ms$brd_dep[match(animal.id, ms$indiv)])),
         ms_wiarr = date(ymd_hms(ms$win_arr[match(animal.id, ms$indiv)])),
         ms_widep = date(ymd_hms(ms$win_dep[match(animal.id, ms$indiv)])),
         ms_brarr = date(ymd_hms(ms$brd_arr[match(animal.id, ms$indiv)])))
dp

dp <- dp %>% mutate(bd_r = yday(dmy(last_bred_day)),
                    wa_r = yday(dmy(first_afr_day)),
                    wd_r = yday(dmy(last_afr_day)),
                    ba_r = yday(dmy(first_day_scot_spring)),
                    
                    bd = yday(ms_brdep),
                    wa = yday(ms_wiarr),
                    wd = yday(ms_widep),
                    ba = yday(ms_brarr))


cor.test(y = dp$bd, x = dp$bd_r)
cor.test(y = dp$wa, x = dp$wa_r)
cor.test(y = dp$wd, x = dp$wd_r, use =  "complete.obs")
cor.test(y = dp$ba, x = dp$ba_r, use =  "complete.obs")




ggplot(data = dp, aes(x = bd_r, y = bd)) + geom_point()

