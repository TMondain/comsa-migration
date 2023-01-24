rm(list = ls())

########################################
#####     Wind analysis models     #####
########################################

library(tidyverse)
library(gridExtra)
library(lme4)
library(MuMIn)
library(cowplot)
library(readr)
library(lmerTest)

com_d <- read.csv('data/wind_analyses/combined_real_simulated_costs.csv')


#################################################################
#####     Table 1 : Does cost differ between migrations     #####
#################################################################

# real birds
mr <- lm(c_ind~mig, data = subset(com_d, type == "real"))
summary(mr)
hist(resid(mr))
par(mfrow = c(2,2))
plot(mr)
par(mfrow = c(1,1))

# simulated birds
ms <- lmer(c_ind~mig + (1|indiv), data = subset(com_d, type == "sim"))
summary(ms)
hist(resid(ms))
plot(ms)

# to get chi2 and p-value
msi <-lmer(c_ind~1 + (1|indiv), data = subset(com_d, type == "sim"))  
anova(ms, msi)


data.frame(Bird_type = c("Observed", "Simulated"), 
           Autumn = c(mr$coefficients[1], summary(ms)$coefficients[1,1]), 
           Spring = c(mr$coefficients[1] + mr$coefficients[2],
                      summary(ms)$coefficients[1,1] + summary(ms)$coefficients[2,1]),
           F_T_statistic = c(summary(mr)$fstatistic[1], summary(ms)$coefficients[2,3]),
           P_value = c(anova(mr)$'Pr(>F)'[1], NA),
           Marginal_R2 = c(NA, r.squaredGLMM(ms)[1]),
           Adjusted_Conditional_R2 = c(summary(mr)$adj.r.squared, r.squaredGLMM(ms)[2]))


# real birds, by location does cost differ between autumn and spring?
# it does, but the sample sizes apart from the cumbria combination are pathetic
# lump together and accept it's not perfect
sed_mig_m <- lm(c_ind~mig, data = subset(com_d, type == "real" & loc == "Sedbergh"))
summary(sed_mig_m)
scot_mig_m <- lm(c_ind~mig, data = subset(com_d, type == "real" & loc == "Scotland"))
summary(scot_mig_m)
sene_mig_m <- lm(c_ind~mig, data = subset(com_d, type == "real" & loc == "Senegal"))
summary(sene_mig_m)

ggplot(com_d, aes(x=type, y = c_ind, fill = mig)) +
  geom_boxplot() +
  facet_wrap(~loc) +
  theme_bw()

###################################################################################
#####     MODEL SET 2 : Does migration cost vary between tagging location     #####
###################################################################################


### Real autumn

ra <- subset(com_d, mig == "autumn" & type == "real")
head(ra)

hist(ra$c_ind)
m_aut_r <- lm(c_ind ~ loc, data = ra)
summary(m_aut_r)

### Real spring

rs <- subset(com_d, mig == "spring" & type == "real")
m_spr_r <- lm(c_ind ~ loc, data = rs)
summary(m_spr_r)


### Simulated autumn
sa <- subset(com_d, mig == "autumn" & type == "sim")
m_aut_s <- lmer(c_ind ~ loc + (1|indiv), data = sa)
summary(m_aut_s)$coefficients


### Simulated spring
ss <- subset(com_d, mig == "spring" & type == "sim")
m_spr_s <- lmer(c_ind ~ loc + (1|indiv), data = ss)
summary(m_spr_s)


######################################################################################
#####     Table 2 : Real and Simulated table by location (Not used in thesis)    #####
######################################################################################


data.frame(Bird_type = rep(c("Observed", "Simulated"), each = 2), Migration = rep(c("Autumn", "Spring"),2), 
           Sedbergh = c(m_aut_r$coefficients[1], 
                        m_spr_r$coefficients[1],
                        summary(m_aut_s)$coefficients[1,1],
                        summary(m_spr_s)$coefficients[1,1]), 
           Senegal = c(m_aut_r$coefficients[1] + m_aut_r$coefficients[2], 
                       m_spr_r$coefficients[1] + m_spr_r$coefficients[2],
                       summary(m_aut_s)$coefficients[1,1] + summary(m_aut_s)$coefficients[2,1],
                       summary(m_spr_s)$coefficients[1,1] + summary(m_spr_s)$coefficients[2,1]),
           Scotland = c(m_aut_r$coefficients[1] + m_aut_r$coefficients[3],
                        m_spr_r$coefficients[1] + m_spr_r$coefficients[3],
                        summary(m_aut_s)$coefficients[1,1] + summary(m_aut_s)$coefficients[3,1],
                        summary(m_spr_s)$coefficients[1,1] + summary(m_spr_s)$coefficients[3,1]),
           F_statistic = c(summary(m_aut_r)$fstatistic[1],
                           summary(m_spr_r)$fstatistic[1], 
                           NA, NA),
           P_value = c(anova(m_aut_r)$'Pr(>F)'[1],
                       anova(m_spr_r)$'Pr(>F)'[1],
                       NA, NA),
           Marginal_R2 = c(NA, NA, r.squaredGLMM(m_aut_s)[1], r.squaredGLMM(m_spr_s)[1]),
           Adjusted_Conditional_R2 = c(summary(m_aut_r)$adj.r.squared,
                                       summary(m_spr_r)$adj.r.squared, 
                                       r.squaredGLMM(m_aut_s)[2], 
                                       r.squaredGLMM(m_spr_s)[2]))



co <- data.frame((summary(m_aut_r))$coefficients)
co_s <- data.frame((summary(m_spr_r))$coefficients)

co_sim <- data.frame((summary(m_aut_s))$coefficients)
co_sim_spr <- data.frame((summary(m_spr_s))$coefficients)

options(scipen = 999)

rst <- data.frame(Model_number = rep(c(1:4), each = 3), 
                  Location = rep(c("Sedbergh", "Senegal", "Scotland"), 4),
                  Bird_type = rep(c("Observed", "Simulated"), each = 6),  
                  Migration = rep(rep(c("Autumn", "Spring"), each = 3), 2),
                  Estimate = c(co$Estimate, co_s$Estimate,
                               co_sim$Estimate, co_sim_spr$Estimate), 
                  Standard_error = c(co$Std..Error,co_s$Std..Error,
                                     co_sim$Std..Error, co_sim_spr$Std..Error), 
                  T_value = c(co$t.value, co_s$t.value, co_sim$t.value, co_sim_spr$t.value), 
                  P_value = c(co$Pr...t..,co_s$Pr...t.., NA, NA, NA, NA, NA, NA),
                  Marginal_R2 = c(rep(NA, 6), 
                                  r.squaredGLMM(m_aut_s)[1], r.squaredGLMM(m_aut_s)[1], r.squaredGLMM(m_aut_s)[1], 
                                  r.squaredGLMM(m_spr_s)[1], r.squaredGLMM(m_spr_s)[1], r.squaredGLMM(m_spr_s)[1]),
                  Adjusted_Conditional_R2 = c(rep(c(summary(m_aut_r)$adj.r.squared, summary(m_spr_r)$adj.r.squared), each = 3),
                                              r.squaredGLMM(m_aut_s)[2], r.squaredGLMM(m_aut_s)[2], r.squaredGLMM(m_aut_s)[2], 
                                              r.squaredGLMM(m_spr_s)[2], r.squaredGLMM(m_spr_s)[2], r.squaredGLMM(m_spr_s)[2]))
rst



########################################################################
#####     Table 3 : Are Real sandpipers better than on average     #####
########################################################################


#####     Using mixed effects modelling      #####

## Observed vs simulated 
sed_a1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Sedbergh" & mig == "autumn")))
sed_a <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Sedbergh" & mig == "autumn")))
sed_a1_i <- (lmer(c_ind ~ 1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Sedbergh" & mig == "autumn")))
ansed_a1 <- anova(sed_a1, sed_a1_i)
ansed_a1$`Pr(>Chisq)`

sen_a1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Senegal" & mig == "autumn")))
sen_a <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Senegal" & mig == "autumn")))
sen_a1_i <- (lmer(c_ind ~ 1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Senegal" & mig == "autumn")))
ansen_a1 <- anova(sen_a1, sen_a1_i)

scot_a1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Scotland" & mig == "autumn")))
scot_a <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Scotland" & mig == "autumn")))
scot_a1_i <- (lmer(c_ind ~ 1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Scotland" & mig == "autumn")))
anscot_a1 <- anova(scot_a1, scot_a1_i)

sed_s1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Sedbergh" & mig == "spring")))
sed_s <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Sedbergh" & mig == "spring")))
sed_s1_i <- (lmer(c_ind ~ 1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Sedbergh" & mig == "spring")))
ansed_s1 <- anova(sed_s1, sed_s1_i)

sen_s1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Senegal" & mig == "spring")))
sen_s <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Senegal" & mig == "spring")))
sen_s1_i <- (lmer(c_ind ~ 1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Senegal" & mig == "spring")))
ansen_s1 <- anova(sen_s1, sen_s1_i)

scot_s1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Scotland" & mig == "spring")))
scot_s <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Scotland" & mig == "spring")))
scot_s1_i <- (lmer(c_ind ~ 1 + (1|indiv), data = subset(com_d, type != "real_gls_err" & loc == "Scotland" & mig == "spring")))
anscot_s1 <- anova(scot_s1, scot_s1_i)


# rbind(sed_a$coefficients,
#       sen_a$coefficients,
#       scot_a$coefficients)

# With two groups for cost reporting
obvs_sim <- data.frame(Location = rep(c("Sedbergh", "Senegal", "Scotland"), 2),
                       Migration = rep(c("Autumn", "Spring"), each = 3),
                       Intercept = c(sed_a$coefficients[1,1],sen_a$coefficients[1,1],scot_a$coefficients[1,1],
                                     sed_s$coefficients[1,1],sen_s$coefficients[1,1],scot_s$coefficients[1,1]), 
                       Simulated_estimate = c(sed_a$coefficients[2,1],sen_a$coefficients[2,1],scot_a$coefficients[2,1],
                                              sed_s$coefficients[2,1],sen_s$coefficients[2,1],scot_s$coefficients[2,1]), 
                       Standard_error = c(sed_a$coefficients[2,2],sen_a$coefficients[2,2],scot_a$coefficients[2,2],
                                          sed_s$coefficients[2,2],sen_s$coefficients[2,2],scot_s$coefficients[2,2]),
                       
                       X2_value = c(ansed_a1$Chisq[2], ansen_a1$Chisq[2], anscot_a1$Chisq[2],
                                    ansed_s1$Chisq[2], ansen_s1$Chisq[2], anscot_s1$Chisq[2]),
                       P_value = c(ansed_a1$`Pr(>Chisq)`[2], ansen_a1$`Pr(>Chisq)`[2], anscot_a1$`Pr(>Chisq)`[2],
                                   ansed_s1$`Pr(>Chisq)`[2], ansen_s1$`Pr(>Chisq)`[2], anscot_s1$`Pr(>Chisq)`[2]),
                       
                       T_value = c(sed_a$coefficients[2,3], sen_a$coefficients[2,3], scot_a$coefficients[2,3],
                                   sed_s$coefficients[2,3], sen_s$coefficients[2,3], scot_s$coefficients[2,3]),
                       RE_variance = c(as.data.frame(VarCorr(sed_a1))[1,4], as.data.frame(VarCorr(sen_a1))[1,4],
                                       as.data.frame(VarCorr(scot_a1))[1,4], as.data.frame(VarCorr(sed_s1))[1,4],
                                       as.data.frame(VarCorr(sen_s1))[1,4], as.data.frame(VarCorr(scot_s1))[1,4]),
                       RE_stdev = c(as.data.frame(VarCorr(sed_a1))[1,5], as.data.frame(VarCorr(sen_a1))[1,5],
                                    as.data.frame(VarCorr(scot_a1))[1,5], as.data.frame(VarCorr(sed_s1))[1,5],
                                    as.data.frame(VarCorr(sen_s1))[1,5], as.data.frame(VarCorr(scot_s1))[1,5]),
                       Marginal_R2 = c(r.squaredGLMM(sed_a1)[1], r.squaredGLMM(sen_a1)[1], r.squaredGLMM(scot_a1)[1], 
                                       r.squaredGLMM(sed_s1)[1], r.squaredGLMM(sen_s1)[1], r.squaredGLMM(scot_s1)[1]),
                       Conditional_R2 = c(r.squaredGLMM(sed_a1)[2], r.squaredGLMM(sen_a1)[2], r.squaredGLMM(scot_a1)[2], 
                                          r.squaredGLMM(sed_s1)[2], r.squaredGLMM(sen_s1)[2], r.squaredGLMM(scot_s1)[2]))
obvs_sim
# write.csv(obvs_sim, file = "AnalysesAndPlots/Outputs/Models_obs_vs_sim_Chi2.csv")



summary(lmer(c_ind ~ type + (1|indiv), data = subset(com_d, mig == "autumn" )))
summary(lmer(c_ind ~ type + (1|indiv), data = subset(com_d, mig == "spring" )))


# with three groups for cost reporting
data.frame(Location = rep(c("Sedbergh", "Senegal", "Scotland"), 2),
           Migration = rep(c("Autumn", "Spring"), each = 3),
           Intercept_observed_bird = c(sed_a$coefficients[1,1],sen_a$coefficients[1,1],scot_a$coefficients[1,1],
                                       sed_s$coefficients[1,1],sen_s$coefficients[1,1],scot_s$coefficients[1,1]), 
           Simulated_estimate = c(sed_a$coefficients[2,1],sen_a$coefficients[2,1],scot_a$coefficients[2,1],
                                  sed_s$coefficients[2,1],sen_s$coefficients[2,1],scot_s$coefficients[2,1]),
           Standard_error = c(sed_a$coefficients[2,2],sen_a$coefficients[2,2],scot_a$coefficients[2,2],
                              sed_s$coefficients[2,2],sen_s$coefficients[2,2],scot_s$coefficients[2,2]),
           T_value = c(sed_a$coefficients[2,3], sen_a$coefficients[2,3], scot_a$coefficients[2,3],
                       sed_s$coefficients[2,3], sen_s$coefficients[2,3], scot_s$coefficients[2,3]),
           RE_variance = c(as.data.frame(VarCorr(sed_a1))[1,4], as.data.frame(VarCorr(sen_a1))[1,4],
                           as.data.frame(VarCorr(scot_a1))[1,4], as.data.frame(VarCorr(sed_s1))[1,4],
                           as.data.frame(VarCorr(sen_s1))[1,4], as.data.frame(VarCorr(scot_s1))[1,4]),
           RE_stdev = c(as.data.frame(VarCorr(sed_a1))[1,5], as.data.frame(VarCorr(sen_a1))[1,5],
                        as.data.frame(VarCorr(scot_a1))[1,5], as.data.frame(VarCorr(sed_s1))[1,5],
                        as.data.frame(VarCorr(sen_s1))[1,5], as.data.frame(VarCorr(scot_s1))[1,5]),
           Marginal_R2 = c(r.squaredGLMM(sed_a1)[1], r.squaredGLMM(sen_a1)[1], r.squaredGLMM(scot_a1)[1], 
                           r.squaredGLMM(sed_s1)[1], r.squaredGLMM(sen_s1)[1], r.squaredGLMM(scot_s1)[1]),
           Conditional_R2 = c(r.squaredGLMM(sed_a1)[2], r.squaredGLMM(sen_a1)[2], r.squaredGLMM(scot_a1)[2], 
                              r.squaredGLMM(sed_s1)[2], r.squaredGLMM(sen_s1)[2], r.squaredGLMM(scot_s1)[2]))



#####     Observed vs geolocation error vs simulated     ##### 
sed_a1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Sedbergh" & mig == "autumn")))
sed_a <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Sedbergh" & mig == "autumn")))
sen_a1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Senegal" & mig == "autumn")))
sen_a <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Senegal" & mig == "autumn")))
scot_a1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Scotland" & mig == "autumn")))
scot_a <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Scotland" & mig == "autumn")))

sed_s1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Sedbergh" & mig == "spring")))
sed_s <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Sedbergh" & mig == "spring")))
sen_s1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Senegal" & mig == "spring")))
sen_s <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Senegal" & mig == "spring")))
scot_s1 <- (lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Scotland" & mig == "spring")))
scot_s <- summary(lmer(c_ind ~ type-1 + (1|indiv), data = subset(com_d, loc == "Scotland" & mig == "spring")))



# rbind(sed_a$coefficients,
#       sen_a$coefficients,
#       scot_a$coefficients)

# With two groups for cost reporting
gls_err <- data.frame(Location = rep(c("Sedbergh", "Senegal", "Scotland"), 2),
                      Migration = rep(c("Autumn", "Spring"), each = 3),
                      Intercept = c(sed_a$coefficients[1,1],sen_a$coefficients[1,1],scot_a$coefficients[1,1],
                                    sed_s$coefficients[1,1],sen_s$coefficients[1,1],scot_s$coefficients[1,1]), 
                      Geolocation_error = c(sed_a$coefficients[2,1],sen_a$coefficients[2,1],scot_a$coefficients[2,1],
                                            sed_s$coefficients[2,1],sen_s$coefficients[2,1],scot_s$coefficients[2,1]), 
                      Simulated_estimate = c(sed_a$coefficients[3,1],sen_a$coefficients[3,1],scot_a$coefficients[3,1],
                                             sed_s$coefficients[3,1],sen_s$coefficients[3,1],scot_s$coefficients[3,1]), 
                      Standard_error = c(sed_a$coefficients[2,2],sen_a$coefficients[2,2],scot_a$coefficients[2,2],
                                         sed_s$coefficients[2,2],sen_s$coefficients[2,2],scot_s$coefficients[2,2]),
                      T_value = c(sed_a$coefficients[2,3], sen_a$coefficients[2,3], scot_a$coefficients[2,3],
                                  sed_s$coefficients[2,3], sen_s$coefficients[2,3], scot_s$coefficients[2,3]),
                      RE_variance = c(as.data.frame(VarCorr(sed_a1))[1,4], as.data.frame(VarCorr(sen_a1))[1,4],
                                      as.data.frame(VarCorr(scot_a1))[1,4], as.data.frame(VarCorr(sed_s1))[1,4],
                                      as.data.frame(VarCorr(sen_s1))[1,4], as.data.frame(VarCorr(scot_s1))[1,4]),
                      RE_stdev = c(as.data.frame(VarCorr(sed_a1))[1,5], as.data.frame(VarCorr(sen_a1))[1,5],
                                   as.data.frame(VarCorr(scot_a1))[1,5], as.data.frame(VarCorr(sed_s1))[1,5],
                                   as.data.frame(VarCorr(sen_s1))[1,5], as.data.frame(VarCorr(scot_s1))[1,5]),
                      Marginal_R2 = c(r.squaredGLMM(sed_a1)[1], r.squaredGLMM(sen_a1)[1], r.squaredGLMM(scot_a1)[1], 
                                      r.squaredGLMM(sed_s1)[1], r.squaredGLMM(sen_s1)[1], r.squaredGLMM(scot_s1)[1]),
                      Conditional_R2 = c(r.squaredGLMM(sed_a1)[2], r.squaredGLMM(sen_a1)[2], r.squaredGLMM(scot_a1)[2], 
                                         r.squaredGLMM(sed_s1)[2], r.squaredGLMM(sen_s1)[2], r.squaredGLMM(scot_s1)[2]))

gls_err

# write.csv(gls_err, file = "AnalysesAndPlots/Outputs/Models_with_gls_err.csv")


summary(lmer(c_ind ~ type + (1|indiv), data = subset(com_d, mig == "autumn" )))
summary(lmer(c_ind ~ type + (1|indiv), data = subset(com_d, mig == "spring" )))


# with three groups for cost reporting
data.frame(Location = rep(c("Sedbergh", "Senegal", "Scotland"), 2),
           Migration = rep(c("Autumn", "Spring"), each = 3),
           Intercept_observed_bird = c(sed_a$coefficients[1,1],sen_a$coefficients[1,1],scot_a$coefficients[1,1],
                                       sed_s$coefficients[1,1],sen_s$coefficients[1,1],scot_s$coefficients[1,1]), 
           Geolocator_error_estimate = c(sed_a$coefficients[2,1],sen_a$coefficients[2,1],scot_a$coefficients[2,1],
                                         sed_s$coefficients[2,1],sen_s$coefficients[2,1],scot_s$coefficients[2,1]),
           Simulated_estimate = c(sed_a$coefficients[3,1],sen_a$coefficients[3,1],scot_a$coefficients[3,1],
                                  sed_s$coefficients[3,1],sen_s$coefficients[3,1],scot_s$coefficients[3,1]), 
           Standard_error = c(sed_a$coefficients[2,2],sen_a$coefficients[2,2],scot_a$coefficients[2,2],
                              sed_s$coefficients[2,2],sen_s$coefficients[2,2],scot_s$coefficients[2,2]),
           T_value = c(sed_a$coefficients[2,3], sen_a$coefficients[2,3], scot_a$coefficients[2,3],
                       sed_s$coefficients[2,3], sen_s$coefficients[2,3], scot_s$coefficients[2,3]),
           RE_variance = c(as.data.frame(VarCorr(sed_a1))[1,4], as.data.frame(VarCorr(sen_a1))[1,4],
                           as.data.frame(VarCorr(scot_a1))[1,4], as.data.frame(VarCorr(sed_s1))[1,4],
                           as.data.frame(VarCorr(sen_s1))[1,4], as.data.frame(VarCorr(scot_s1))[1,4]),
           RE_stdev = c(as.data.frame(VarCorr(sed_a1))[1,5], as.data.frame(VarCorr(sen_a1))[1,5],
                        as.data.frame(VarCorr(scot_a1))[1,5], as.data.frame(VarCorr(sed_s1))[1,5],
                        as.data.frame(VarCorr(sen_s1))[1,5], as.data.frame(VarCorr(scot_s1))[1,5]),
           Marginal_R2 = c(r.squaredGLMM(sed_a1)[1], r.squaredGLMM(sen_a1)[1], r.squaredGLMM(scot_a1)[1], 
                           r.squaredGLMM(sed_s1)[1], r.squaredGLMM(sen_s1)[1], r.squaredGLMM(scot_s1)[1]),
           Conditional_R2 = c(r.squaredGLMM(sed_a1)[2], r.squaredGLMM(sen_a1)[2], r.squaredGLMM(scot_a1)[2], 
                              r.squaredGLMM(sed_s1)[2], r.squaredGLMM(sen_s1)[2], r.squaredGLMM(scot_s1)[2]))



###############################################################################
#####     Table 4 : Comparing duration of autumn and spring migration     #####
###############################################################################

gm <- read_csv("data/movement_data/GLS_mig_R.csv")
head(gm)

md <- gm %>% dplyr::select(aut = aut_dur, spr = spr_dur, loc = loc) %>% 
  pivot_longer(cols = c(aut, spr), names_to = "mig", values_to = "dur")
head(md)

dm <- lm(dur ~ mig*loc, data=md)
summary(dm)
hist(resid(dm))
par(mfrow=c(2,2))
plot(dm)
par(mfrow = c(1,1))

res.aov2 <- aov(dur ~ mig + loc +
                  mig:loc, data = md)
summary(res.aov2)
as.data.frame(summary(res.aov2)[[1]]) %>% 
  rownames_to_column(var = "Estimate")


TUKEY <- TukeyHSD(x=res.aov2, conf.level=0.95)

# Tuckey test representation :
plot(TUKEY , las=1 , col="brown")


dms <- lm(dur ~ mig, data=md)
summary(dms)
hist(resid(dms))
par(mfrow=c(2,2))
plot(dms)
par(mfrow = c(1,1))

# sample size per location
na.omit(subset(md, loc == "Sedbergh")) %>% group_by(mig) %>% tally()
na.omit(subset(md, loc == "Scotland")) %>% group_by(mig) %>% tally()
na.omit(subset(md, loc == "Senegal")) %>% group_by(mig) %>% tally()

dm_sed <- lm(dur ~ mig, data=subset(md, loc == "Sedbergh"))
s_sed <- summary(dm_sed)
s_sed$coefficients
dm_scot <- lm(dur ~ mig, data=subset(md, loc == "Scotland"))
s_scot <- summary(dm_scot)
s_scot$coefficients
dm_sen <- lm(dur ~ mig, data=subset(md, loc == "Senegal"))
s_sen <- summary(dm_sen)
hist(resid(dm_sen))
plot(dm_sen)

s_sen$coefficients

data.frame(Model = c("Migration", "Migration x Location", "Migration x Location", "Migration x Location"), 
           Location = c("All", "Sedbergh", "Scotland", "Senegal"), 
           Autumn = c(dms$coefficients[1], s_sed$coefficients[1,1], s_scot$coefficients[1,1], s_sen$coefficients[1,1]), 
           Spring = c(dms$coefficients[1] + dms$coefficients[2], 
                      s_sed$coefficients[1,1] + s_sed$coefficients[2,1], 
                      s_scot$coefficients[1,1] + s_scot$coefficients[2,1], 
                      s_sen$coefficients[1,1] + s_sen$coefficients[2,1]),
           F_statistic = c(summary(dms)$fstatistic[1], (s_sed)$fstatistic[1], (s_scot)$fstatistic[1], (s_sen)$fstatistic[1]),
           p_value = c(anova(dms)$'Pr(>F)'[1], anova(dm_sed)$'Pr(>F)'[1], 
                       anova(dm_scot)$'Pr(>F)'[1], anova(dm_sen)$'Pr(>F)'[1]),
           Adjusted_R2 = c(summary(dms)$adj.r.squared, (s_sed)$adj.r.squared, (s_scot)$adj.r.squared, (s_sen)$adj.r.squared))


mdur <- ggplot(data = md, aes(x = mig, y = dur)) + geom_boxplot() +
  xlab("Migration") + ylab("Duration (days)") +
  theme_classic() +
  theme(text = element_text(size = 15)) + scale_x_discrete(labels = c("Autumn", "Spring")) 

mdur_loc <- ggplot(data = md, aes(x = mig, y = dur, colour = loc)) + geom_boxplot() +
  scale_color_discrete(labels = c("Scotland", "Sedbergh", "Senegal")) +
  labs(colour = "Location") +
  xlab("Migration") + ylab("Duration (days)") +
  theme_classic() +
  theme(text = element_text(size = 15)) + scale_x_discrete(labels = c("Autumn", "Spring"))

# c_dur <- plot_grid(mdur, mdur_loc, ncol = 2, labels = c("a", "b"))

getwd()
# ggsave(mdur_loc, file = ""AnalysesAndPlots/Outputs/Migration_duration_loc.tiff", device = "tiff", dpi = 600,
#        height=5, width=5)
# ggsave(mdur, file = "AnalysesAndPlots/Outputs/Migration_duration.tiff", device = "tiff", dpi = 600,
#        height=5, width=5)




