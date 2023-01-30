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


tab_3 <- data.frame(Bird_type = c("Observed", "Simulated"), 
                    Autumn = round(c(mr$coefficients[1], summary(ms)$coefficients[1,1]),2), 
                    Spring = round(c(mr$coefficients[1] + mr$coefficients[2],
                                     summary(ms)$coefficients[1,1] + summary(ms)$coefficients[2,1]), 2),
                    Test_statistic = round(c(summary(mr)$fstatistic[1], summary(ms)$coefficients[2,3]), 2),
                    P_value = round(c(anova(mr)$'Pr(>F)'[1], NA), 2),
                    Marginal_R2 = round(c(NA, r.squaredGLMM(ms)[1]), 2),
                    Adjusted_Conditional_R2 = round(c(summary(mr)$adj.r.squared, r.squaredGLMM(ms)[2]), 2))
tab_3

write.csv(tab_3, 
          file = "outputs/wind_analyses/tab_3_costs_between_migrations.csv")


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


obvs_sim <- data.frame(Location = rep(c("Cumbria", "Senegal", "Scotland"), 2),
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
write.csv(obvs_sim, 
          file = "outputs/wind_analyses/tab_4_real_vs_random_costs.csv")


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
gls_err <- data.frame(Location = rep(c("Cumbria", "Senegal", "Scotland"), 2),
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

# write.csv(gls_err, file = "outputs/wind_analyses/tab_4S_migration_costs_gls_error.csv")



