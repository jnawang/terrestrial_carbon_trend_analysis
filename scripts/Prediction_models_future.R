# Develop ensemble models to project future NEE trend. 
# Junna Wang, March 25, 2026
# Tomorrow's task: I need to double-check disturbance effects; 
# I need assess model accuracy using cross-validation;


rm(list=ls())

library(librarian)
shelf(tidyverse, ggpubr, lme4, dagitty, lmerTest, brms, randomForest, glmmTMB, performance, ranger, caret, Metrics, rsample)

trend_driver <- read.csv(file.path("data", "potential_trend_driver_sen.csv"))
sites_long <- read.csv(file.path('data', 'long_term_ec_sites.csv'))

# adding IGBP, climate, and category into trend driver
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source")], by = "site_ID")
}

#-------Models for evergreen forests
data_cold_evergreen_forest <- trend_driver %>% 
  filter(IGBP %in% c('ENF', 'EBF') & category == "cold")

# use leaving one out method; 
# I need to use leave-one-out CV for both methods
# training performance; testing performance
df_performance <- data.frame()
#
sites <- unique(data_cold_evergreen_forest$site_ID)
#
set.seed(123)
nfolder <- 5
# folds <- createFolds(1:length(sites), k = nfolder, list = TRUE)
group_cv <- rsample::group_vfold_cv(data = data_cold_evergreen_forest, group = site_ID, v = nfolder)  #, strata = "NEE"
# strata means that the distribution of that variable similar across folds. 
#
weight_rf <- 0.5

# parameters of rf models
m_try = 4
min_node_size = 5

#
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

for (i in 1:nfolder) {
  sites_fold <- unique(assessment(group_cv$splits[[i]])$site_ID)
  data_cold_evergreen_forest_cv <- data_cold_evergreen_forest %>% filter(!site_ID %in% sites_fold)
  #----------------------------------Linear mixed models only-------------------
  # train the model
  fit <- brm(
    bf(
      formula = NEE_at ~ CO2 + LE_gs_at + TA_gs_at + TA_ngs_at + MAP + MAT + NEE + GPP + LE + ER + (1 | IGBP) + (1 | Climate.class) + (1 | source),
      autocor = cor_ar(~ year_mid | site_ID, p = 1)
    ),
    data = data_cold_evergreen_forest_cv,
    backend = "cmdstanr",
    iter = 2000,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    silent = TRUE
  )
  # summary(fit)
  # do the projection
  data_new <- data_cold_evergreen_forest
  data_new$site_ID = data_cold_evergreen_forest_cv$site_ID[1]          # do not use site_ID information during projection. 
  ep <- posterior_epred(fit, newdata = data_new, allow_new_levels = TRUE, incl_autocor = FALSE)
  pred_lmm <- apply(ep, 2, mean)
  plot(pred_lmm, data_cold_evergreen_forest$NEE_at)
  
  
  #----------------------------------Random forests only------------------------
  rf_model <- ranger(
    NEE_at ~ CO2 + LE_gs_at + TA_gs_at + TA_ngs_at + MAP + MAT + NEE + GPP + LE + ER + IGBP + Climate.class + source, 
    data = data_cold_evergreen_forest_cv, 
    mtry = m_try, 
    min.node.size = min_node_size
  )
  pred_rf <- predict(rf_model, data = data_cold_evergreen_forest)
  plot(pred_rf$predictions, data_cold_evergreen_forest$NEE_at)
  #-------------------------------combination of lmm and RF---------------------
  df_performance <- bind_rows(df_performance, data.frame(climate_veg = 'cold_evergreen_forest', site_ID = data_cold_evergreen_forest$site_ID, 
                                                         obs = data_cold_evergreen_forest$NEE_at, 
                                                         lmm = pred_lmm, rf = pred_rf$predictions, 
                                                         lmm_rf = pred_lmm * (1 - weight_rf) + pred_rf$predictions * weight_rf, 
                                                         test = data_cold_evergreen_forest$site_ID %in% sites_fold))
}

# I can check results now. 
df_performance %>% filter(test) %>% group_by(climate_veg) %>% summarise(rmse_lmm = sqrt(mean((obs-lmm)^2)),
                                                                        rmse_rf = sqrt(mean((obs-rf)^2)),
                                                                        rmse_lmm_rf = sqrt(mean((obs-lmm_rf)^2)))
df_performance %>% filter(!test) %>% group_by(climate_veg) %>% summarise(rmse_lmm = sqrt(mean((obs-lmm)^2)),
                                                                        rmse_rf = sqrt(mean((obs-rf)^2)),
                                                                        rmse_lmm_rf = sqrt(mean((obs-lmm_rf)^2)))

cor(df_performance$obs[df_performance$test], df_performance$rf[df_performance$test])^2
cor(df_performance$obs[df_performance$test], df_performance$lmm[df_performance$test])^2
cor(df_performance$obs[df_performance$test], df_performance$lmm_rf[df_performance$test])^2

cor(df_performance$obs[!df_performance$test], df_performance$rf[!df_performance$test])^2
cor(df_performance$obs[!df_performance$test], df_performance$lmm[!df_performance$test])^2
cor(df_performance$obs[!df_performance$test], df_performance$lmm_rf[!df_performance$test])^2

plot(df_performance$obs[df_performance$test], df_performance$rf[df_performance$test])
plot(df_performance$obs[df_performance$test], df_performance$lmm[df_performance$test])

