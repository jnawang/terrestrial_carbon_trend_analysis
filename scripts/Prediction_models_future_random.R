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

climate_vegs <- c('cold_evergreen_forest', 'cold_deciduous_forest', 'cold_non_forest', 'dry', 'wet-warm')
key_vars <- c('CO2 + LE_gs_at + TA_gs_at + TA_ngs_at', 
              'CO2 + SW_gs_at + Lgs_at + TA_ngs_at',        # we could remove CO2 and TA_ngs_at
              'P_gs_at + CO2',                              # we should remove CO2
              'LE_gs_at + TA_gs_at + SW_gs_at + P_ngs_at + SWC_gs_at', 
              'CO2 + TA_at + SW_at + DIST_at')

# training performance; testing performance
df_performance <- data.frame()
# initial weight of rf method
weight_rf <- 0.5

# parameters of rf models
m_try = 4
min_node_size = 5

# 5 fold cross validation
set.seed(123)
nfolder <- 5

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

#
for (iclimate_veg in 1:length(climate_vegs)) {
  climate_veg <- climate_vegs[iclimate_veg]
  key_var <- key_vars[iclimate_veg]
  #
  print(climate_veg)
  #
  if (climate_veg == 'cold_evergreen_forest') {
    data_cold_evergreen_forest <- trend_driver %>% 
      filter(IGBP %in% c('ENF', 'EBF') & category == "cold")
  } else if (climate_veg == 'cold_deciduous_forest') {
    data_cold_evergreen_forest <- trend_driver %>% 
      filter(IGBP %in% c('DBF', 'DNF', 'MF') & category == "cold")  %>% filter(!is.na(LE_gs_at))
  } else if (climate_veg == 'cold_non_forest') {
    data_cold_evergreen_forest <- trend_driver %>% 
      filter(!IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold")  %>% dplyr::filter(!is.na(LE_gs_at))
  } else if (climate_veg == 'dry') {
    data_cold_evergreen_forest <- trend_driver %>% 
      filter(category == "dry") %>% dplyr::filter(!is.na(LE_gs_at))
  } else if (climate_veg == 'wet-warm') {
    data_cold_evergreen_forest <- trend_driver %>% 
      filter(category == "wet") %>% dplyr::filter(!is.na(LE_gs_at))
  }
  #
  # I need to use five-fold CV for both methods
  sites <- unique(data_cold_evergreen_forest$site_ID)
  #
  folds <- createFolds(1:length(data_cold_evergreen_forest), k = nfolder, list = TRUE)
  # group_cv <- rsample::group_vfold_cv(data = data_cold_evergreen_forest, group = site_ID, v = nfolder)  #, strata = "NEE"
  # strata means that the distribution of that variable similar across folds. 
  
  for (i in 1:nfolder) {
    data_cold_evergreen_forest_cv <- data_cold_evergreen_forest[-folds[[i]], ]
    #----------------------------------Linear mixed models only-------------------
    # train the model
    fit <- brm(
      bf(
        formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + NEE + GPP + LE + ER + (1 | IGBP) + (1 | Climate.class) + (1 | source)"),
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
      formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + NEE + GPP + LE + ER + IGBP + Climate.class + source"), 
      data = data_cold_evergreen_forest_cv,
      mtry = m_try, 
      min.node.size = min_node_size
    )
    pred_rf <- predict(rf_model, data = data_cold_evergreen_forest)
    plot(pred_rf$predictions, data_cold_evergreen_forest$NEE_at)
    
    
    # weight should be given based on performance
    
    
    #-------------------------------combination of lmm and RF---------------------
    df_performance <- bind_rows(df_performance, data.frame(climate_veg = climate_veg, site_ID = data_cold_evergreen_forest$site_ID, 
                                                           obs = data_cold_evergreen_forest$NEE_at, 
                                                           lmm = pred_lmm, rf = pred_rf$predictions, 
                                                           lmm_rf = pred_lmm * (1 - weight_rf) + pred_rf$predictions * weight_rf, 
                                                           test = (1:nrow(data_cold_evergreen_forest)) %in% folds[[i]]))
    
    
    
    
  }  

}

# I can check results now. 
df_performance %>% filter(test) %>% group_by(climate_veg) %>% summarise(rmse_lmm = sqrt(mean((obs-lmm)^2)),
                                                                        rmse_rf = sqrt(mean((obs-rf)^2)),
                                                                        rmse_lmm_rf = sqrt(mean((obs-lmm_rf)^2)))
df_performance %>% filter(!test) %>% group_by(climate_veg) %>% summarise(rmse_lmm = sqrt(mean((obs-lmm)^2)),
                                                                        rmse_rf = sqrt(mean((obs-rf)^2)),
                                                                        rmse_lmm_rf = sqrt(mean((obs-lmm_rf)^2)))

for (i in 1:length(climate_vegs)) {
  print(climate_vegs[i])
  cor(df_performance$obs[df_performance$test & df_performance$climate_veg == climate_vegs[i]], df_performance$rf[df_performance$test & df_performance$climate_veg == climate_vegs[i]])^2      # & df_performance$rf < 15
  cor(df_performance$obs[df_performance$test & df_performance$climate_veg == climate_vegs[i]], df_performance$lmm[df_performance$test & df_performance$climate_veg == climate_vegs[i]])^2
  cor(df_performance$obs[df_performance$test & df_performance$climate_veg == climate_vegs[i]], df_performance$lmm_rf[df_performance$test & df_performance$climate_veg == climate_vegs[i]])^2
  
  cor(df_performance$obs[!df_performance$test & df_performance$climate_veg == climate_vegs[i]], df_performance$rf[!df_performance$test & df_performance$climate_veg == climate_vegs[i]])^2
  cor(df_performance$obs[!df_performance$test & df_performance$climate_veg == climate_vegs[i]], df_performance$lmm[!df_performance$test  & df_performance$climate_veg == climate_vegs[i]])^2
  cor(df_performance$obs[!df_performance$test & df_performance$climate_veg == climate_vegs[i]], df_performance$lmm_rf[!df_performance$test  & df_performance$climate_veg == climate_vegs[i]])^2
  
  plot(df_performance$obs[df_performance$test], df_performance$rf[df_performance$test])
  plot(df_performance$obs[df_performance$test], df_performance$lmm[df_performance$test])
}




