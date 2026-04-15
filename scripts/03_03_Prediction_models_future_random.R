# Develop ensemble models to project future NEE trend and use cross-validation to test model accuracy. 
# Junna Wang, March 25, 2026
# Tomorrow's task: I need to double-check disturbance effects; 
# I need assess model accuracy using cross-validation;

# I still need to update these models with new driver analysis. Junna

rm(list=ls())

library(librarian)
shelf(tidyverse, ggpubr, lme4, dagitty, lmerTest, brms, randomForest, glmmTMB, performance, ranger, caret, Metrics, rsample, terra)

trend_driver <- read.csv(file.path("data", "potential_trend_driver.csv"))
sites_long <- read.csv(file.path('data', 'long_term_ec_sites.csv'))

# adding IGBP, climate, and category into trend driver
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "AI")], by = "site_ID")
}

#####--------
# bio4 <- terra::rast('/Users/junnawang/UCDLab/Biodiversity/data/climate/current/bio4.tif')
# bio15 <- terra::rast('/Users/junnawang/UCDLab/Biodiversity/data/climate/current/bio15.tif')
# 
# bio4_interannual <- terra::rast('/Users/junnawang/UCDLab/Biodiversity/data/climate/current/inter_annual_variance_bio4.tif')
# bio15_interannual <- terra::rast('/Users/junnawang/UCDLab/Biodiversity/data/climate/current/inter_annual_variance_bio15.tif')
# 
# 
# bio4  <- terra::project(bio4, "EPSG:4326")
# bio15 <- terra::project(bio15, "EPSG:4326")
# 
# bio4_interannual <- terra::project(bio4_interannual, "EPSG:4326")
# bio15_interannual <- terra::project(bio4_interannual, "EPSG:4326")
# 
# sites_long$bio4_ia <- terra::extract(bio4_interannual, data.frame(x = sites_long$Long, y = sites_long$Lat))[, 2]
# sites_long$bio15_ia <- terra::extract(bio15_interannual, data.frame(x = sites_long$Long, y = sites_long$Lat))[, 2]
# 
# sites_long$bio4 <- terra::extract(bio4, data.frame(x = sites_long$Long, y = sites_long$Lat))[, 2]
# sites_long$bio15 <- terra::extract(bio15, data.frame(x = sites_long$Long, y = sites_long$Lat))[, 2]
# 
# 
# if (!'bio4' %in% names(trend_driver)) {
#   trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "bio4", "bio15", "bio4_ia", "bio15_ia")], by = "site_ID")
# }


# we only use significant variables
climate_vegs <- c('cold_evergreen_forest', 'cold_deciduous_forest', 'cold_non_forest', 'dry', 'wet-warm')
key_vars <- c('CO2 + LE_gs_at + TA_gs_at + TA_ngs_at', 
              'SW_gs_at + Lgs_at + TA_ngs_at + CO2', 
              'P_gs_at',                          
              'LE_gs_at + TA_gs_at + SW_gs_at + P_ngs_at + SWC_gs_at', 
              'CO2 + TA_at + SW_at + DIST_at')

# training performance; testing performance
df_performance <- data.frame()

# parameters of rf models
m_try = 4              # m_try = 4 is better than m_try = 3 and 5 for cold non-forest first. 
min_node_size = 5

# 5 fold cross validation
set.seed(128)
nfolder <- 5

options(mc.cores = parallel::detectCores())

#
for (iclimate_veg in 1:length(climate_vegs)) {
  climate_veg <- climate_vegs[iclimate_veg]
  key_var <- key_vars[iclimate_veg]
  #
  print(climate_veg)
  #
  if (climate_veg == 'cold_evergreen_forest') {
    data_group <- trend_driver %>% 
      filter(IGBP %in% c('ENF', 'EBF') & category == "cold") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'cold_deciduous_forest') {
    data_group <- trend_driver %>% 
      filter(IGBP %in% c('DBF', 'DNF', 'MF') & category == "cold") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'cold_non_forest') {
    data_group <- trend_driver %>% 
      filter(!IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'dry') {
    data_group <- trend_driver %>% 
      filter(category == "dry") %>% dplyr::filter(!is.na(LE_gs_at)) 
  } else if (climate_veg == 'wet-warm') {
    data_group <- trend_driver %>% 
      filter(category == "wet") %>% dplyr::filter(!is.na(LE_gs_at)) 
  }
  
  #
  # I need to use five-fold CV for both methods
  sites <- unique(data_group$site_ID)
  #
  folds <- createFolds(1:length(data_group), k = nfolder, list = TRUE)
  # group_cv <- rsample::group_vfold_cv(data = data_group, group = site_ID, v = nfolder)  #, strata = "NEE"
  # strata means that the distribution of that variable similar across folds. 
  
  for (i in 1:nfolder) {
    data_group_cv <- data_group[-folds[[i]], ]
    #--------------------------------------To predict NEE_at----------------------
    #----------------------------------Linear mixed models only-------------------
    # train the model
    # fit <- brm(
    #   bf(
    #     formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + NEE + GPP + LE + (1 | IGBP) + (1 | Climate.class) + (1 | source)"),
    #     autocor = cor_ar(~ year_mid | site_ID, p = 1)
    #   ),
    #   data = data_group_cv,
    #   backend = "cmdstanr",
    #   iter = 1000,
    #   control = list(adapt_delta = 0.99, max_treedepth = 15),
    #   silent = TRUE
    # )
    # # summary(fit)
    # # do the projection
    # data_new <- data_group
    # data_new$site_ID = data_group_cv$site_ID[1]          # do not use site_ID information during projection.
    # ep <- posterior_epred(fit, newdata = data_new, allow_new_levels = TRUE, incl_autocor = FALSE)
    # pred_lmm <- apply(ep, 2, mean)
    # plot(pred_lmm, data_group$NEE_at)
    
    #----------------------------------Random forests only------------------------bio4 + bio15 + 
    rf_model <- ranger(
      formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + NEE + GPP + LE + IGBP + Climate.class + source + AI"), 
      data = data_group_cv,
      mtry = m_try, 
      min.node.size = min_node_size,
      importance = "permutation"
    )
    pred_rf <- predict(rf_model, data = data_group)
    plot(pred_rf$predictions, data_group$NEE_at, main = climate_veg)
    # print(importance(rf_model))
    
    #-------------------------------combination of lmm and RF---------------------
    df_performance <- bind_rows(df_performance, data.frame(climate_veg = climate_veg, site_ID = data_group$site_ID, 
                                                           var = 'NEE_at', 
                                                           obs = data_group$NEE_at, 
                                                           rf = pred_rf$predictions, 
                                                           # lmm = pred_lmm, 
                                                           # lmm_rf = pred_lmm * (1 - weight_rf) + pred_rf$predictions * weight_rf, 
                                                           test = (1:nrow(data_group)) %in% folds[[i]]))
    #
    if (climate_veg == 'cold_deciduous_forest') {
      rf_model <- ranger(
        formula = paste0("Lgs_at ~ TA_ngs_at + SW_gs_at + P_gs_at", " + MAP + MAT + NEE + GPP + LE + IGBP + Climate.class + source + AI"), # 
        data = data_group_cv, 
        mtry = m_try, 
        min.node.size = min_node_size,
        importance = "permutation"
      )
      # rf_model$variable.importance
      pred_rf <- predict(rf_model, data = data_group)
      plot(pred_rf$predictions, data_group$Lgs_at, main = climate_veg)
      #
      df_performance <- bind_rows(df_performance, data.frame(climate_veg = climate_veg, site_ID = data_group$site_ID, 
                                                             var = 'Lgs_at', 
                                                             obs = data_group$Lgs_at, 
                                                             rf = pred_rf$predictions, 
                                                             test = (1:nrow(data_group)) %in% folds[[i]]))
    }
    #
    if (climate_veg == 'dry') {
      rf_model <- ranger(
        formula = paste0("SWC_gs_at ~ TA_gs_at + SW_gs_at + P_gs_at + LE_gs_at", " + MAP + MAT + NEE + GPP + LE + Climate.class + AI"), 
        data = data_group_cv, 
        mtry = m_try, 
        min.node.size = min_node_size,
        importance = "permutation"
      )
      # rf_model$variable.importance
      pred_rf <- predict(rf_model, data = data_group)
      plot(pred_rf$predictions, data_group$SWC_gs_at, main = climate_veg)
      #
      df_performance <- bind_rows(df_performance, data.frame(climate_veg = climate_veg, site_ID = data_group$site_ID, 
                                                             var = 'SWC_gs_at', 
                                                             obs = data_group$SWC_gs_at, 
                                                             rf = pred_rf$predictions, 
                                                             test = (1:nrow(data_group)) %in% folds[[i]]))
    }
    # 
    if (climate_veg == 'wet-warm') {
      rf_model <- ranger(
        formula = paste0("DIST_at ~ TA_at + SW_at + CO2 + P_at", " + MAP + MAT + NEE + GPP + LE + IGBP + Climate.class + source + AI"), 
        data = data_group_cv, 
        mtry = m_try, 
        min.node.size = min_node_size,
        importance = "permutation"
      )
      # rf_model$variable.importance
      pred_rf <- predict(rf_model, data = data_group)
      plot(pred_rf$predictions, data_group$DIST_at, main = climate_veg)
      #
      df_performance <- bind_rows(df_performance, data.frame(climate_veg = climate_veg, site_ID = data_group$site_ID, 
                                                             var = 'DIST_at', 
                                                             obs = data_group$DIST_at, 
                                                             rf = pred_rf$predictions, 
                                                             test = (1:nrow(data_group)) %in% folds[[i]]))
    }
  }  
}

# I can check results now. 
df_performance %>% filter(var == 'NEE_at') %>% filter(test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# climate_veg           rmse_rf    r2
# 1 cold_deciduous_forest    9.53 0.655
# 2 cold_evergreen_forest   12.6  0.447
# 3 cold_non_forest          8.75 0.412
# 4 dry                      7.90 0.493
# 5 wet-warm                13.5  0.605
 
df_performance %>% filter(var == 'NEE_at') %>% filter(!test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# climate_veg           rmse_rf    r2
# 1 cold_deciduous_forest    4.59 0.942
# 2 cold_evergreen_forest    5.62 0.935
# 3 cold_non_forest          3.91 0.918
# 4 dry                      3.77 0.932
# 5 wet-warm                 6.30 0.939

df_performance %>% filter(var == 'DIST_at') %>% filter(test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
#  wet-warm    0.00302 0.0719

df_performance %>% filter(var == 'DIST_at') %>% filter(!test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# wet-warm    0.00139 0.864

df_performance %>% filter(var == 'Lgs_at') %>% filter(test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# cold_deciduous_forest   0.819 0.372

df_performance %>% filter(var == 'Lgs_at') %>% filter(!test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# cold_deciduous_forest   0.433 0.864

df_performance %>% filter(var == 'SWC_gs_at') %>% filter(test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# dry 0.197 0.346

df_performance %>% filter(var == 'SWC_gs_at') %>% filter(!test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# dry 0.0878 0.941
