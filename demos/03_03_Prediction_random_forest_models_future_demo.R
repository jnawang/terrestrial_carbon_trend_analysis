# Develop ensemble models to project future NEE trend and use cross-validation to test model accuracy. 
# Junna Wang

rm(list=ls())

library(librarian)
shelf(tidyverse, ggpubr, lme4, dagitty, lmerTest, brms, randomForest, glmmTMB, performance, ranger, caret, Metrics, rsample, terra)

# please set your working directory
# setwd("/Users/junnawang/Downloads/demos")

trend_driver <- read.csv("potential_trend_driver.csv")
sites_long <- read.csv('long_term_ec_sites.csv')
#
sites_long <- sites_long %>% mutate(source = ifelse(source == "AmeriFlux_BASE", "Eddyproc", "Oneflux"))


# adding IGBP, climate, and category into trend driver
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "AI")], by = "site_ID")
}

#####--------
# we only use significant variables
climate_vegs <- c('cold_evergreen_forest', 'cold_deciduous_forest', 'cold_non_forest', 'dry', 'wet-warm')
key_vars <- c('CO2 + TA_gs_at',             #  + Lgs_at
              'CO2 + TA_ngs_at',            
              'SWC_gs_at',
              'CO2 + TA_gs_at + P_gs_at + SWC_gs_at', 
              'CO2 + TA_at + LE_at + DIST_at') 


# training performance; testing performance
df_performance <- data.frame()
df_partial_data <- data.frame()
df_rf_importance <- data.frame()
df_performance2 <- data.frame()       # This should focused on whether we can project past and future changes. 

# parameters of rf models
m_try = 2              # m_try = 4 is better than m_try = 3 and 5 for cold non-forest first. 
min_node_size = 5

# 5 fold cross validation
set.seed(128)
nfolder <- 5

#
for (iclimate_veg in 1:length(climate_vegs)) {
  climate_veg <- climate_vegs[iclimate_veg]
  key_var <- key_vars[iclimate_veg]
  #
  print(climate_veg)
  #
  if (climate_veg == 'cold_evergreen_forest') {
    data_group <- trend_driver %>% 
      filter(IGBP %in% c('ENF', 'EBF') & category == "cold")  
  } else if (climate_veg == 'cold_deciduous_forest') {
    data_group <- trend_driver %>% 
      filter(IGBP %in% c('DBF', 'DNF', 'MF') & category == "cold")
  } else if (climate_veg == 'cold_non_forest') {
    data_group <- trend_driver %>% 
      filter(!IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold") 
  } else if (climate_veg == 'dry') {
    data_group <- trend_driver %>% 
      filter(category == "dry") 
  } else if (climate_veg == 'wet-warm') {
    data_group <- trend_driver %>% 
      filter(category == "wet") 
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
    #----------------------------------Random forests only------------------------bio4 + bio15 + 
    rf_model <- ranger(
      formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + AI + IGBP + Climate.class + GPP_at_mean"),  #   + MAP + MAT + AI + IGBP + Climate.class + GPP_at_mean + source
      data = data_group_cv,
      mtry = m_try, 
      min.node.size = min_node_size,
      importance = "permutation", 
      seed = 189,
      num.threads = 1
    )
    pred_rf <- predict(rf_model, data = data_group, num.threads = 1)
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

      if (climate_veg == 'cold_evergreen_forest') {
        rf_model <- ranger(
          formula = paste0("Lgs_at ~ TA_ngs_at + LE_gs_at + P_gs_at + CO2", " + MAP + MAT + AI + IGBP + Climate.class"), # 
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
    if (climate_veg %in% c('cold_non_forest', "dry")) {
      rf_model <- ranger(
        formula = paste0("SWC_gs_at ~ TA_gs_at + SW_gs_at + P_gs_at + CO2", " + MAP + MAT + AI + IGBP + Climate.class"), 
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
    # Disturbance projections are not used. 
    if (climate_veg %in% c('wet-warm', "dry")) {
      if (climate_veg == "dry") {
        rf_model <- ranger(
          formula = paste0("DIST_at ~ CO2 + P_gs_at + P_ngs_at + TA_gs_at + SW_gs_at", " + MAP + MAT + AI + IGBP + Climate.class"), 
          data = data_group_cv, 
          mtry = m_try, 
          min.node.size = min_node_size,
          importance = "permutation"
        )
      } else {
        rf_model <- ranger(
          formula = paste0("DIST_at ~ CO2 + TA_at + SW_at + P_at", " + MAP + MAT + AI + IGBP + Climate.class"), 
          data = data_group_cv, 
          mtry = m_try, 
          min.node.size = min_node_size, 
          importance = "permutation"
        )
        # print(rf_model$variable.importance)
      }
      
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
  
  # A RF model for each climate-veg class
  rf_model <- ranger(
    formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + AI + IGBP + Climate.class + GPP_at_mean"), 
    data = data_group,
    mtry = m_try, 
    min.node.size = min_node_size,
    importance = "permutation"
  )
  
  pred_rf <- predict(rf_model, data = data_group)
  df_performance2 <- bind_rows(df_performance2, data.frame(climate_veg = climate_veg, site_ID = data_group$site_ID, var = "NEE_at", 
                                                           year_mid = data_group$year_mid, obs = data_group$NEE_at, rf = pred_rf$predictions))
  
  # Get overall importance of each climate-veg class
  rf_importance <- enframe(ranger::importance(rf_model), name = "variable", value = "importance") %>%
    mutate(
      relative_importance = 100 * abs(importance) / sum(abs(importance), na.rm = TRUE),
      climate_veg = climate_veg
    ) %>%
    arrange(desc(relative_importance))
  df_rf_importance <- bind_rows(df_rf_importance, rf_importance)
  
  #
  pred_fun <- function(object, newdata) {
    # predict(object, data = newdata)$predictions
    mean(predict(object, data = newdata)$predictions, na.rm = TRUE)
  }
  
  key_variable <- trimws(unlist(strsplit(key_var, "\\+")))
  partial_data <- map_dfr(key_variable, function(v) {
    pd <- pdp::partial(
      object = rf_model,
      pred.var = v,
      train = data_group,
      pred.fun = pred_fun,
      grid.resolution = 50,
      ice = FALSE,
      progress = "none"
    ) %>%
      as_tibble()
    
    pd %>%
      mutate(
        variable = v,
        x_value = as.character(.data[[v]]),
        x_numeric = suppressWarnings(as.numeric(as.character(.data[[v]]))),
        climate_veg = climate_veg
      ) %>%
      dplyr::select(climate_veg, variable, x_value, x_numeric, yhat)
  })
  
  df_partial_data <- bind_rows(df_partial_data, partial_data)
}

# please note that: your result may be different from mine, as different random seeds were used. 
df_performance %>% filter(var == 'NEE_at') %>% filter(test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2, mae = mean(abs(obs-rf)))
# 1 cold_deciduous_forest    9.69 0.659  7.33
# 2 cold_evergreen_forest   10.3  0.712  7.20
# 3 cold_non_forest          8.75 0.404  4.49
# 4 dry                      7.35 0.531  4.84
# 5 wet-warm                13.9  0.436  8.90
 
df_performance %>% filter(var == 'NEE_at') %>% filter(!test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# climate_veg           rmse_rf    r2
# 1 cold_deciduous_forest    5.43 0.924
# 2 cold_evergreen_forest    4.77 0.941
# 3 cold_non_forest          4.13 0.893
# 4 dry                      3.48 0.924
# 5 wet-warm                 6.48 0.926


df_performance %>% filter(var == 'SWC_gs_at') %>% filter(test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# 1 cold_non_forest   0.111 0.483
# 2 dry               0.129 0.697

df_performance %>% filter(var == 'SWC_gs_at') %>% filter(!test) %>% group_by(climate_veg) %>% 
  summarise(rmse_rf = sqrt(mean((obs-rf)^2)), r2 = cor(obs, rf)^2)
# 1 cold_non_forest  0.0461 0.937
# 2 dry              0.0461 0.964

