# Tune a weighted RF using leave-one-site-out performance of NEE-trend change.
# Model fitting is weighted; site-level decadal trends and metrics are unweighted.

rm(list=ls())

library(librarian)
shelf(tidyverse,ranger)

set.seed(128)
data_dir <- "data2026"
outdir <- file.path("output","RF_LOSO_tuning")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

window_threshold <- 5
num_trees_tuning <- 500
num_threads <- 1

# 4 x 4 x 3 x 2 = 96 parameter combinations. mtry = 1, min.node.size = 10, fraction = 0.8; give the best LOSO result!
tuning_grid <- tidyr::expand_grid(
  mtry=c(1L,2L,3L),
  min.node.size=c(2L,3L,4L,5L,6L,7L,8L,9L,10L,12L),
  sample.fraction=c(0.8,1),
  splitrule=c("variance","extratrees")
) %>% mutate(parameter_id=row_number(),.before=1)

trend_driver <- read.csv(file.path(data_dir,"potential_trend_driver.csv"))
sites_long <- read.csv(file.path(data_dir,"long_term_ec_sites.csv")) %>%
  mutate(source=ifelse(source=="AmeriFlux_BASE","Eddyproc","Oneflux"))

if(!"category" %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(
    sites_long[,c("site_ID","IGBP","Climate.class","category","NEEpartition","source","AI")],
    by="site_ID"
  )
}

# Same uncertainty weights used in the current 03_03 models.
site_uncertainty <- trend_driver %>% group_by(site_ID) %>% summarise(
  NEE_at_sd_min=if(all(is.na(NEE_at_sd))) NA_real_ else min(NEE_at_sd,na.rm=TRUE),
  .groups="drop"
)

trend_driver <- trend_driver %>% left_join(site_uncertainty,by="site_ID") %>% mutate(
  NEE_at_sd_scaled=NEE_at_sd/(NEE_at_sd_min*sqrt(nyear)),
  model_weight=1/NEE_at_sd_scaled^2,
  .obs_id=row_number()
) %>% filter(is.finite(model_weight),model_weight>0)

climate_vegs <- c(
  "cold_evergreen_forest","cold_deciduous_forest",
  "cold_non_forest","dry","wet-warm"
)

key_vars <- c(
  "CO2 + TA_gs_at + SWC_gs_at",
  "CO2 + TA_gs_at + DIST_at",
  "CO2 + SWC_gs_at",
  "CO2 + TA_gs_at + P_gs_at + SWC_gs_at",
  "CO2 + TA_at + LE_at"
)

nee_common_predictors <- c("MAP","MAT","IGBP","GPP_at_mean","NEE")
get_nee_predictors <- function(key_var) unique(c(
  all.vars(as.formula(paste("~",key_var))),nee_common_predictors
))

get_group <- function(d,g) {
  if(g=="cold_evergreen_forest") filter(d,IGBP %in% c("ENF","EBF"),category=="cold")
  else if(g=="cold_deciduous_forest") filter(d,IGBP %in% c("DBF","DNF","MF"),category=="cold")
  else if(g=="cold_non_forest") filter(d,!IGBP %in% c("ENF","EBF","DBF","DNF","MF"),category=="cold")
  else if(g=="dry") filter(d,category=="dry")
  else filter(d,category=="wet")
}

usable_predictors <- function(train,test,predictors) {
  keep <- vapply(predictors,function(v) {
    if(dplyr::n_distinct(train[[v]])<2) return(FALSE)
    if(is.factor(train[[v]]))
      return(all(as.character(test[[v]]) %in% as.character(train[[v]])))
    TRUE
  },logical(1))
  predictors[keep]
}

change_per_decade <- function(year,y) {
  ok <- is.finite(year) & is.finite(y)
  if(dplyr::n_distinct(year[ok])<window_threshold) return(NA_real_)
  10*unname(coef(lm(y[ok]~year[ok]))[2])
}

cor_test_p <- function(x,y) {
  ok <- is.finite(x) & is.finite(y)
  if(sum(ok)<3 || sd(x[ok])==0 || sd(y[ok])==0) return(NA_real_)
  cor.test(x[ok],y[ok])$p.value
}

# Prepare the identical LOSO datasets once, before tuning.
analysis_groups <- purrr::map2(climate_vegs,key_vars,function(climate_veg,key_var) {
  predictors <- get_nee_predictors(key_var)
  d <- get_group(trend_driver,climate_veg)
  needed <- c("NEE_at","site_ID","year_mid","model_weight",predictors)
  d <- d[complete.cases(d[,needed,drop=FALSE]),,drop=FALSE]
  d <- d[is.finite(d$model_weight) & d$model_weight>0,,drop=FALSE]
  d$IGBP <- factor(d$IGBP)
  test_sites <- d %>% group_by(site_ID) %>%
    summarise(nwindow=n_distinct(year_mid),.groups="drop") %>%
    filter(nwindow>=window_threshold) %>% pull(site_ID)
  list(climate_veg=climate_veg,predictors=predictors,data=d,test_sites=test_sites)
})

fit_one_loso <- function(group_object,parameter,group_index) {
  d <- group_object$data
  predictors <- group_object$predictors
  climate_veg <- group_object$climate_veg

  purrr::imap_dfr(group_object$test_sites,function(site,site_index) {
    test_id <- d$site_ID==site
    train <- droplevels(d[!test_id,,drop=FALSE])
    test <- d[test_id,,drop=FALSE]
    train$.model_weight <- train$model_weight/mean(train$model_weight)

    p <- usable_predictors(train,test,predictors)
    factor_vars <- p[vapply(train[p],is.factor,logical(1))]
    for(v in factor_vars)
      test[[v]] <- factor(as.character(test[[v]]),levels=levels(train[[v]]))

    result <- tryCatch({
      rf_model <- ranger::ranger(
        reformulate(p,response="NEE_at"),data=train,
        num.trees=num_trees_tuning,
        mtry=min(parameter$mtry,length(p)),
        min.node.size=parameter$min.node.size,
        sample.fraction=parameter$sample.fraction,
        replace=TRUE,
        splitrule=parameter$splitrule,
        case.weights=train$.model_weight,
        importance="none",
        # seed=100000L+parameter$parameter_id*1000L+group_index*100L+site_index,
        seed=100000L+group_index*100L+site_index, 
        num.threads=num_threads
      )
      pred <- as.numeric(predict(rf_model,data=test)$predictions)
      tibble(
        obs_change=change_per_decade(test$year_mid,test$NEE_at),
        pred_change=change_per_decade(test$year_mid,pred),
        error=NA_character_
      )
    },error=function(e) tibble(
      obs_change=change_per_decade(test$year_mid,test$NEE_at),
      pred_change=NA_real_,error=conditionMessage(e)
    ))

    result %>% mutate(
      parameter_id=parameter$parameter_id,
      mtry=parameter$mtry,
      min.node.size=parameter$min.node.size,
      sample.fraction=parameter$sample.fraction,
      splitrule=parameter$splitrule,
      climate_veg=climate_veg,site_ID=site,
      .before=1
    )
  })
}

site_results <- vector("list",nrow(tuning_grid))
checkpoint_file <- file.path(outdir,"RF_LOSO_tuning_site_changes_checkpoint.csv")

for(j in seq_len(nrow(tuning_grid))) {
  parameter <- tuning_grid[j,,drop=FALSE]
  message(
    "Combination ",j,"/",nrow(tuning_grid),": mtry=",parameter$mtry,
    ", min.node.size=",parameter$min.node.size,
    ", sample.fraction=",parameter$sample.fraction,
    ", splitrule=",parameter$splitrule
  )
  site_results[[j]] <- purrr::imap_dfr(
    analysis_groups,
    ~fit_one_loso(.x,parameter,.y)
  )
  write.csv(bind_rows(site_results[seq_len(j)]),checkpoint_file,row.names=FALSE)
}

df_rf_loso_tuning_sites <- bind_rows(site_results)

summarise_performance <- function(d) {
  valid <- is.finite(d$obs_change) & is.finite(d$pred_change)
  tibble(
    nsite=sum(valid),
    n_failed=sum(!valid),
    RMSE=sqrt(mean((d$obs_change[valid]-d$pred_change[valid])^2)),
    r=cor(d$obs_change[valid],d$pred_change[valid]),
    R2=r^2,
    p=cor_test_p(d$obs_change[valid],d$pred_change[valid]),
    MAE=mean(abs(d$obs_change[valid]-d$pred_change[valid])),
    bias=mean(d$pred_change[valid]-d$obs_change[valid])
  )
}

df_rf_loso_tuning_overall <- df_rf_loso_tuning_sites %>%
  group_by(parameter_id,mtry,min.node.size,sample.fraction,splitrule) %>%
  group_modify(~summarise_performance(.x)) %>% ungroup() %>%
  mutate(
    eligible=nsite==max(nsite),
    rank_RMSE=if_else(eligible,rank(RMSE,ties.method="min"),NA_real_),
    rank_MAE=if_else(eligible,rank(MAE,ties.method="min"),NA_real_),
    rank_R2=if_else(eligible,rank(-R2,ties.method="min"),NA_real_)
  ) %>% arrange(desc(eligible),RMSE,MAE,desc(R2))

df_rf_loso_tuning_by_climate <- df_rf_loso_tuning_sites %>%
  group_by(parameter_id,mtry,min.node.size,sample.fraction,splitrule,climate_veg) %>%
  group_modify(~summarise_performance(.x)) %>% ungroup() %>%
  group_by(climate_veg) %>% mutate(
    eligible=nsite==max(nsite),
    rank_RMSE=if_else(eligible,rank(RMSE,ties.method="min"),NA_real_),
    rank_MAE=if_else(eligible,rank(MAE,ties.method="min"),NA_real_),
    rank_R2=if_else(eligible,rank(-R2,ties.method="min"),NA_real_)
  ) %>% ungroup() %>% arrange(climate_veg,desc(eligible),RMSE,MAE,desc(R2))

eligible_overall <- df_rf_loso_tuning_overall %>% filter(eligible)
best_parameters <- bind_rows(
  eligible_overall %>% slice_min(RMSE,n=1,with_ties=FALSE) %>% mutate(selected_by="minimum RMSE",.before=1),
  eligible_overall %>% slice_min(MAE,n=1,with_ties=FALSE) %>% mutate(selected_by="minimum MAE",.before=1),
  eligible_overall %>% slice_max(R2,n=1,with_ties=FALSE) %>% mutate(selected_by="maximum R2",.before=1)
)

best_parameters_by_climate <- df_rf_loso_tuning_by_climate %>% filter(eligible) %>%
  group_by(climate_veg) %>% slice_min(RMSE,n=1,with_ties=FALSE) %>% ungroup() %>%
  mutate(selected_by="minimum RMSE",.before=1)

write.csv(df_rf_loso_tuning_sites,file.path(outdir,"RF_LOSO_tuning_site_changes.csv"),row.names=FALSE)
write.csv(df_rf_loso_tuning_overall,file.path(outdir,"RF_LOSO_tuning_metrics_overall.csv"),row.names=FALSE)
write.csv(df_rf_loso_tuning_by_climate,file.path(outdir,"RF_LOSO_tuning_metrics_by_climate.csv"),row.names=FALSE)
write.csv(best_parameters,file.path(outdir,"RF_LOSO_best_parameters_overall.csv"),row.names=FALSE)
write.csv(best_parameters_by_climate,file.path(outdir,"RF_LOSO_best_parameters_by_climate.csv"),row.names=FALSE)

cat("\nBest parameter combinations across all 52 sites:\n")
print(best_parameters)
cat("\nTop 20 combinations ranked by RMSE:\n")
print(df_rf_loso_tuning_overall %>% filter(eligible) %>% slice_head(n=20))
message("Finished. Tuning results are in: ",outdir)
