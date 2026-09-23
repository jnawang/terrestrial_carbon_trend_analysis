# Weighted LM, RF, CO2-linear + RF hybrid, and an equal-weight LM-RF ensemble
# Assessments: (1) LOSO temporal change in NEE trend; (2) five-fold direct prediction
# what I learned, the way of giving random seeds affect RF performance. 

rm(list=ls())

library(librarian)
shelf(tidyverse,ranger,caret,ggpubr)

trend_driver <- read.csv(file.path("data2026","potential_trend_driver.csv"))
sites_long <- read.csv(file.path("data2026","long_term_ec_sites.csv")) %>%
  mutate(source=ifelse(source=="AmeriFlux_BASE","Eddyproc","Oneflux"))

if(!"category" %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(
    sites_long[,c("site_ID","IGBP","Climate.class","category","NEEpartition","source","AI")],
    by="site_ID"
  )
}

trend_driver_site <- trend_driver %>% group_by(site_ID) %>% summarise(
  NEE_at_sd_min=if(all(is.na(NEE_at_sd))) NA_real_ else min(NEE_at_sd,na.rm=TRUE),
  .groups="drop"
)
trend_driver <- trend_driver %>% left_join(trend_driver_site,by="site_ID") %>% mutate(
  NEE_at_sd_scaled=NEE_at_sd/(NEE_at_sd_min*sqrt(nyear)),
  model_weight=1/NEE_at_sd_scaled^2
) %>% filter(is.finite(model_weight),model_weight>0)
trend_driver$.obs_id <- seq_len(nrow(trend_driver))
variable_weight <- "model_weight"
climate_vegs <- c("cold_evergreen_forest","cold_deciduous_forest","cold_non_forest","dry","wet-warm")
key_vars <- c(
  "CO2 + TA_gs_at + SWC_gs_at",
  "CO2 + TA_gs_at + DIST_at",
  "CO2 + SWC_gs_at",
  "CO2 + TA_gs_at + P_gs_at + SWC_gs_at",
  "CO2 + TA_at + LE_at"
)
nee_common_predictors <- c("MAP","MAT","IGBP","GPP_at_mean","NEE")
swc_predictors <- c("CO2","TA_gs_at","SW_gs_at","P_gs_at","MAP","MAT","NEE","IGBP")
swc_groups <- c("cold_evergreen_forest","cold_non_forest","dry")
base_methods <- c("LM","RF","LM + RF")
method_order <- c("Ensemble","LM","RF","LM + RF")
m_try <- 2
min_node_size <- 3
num_trees <- 500
sample_fraction <- 0.8
replace_samples <- TRUE
split_rule <- "extratrees"

window_threshold <- 5

get_nee_predictors <- function(key_var) unique(c(all.vars(as.formula(paste("~",key_var))),nee_common_predictors))

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
    if(is.factor(train[[v]])) return(all(as.character(test[[v]]) %in% as.character(train[[v]])))
    TRUE
  },logical(1))
  predictors[keep]
}

fit_co2 <- function(y,data) {
  d <- data
  d$.y <- y
  lm(.y~CO2,data=d,weights=.model_weight)
}

fit_hybrid <- function(train,test,rf_predictors,response="NEE_at",mtry=4,min.node.size=1,
                       num.trees=500,max_iter=30,tol=2e-3,seed=128) {
  stopifnot(!"CO2" %in% rf_predictors,length(rf_predictors)>0)
  y <- train[[response]]
  w <- train$.model_weight
  co2_model <- fit_co2(y,train)
  co2_train <- as.numeric(predict(co2_model,newdata=train))
  y_mean <- weighted.mean(y,w)
  y_scale <- max(sqrt(weighted.mean((y-y_mean)^2,w)),.Machine$double.eps)
  converged <- FALSE

  for(iteration in seq_len(max_iter)) {
    rf_data <- train
    rf_data$.resid <- y-co2_train
    rf_model <- ranger(
      reformulate(rf_predictors,response=".resid"),data=rf_data,num.trees=num.trees,
      mtry=min(mtry,length(rf_predictors)),min.node.size=min.node.size,
      sample.fraction=sample_fraction,replace=replace_samples, splitrule=split_rule,
      case.weights=rf_data$.model_weight,oob.error=TRUE,seed=seed+iteration
    )
    rf_oob_raw <- rf_model$predictions
    if(anyNA(rf_oob_raw)) stop("Missing OOB predictions in LM + RF.")
    rf_center <- weighted.mean(rf_oob_raw,w)
    rf_oob <- rf_oob_raw-rf_center
    new_co2_model <- fit_co2(y-rf_oob,train)
    new_co2_train <- as.numeric(predict(new_co2_model,newdata=train))
    change <- sqrt(weighted.mean((new_co2_train-co2_train)^2,w))/y_scale
    co2_model <- new_co2_model
    co2_train <- new_co2_train
    if(change<tol) {
      converged <- TRUE
      break
    }
  }
  if(!converged) warning("LM + RF did not converge after ",max_iter," iterations.")

  rf_data <- train
  rf_data$.resid <- y-co2_train
  rf_model <- ranger(
    reformulate(rf_predictors,response=".resid"),data=rf_data,num.trees=num.trees,
    mtry=min(mtry,length(rf_predictors)),min.node.size=min.node.size,
    sample.fraction=sample_fraction,replace=replace_samples,splitrule=split_rule,
    case.weights=rf_data$.model_weight,oob.error=TRUE,seed=seed+max_iter+1L
  )
  rf_center <- weighted.mean(rf_model$predictions,w)

  as.numeric(predict(co2_model,newdata=test))+
    as.numeric(predict(rf_model,data=test)$predictions)-rf_center
}

fit_three_models <- function(train,test,response,predictors,seed) {
  train$.model_weight <- train[[variable_weight]]/mean(train[[variable_weight]])
  p <- usable_predictors(train,test,predictors)
  if(!"CO2" %in% p) stop("CO2 is not estimable for response ",response,".")
  factor_vars <- p[vapply(train[p],is.factor,logical(1))]
  for(v in factor_vars) test[[v]] <- factor(as.character(test[[v]]),levels=levels(train[[v]]))
  p_no_co2 <- setdiff(p,"CO2")

  rf_model <- ranger(
    reformulate(p,response=response),data=train,num.trees=num_trees,
    mtry=min(m_try,length(p)),min.node.size=min_node_size,
    sample.fraction=sample_fraction,replace=replace_samples,splitrule=split_rule,
    case.weights=train$.model_weight,seed=seed
  )
  lm_model <- lm(reformulate(p,response=response),data=train,weights=.model_weight)

  list(
    LM=as.numeric(predict(lm_model,newdata=test)),
    RF=as.numeric(predict(rf_model,data=test)$predictions),
    `LM + RF`=fit_hybrid(
      train,test,p_no_co2,response=response,mtry=m_try,min.node.size=min_node_size,
      num.trees=num_trees,seed=seed+1000L
    )
  )
}

add_ensemble <- function(df) {
  ensemble <- df %>% filter(method %in% c("LM","RF")) %>%
    group_by(obs_id,fold,climate_veg,site_ID,var,year_mid,obs) %>%
    summarise(pred=mean(pred),.groups="drop") %>% mutate(method="Ensemble")
  bind_rows(df,ensemble)
}

change_per_decade <- function(year,y) {
  ok <- is.finite(year) & is.finite(y)
  if(n_distinct(year[ok])<window_threshold) return(NA_real_)
  10*unname(coef(lm(y[ok]~year[ok]))[2])
}

cor_test_p <- function(x,y) {
  ok <- is.finite(x) & is.finite(y)
  if(sum(ok)<3 || sd(x[ok])==0 || sd(y[ok])==0) return(NA_real_)
  cor.test(x[ok],y[ok])$p.value
}

# =============================================================================
# 1. Leave-one-site-out prediction and temporal-change assessment
# =============================================================================

df_performance_loso_base <- data.frame()

for(i in seq_along(climate_vegs)) {
  climate_veg <- climate_vegs[i]
  predictors <- get_nee_predictors(key_vars[i])
  data_group <- get_group(trend_driver,climate_veg)
  needed <- c("NEE_at","site_ID","year_mid",variable_weight,predictors)
  data_group <- data_group[complete.cases(data_group[,needed,drop=FALSE]),,drop=FALSE]
  data_group <- data_group[is.finite(data_group[[variable_weight]]) & data_group[[variable_weight]]>0,,drop=FALSE]
  data_group$IGBP <- factor(data_group$IGBP)

  test_sites <- data_group %>% group_by(site_ID) %>%
    summarise(nwindow=n_distinct(year_mid),.groups="drop") %>%
    filter(nwindow>=window_threshold) %>% pull(site_ID)

  for(site_index in seq_along(test_sites)) {
    site <- test_sites[[site_index]]
    message("LOSO ",climate_veg,": ",site)
    
    test_id <- data_group$site_ID==site
    train <- droplevels(data_group[!test_id,,drop=FALSE])
    test <- data_group[test_id,,drop=FALSE]
    
    rf_seed <- 100000L+i*100L+site_index
    predictions <- fit_three_models(
      train,test,"NEE_at",predictors,rf_seed
    )
    
    df_performance_loso_base <- bind_rows(
      df_performance_loso_base,
      bind_rows(lapply(names(predictions),function(method_name) data.frame(
        obs_id=test$.obs_id,fold=NA_integer_,method=method_name,
        climate_veg=climate_veg,site_ID=test$site_ID,var="NEE_at",
        year_mid=test$year_mid,obs=test$NEE_at,
        pred=as.numeric(predictions[[method_name]])
      )))
    )
  }
}

df_performance_loso <- add_ensemble(df_performance_loso_base)
df_performance2 <- df_performance_loso

df_change_loso <- df_performance_loso %>%
  filter(is.finite(year_mid),is.finite(obs),is.finite(pred)) %>%
  group_by(method,climate_veg,site_ID) %>% filter(n_distinct(year_mid)>=window_threshold) %>% summarise(
    n_year=n_distinct(year_mid),year_start=min(year_mid),year_end=max(year_mid),
    obs_change=change_per_decade(year_mid,obs),
    pred_change=change_per_decade(year_mid,pred),.groups="drop"
  ) 
df_change_cv <- df_change_loso

df_metrics_loso <- df_change_loso %>% filter(is.finite(obs_change),is.finite(pred_change)) %>%
  group_by(method) %>% summarise(
    nsite=n(),RMSE=sqrt(mean((obs_change-pred_change)^2)),
    r=cor(obs_change,pred_change),R2=r^2,p=cor_test_p(obs_change,pred_change),
    MAE=mean(abs(obs_change-pred_change)),bias=mean(pred_change-obs_change),.groups="drop"
  ) %>% mutate(p_adj=p.adjust(p,method="BH"))
df_metrics_cv <- df_metrics_loso
print(df_metrics_loso)

plot_loso <- df_change_loso %>% ggplot(aes(pred_change,obs_change,color=method))+
  geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~method,scales="free")+
  theme_classic()+labs(x="Predicted change per decade",y="Observed change per decade")
print(plot_loso)

# =============================================================================
# 2. Five-fold row-wise CV for direct NEE_at and SWC_gs_at predictions
# =============================================================================

df_performance_5fold_base <- data.frame()
nfolder <- 5

for(i in seq_along(climate_vegs)) {
  climate_veg <- climate_vegs[i]
  data_climate <- get_group(trend_driver,climate_veg)
  target_specs <- list(NEE_at=get_nee_predictors(key_vars[i]))
  if(climate_veg %in% swc_groups) target_specs$SWC_gs_at <- swc_predictors

  for(response in names(target_specs)) {
    predictors <- unique(target_specs[[response]])
    needed <- c(response,"site_ID","year_mid",variable_weight,predictors)
    data_group <- data_climate[complete.cases(data_climate[,needed,drop=FALSE]),,drop=FALSE]
    data_group <- data_group[is.finite(data_group[[variable_weight]]) & data_group[[variable_weight]]>0,,drop=FALSE]
    data_group$IGBP <- factor(data_group$IGBP)
    if(nrow(data_group)<nfolder) next

    set.seed(5000L+100L*i+match(response,names(target_specs)))
    folds <- caret::createFolds(data_group[[response]],k=nfolder,list=TRUE,returnTrain=FALSE)

    for(fold_id in seq_along(folds)) {
      test_index <- folds[[fold_id]]
      train <- droplevels(data_group[-test_index,,drop=FALSE])
      test <- data_group[test_index,,drop=FALSE]
      predictions <- fit_three_models(train,test,response,predictors,6000L+100L*i+fold_id)

      df_performance_5fold_base <- bind_rows(
        df_performance_5fold_base,
        bind_rows(lapply(names(predictions),function(method_name) data.frame(
          obs_id=test$.obs_id,fold=fold_id,method=method_name,climate_veg=climate_veg,
          site_ID=test$site_ID,var=response,year_mid=test$year_mid,obs=test[[response]],
          pred=as.numeric(predictions[[method_name]])
        )))
      )
    }
  }
}

df_performance_5fold <- add_ensemble(df_performance_5fold_base)

df_metrics_5fold <- df_performance_5fold %>% filter(is.finite(obs),is.finite(pred)) %>%
  group_by(var,climate_veg,method) %>% summarise(
    nobs=n(),RMSE=sqrt(mean((obs-pred)^2)),r=cor(obs,pred),R2=r^2,
    p=cor_test_p(obs,pred),MAE=mean(abs(obs-pred)),bias=mean(pred-obs),.groups="drop"
  ) %>% group_by(var,climate_veg) %>% mutate(p_adj=p.adjust(p,method="BH")) %>% ungroup()

df_metrics_5fold_overall <- df_performance_5fold %>% filter(is.finite(obs),is.finite(pred)) %>%
  group_by(var,method) %>% summarise(
    nobs=n(),RMSE=sqrt(mean((obs-pred)^2)),r=cor(obs,pred),R2=r^2,
    p=cor_test_p(obs,pred),MAE=mean(abs(obs-pred)),bias=mean(pred-obs),.groups="drop"
  ) %>% group_by(var) %>% mutate(p_adj=p.adjust(p,method="BH")) %>% ungroup()

print(df_metrics_5fold)
print(df_metrics_5fold_overall)

climate_order <- c("cold_evergreen_forest","cold_deciduous_forest","cold_non_forest","wet-warm","dry")

table_5fold_publication <- df_metrics_5fold %>%
  filter(var=="NEE_at",method %in% method_order) %>% mutate(
    climate_veg=factor(climate_veg,levels=climate_order),
    method=factor(method,levels=method_order),value=sprintf("%.3f (%.2f)",R2,MAE)
  ) %>% select(climate_veg,method,value) %>% arrange(climate_veg,method) %>%
  pivot_wider(names_from=method,values_from=value) %>%
  select(climate_veg,all_of(method_order)) %>% arrange(climate_veg)
print(table_5fold_publication)

table_5fold_overall <- df_metrics_5fold_overall %>% filter(method %in% method_order) %>% mutate(
  var=factor(var,levels=c("NEE_at","SWC_gs_at")),method=factor(method,levels=method_order),
  value=if_else(var=="SWC_gs_at",sprintf("%.3f (%.3f)",r,MAE),sprintf("%.3f (%.2f)",R2,MAE))
) %>% select(var,method,value) %>% arrange(var,method) %>%
  pivot_wider(names_from=method,values_from=value) %>% select(var,all_of(method_order)) %>%
  mutate(var=as.character(var))
print(table_5fold_overall)

table_cv_publication <- df_metrics_loso %>% filter(method %in% method_order) %>% mutate(
  method=factor(method,levels=method_order),value=sprintf("%.3f (%.2f)",R2,MAE)
) %>% select(method,value) %>% arrange(method) %>% pivot_wider(names_from=method,values_from=value) %>%
  select(all_of(method_order)) %>% mutate(assessment="LOSO NEE change per decade",.before=1)
print(table_cv_publication)
