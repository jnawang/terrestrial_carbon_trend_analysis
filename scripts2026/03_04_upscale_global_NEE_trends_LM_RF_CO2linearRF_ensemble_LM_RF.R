# Global projection using weighted LM, weighted RF, linear-CO2 + RF, and LM-RF ensemble.
# The ensemble is strictly (LM + RF)/2; the hybrid is projected independently.
# SWC_gs_at is projected once with weighted RF and used by all NEE models.

rm(list=ls())
library(librarian)
shelf(tidyverse,terra,trend,ranger)
source(file.path("scripts2026","00_functions.R"))
set.seed(128)

data_dir <- "data2026"
output_dir <- data_dir
dist_at_past_file <- file.path(data_dir,"LAI_deficit_trends_DBF.tif")
dist_at_future_file <- dist_at_past_file

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
method_ids <- c("LM","RF","CO2_RF")
method_labels <- c(LM="LM",RF="RF",CO2_RF="Linear CO2 + RF")

m_try <- 2
min_node_size <- 3          # the result similar even if we use min_node_size = 5 here; but 10 coresponds to best LOSO validation.  
num_trees <- 500
sample_fraction <- 0.8
replace_samples <- TRUE
split_rule <- "extratrees"

hybrid_max_iter <- 30
hybrid_tol <- 2e-3

get_nee_predictors <- function(key_var) unique(c(all.vars(as.formula(paste("~",key_var))),nee_common_predictors))

get_group <- function(d,group) {
  if(group=="cold_evergreen_forest") {
    filter(d,IGBP %in% c("ENF","EBF"),category=="cold")
  } else if(group=="cold_deciduous_forest") {
    filter(d,IGBP %in% c("DBF","DNF","MF"),category=="cold")
  } else if(group=="cold_non_forest") {
    filter(d,!IGBP %in% c("ENF","EBF","DBF","DNF","MF"),category=="cold")
  } else if(group=="dry") {
    filter(d,category=="dry")
  } else {
    filter(d,category=="wet")
  }
}

prepare_training_data <- function(data,response,predictors) {
  needed <- c(response,predictors,"model_weight")
  out <- data[complete.cases(data[,needed,drop=FALSE]),,drop=FALSE]
  out <- out[is.finite(out$model_weight) & out$model_weight>0,,drop=FALSE]
  out$IGBP <- factor(out$IGBP)
  out <- droplevels(out)
  predictors <- predictors[vapply(predictors,function(v) dplyr::n_distinct(out[[v]])>1,logical(1))]
  if(!"CO2" %in% predictors) stop("CO2 is not estimable for ",response,".")
  if(nrow(out)<10) stop("Too few complete training rows for ",response,".")
  out$.model_weight <- out$model_weight/mean(out$model_weight)
  list(data=out,predictors=predictors)
}

align_factor_levels <- function(train,newdata,predictors) {
  out <- newdata
  factor_vars <- predictors[vapply(train[predictors],is.factor,logical(1))]
  for(v in factor_vars) {
    lev <- levels(train[[v]])
    fallback <- names(which.max(table(train[[v]])))
    x <- as.character(out[[v]])
    x[!is.na(x) & !x %in% lev] <- fallback
    out[[v]] <- factor(x,levels=lev)
  }
  out
}

fit_linear_co2 <- function(y,data) {
  d <- data
  d$.response <- y
  lm(.response~CO2,data=d,weights=.model_weight)
}

fit_rf_component <- function(data,response,predictors,seed) {
  ranger::ranger(
    reformulate(predictors,response=response),data=data,
    num.trees=num_trees,
    mtry=min(m_try,length(predictors)),
    min.node.size=min_node_size,
    sample.fraction=sample_fraction,
    replace=replace_samples,
    splitrule=split_rule,
    case.weights=data$.model_weight,
    importance="permutation",
    oob.error=TRUE,
    seed=seed
  )
}

fit_co2_rf_model <- function(train,response,predictors,seed) {
  rf_predictors <- setdiff(predictors,"CO2")
  if(!length(rf_predictors)) stop("Linear CO2 + RF requires non-CO2 predictors.")
  y <- train[[response]]
  w <- train$.model_weight
  co2_model <- fit_linear_co2(y,train)
  co2_train <- as.numeric(predict(co2_model,newdata=train))
  y_mean <- weighted.mean(y,w)
  y_scale <- max(sqrt(weighted.mean((y-y_mean)^2,w)),.Machine$double.eps)
  converged <- FALSE

  for(iteration in seq_len(hybrid_max_iter)) {
    rf_data <- train
    rf_data$.resid <- y-co2_train
    rf_model <- fit_rf_component(rf_data,".resid",rf_predictors,seed+iteration)
    rf_oob_raw <- rf_model$predictions
    if(anyNA(rf_oob_raw)) stop("Missing RF OOB predictions in linear CO2 + RF.")
    rf_center <- weighted.mean(rf_oob_raw,w)
    rf_oob <- rf_oob_raw-rf_center
    new_co2_model <- fit_linear_co2(y-rf_oob,train)
    new_co2_train <- as.numeric(predict(new_co2_model,newdata=train))
    relative_change <- sqrt(weighted.mean((new_co2_train-co2_train)^2,w))/y_scale
    co2_model <- new_co2_model
    co2_train <- new_co2_train
    if(relative_change<hybrid_tol) {
      converged <- TRUE
      break
    }
  }
  if(!converged) warning("Linear CO2 + RF for ",response," did not converge after ",hybrid_max_iter," iterations.")

  rf_data <- train
  rf_data$.resid <- y-co2_train
  rf_model <- fit_rf_component(rf_data,".resid",rf_predictors,seed+hybrid_max_iter+1L)
  if(anyNA(rf_model$predictions)) stop("Missing final RF OOB predictions in linear CO2 + RF.")
  rf_center <- weighted.mean(rf_model$predictions,w)

  structure(list(
    co2_model=co2_model,rf_model=rf_model,rf_center=rf_center,
    train=train,predictors=predictors,rf_predictors=rf_predictors,
    response=response,converged=converged,iterations=iteration
  ),class="co2_rf_model")
}

fit_models <- function(data,response,predictors,seed) {
  prepared <- prepare_training_data(data,response,predictors)
  train <- prepared$data
  predictors <- prepared$predictors
  list(
    LM=list(
      model=lm(reformulate(predictors,response=response),data=train,weights=.model_weight),
      train=train,predictors=predictors
    ),
    RF=list(
      model=fit_rf_component(train,response,predictors,seed),
      train=train,predictors=predictors
    ),
    CO2_RF=fit_co2_rf_model(train,response,predictors,seed+100L)
  )
}

predict_one_model <- function(object,newdata,method) {
  out <- rep(NA_real_,nrow(newdata))
  d <- align_factor_levels(object$train,newdata,object$predictors)
  ok <- complete.cases(d[,object$predictors,drop=FALSE])
  if(!any(ok)) return(out)
  if(method=="LM") {
    out[ok] <- as.numeric(predict(object$model,newdata=d[ok,,drop=FALSE]))
  } else if(method=="RF") {
    out[ok] <- as.numeric(predict(object$model,data=d[ok,,drop=FALSE])$predictions)
  } else if(method=="CO2_RF") {
    out[ok] <- as.numeric(predict(object$co2_model,newdata=d[ok,,drop=FALSE]))+
      as.numeric(predict(object$rf_model,data=d[ok,,drop=FALSE])$predictions)-object$rf_center
  }
  out
}

calculate_pdp <- function(models,climate_veg,response,plot_predictors,n_grid=50L) {
  n_grid <- as.integer(n_grid)[1]
  if(!is.finite(n_grid) || n_grid<2) stop("n_grid must be greater than 1.")
  purrr::imap_dfr(models,function(object,method) {
    train <- object$train
    weight <- as.numeric(train$.model_weight)
    predictors <- intersect(plot_predictors,object$predictors)
    baseline <- weighted.mean(predict_one_model(object,train,method),weight,na.rm=TRUE)
    purrr::map_dfr(predictors,function(variable) {
      x <- as.numeric(train[[variable]])
      limits <- as.numeric(stats::quantile(x,c(.02,.98),na.rm=TRUE,names=FALSE))
      if(!all(is.finite(limits))) stop("Non-finite range for ",variable,".")
      grid <- seq.int(limits[1],limits[2],length.out=n_grid)
      purrr::map_dfr(grid,function(value) {
        newdata <- train
        newdata[[variable]] <- value
        prediction <- predict_one_model(object,newdata,method)
        pdp_prediction <- weighted.mean(prediction,weight,na.rm=TRUE)
        tibble(
          climate_veg=climate_veg,response=response,method=method_labels[[method]],
          predictor=variable,value=value,prediction=pdp_prediction,
          partial_effect=pdp_prediction-baseline
        )
      })
    })
  })
}

# Raster preparation
nee_baseline <- rast(file.path(data_dir,"baseline_nee_2005-2021.tif"))
mat_baseline <- rast(file.path(data_dir,"baseline_mat.tif"))
map_baseline <- rast(file.path(data_dir,"baseline_map.tif"))
gpp_at_baseline <- rast(file.path(data_dir,"CEDAR_GPP_Sen_slope_ML.tif"))
land_cover <- rast(file.path(data_dir,"968_Land_Cover_Class_0.25degree_layer2.tif"))
zone_group <- rast(file.path(data_dir,"zones_climate_vegetation_group.tif"))
co2_past <- rast(file.path(data_dir,"co2_mean_2005-2023.tif"))
co2_future <- rast(file.path(data_dir,"co2_mean_2025-2045.tif"))
dist_at_past <- rast(dist_at_past_file)
dist_at_future <- rast(dist_at_future_file)
ta_gs_past <- rast(file.path(data_dir,"t2m_gs_2005_2023.tif"))-273.15
ta_gs_future <- rast(file.path(data_dir,"tas_gs_2026_2050.tif"))-273.15
tp_gs_past <- rast(file.path(data_dir,"tp_gs_2005_2023.tif"))*1000*24
tp_gs_future <- rast(file.path(data_dir,"pr_gs_2026_2050.tif"))*86400
sw_gs_past <- rast(file.path(data_dir,"ssrd_gs_2005_2023.tif"))/3600
sw_gs_future <- rast(file.path(data_dir,"rsds_gs_2026_2050.tif"))
et_gs_past <- rast(file.path(data_dir,"et_gs_2005_2023.tif"))
et_change <- rast(file.path(data_dir,"et_at_past_fut_change.tif"))
le_gs_past <- et_gs_past/3600*2.45e6
le_change <- et_change/3600*2.45e6

sen_slope_raster <- function(x) terra::app(x,function(z) {
  if(anyNA(z)) NA_real_ else as.numeric(trend::sens.slope(z)$estimate)
})

message("Calculating environmental trends...")
ta_gs_at_past <- sen_slope_raster(ta_gs_past)
ta_gs_at_future <- sen_slope_raster(ta_gs_future)
p_gs_at_past <- sen_slope_raster(tp_gs_past)
p_gs_at_future <- sen_slope_raster(tp_gs_future)
sw_gs_at_past <- sen_slope_raster(sw_gs_past)
sw_gs_at_future <- sen_slope_raster(sw_gs_future)
le_gs_at_past <- sen_slope_raster(le_gs_past)
le_gs_at_future <- le_gs_at_past+le_change

continuous_rasters <- list(
  nee_baseline=nee_baseline,mat_baseline=mat_baseline,map_baseline=map_baseline,
  gpp_at_baseline=gpp_at_baseline,co2_past=co2_past,co2_future=co2_future,
  dist_at_past=dist_at_past,dist_at_future=dist_at_future,
  ta_gs_at_past=ta_gs_at_past,ta_gs_at_future=ta_gs_at_future,
  p_gs_at_past=p_gs_at_past,p_gs_at_future=p_gs_at_future,
  sw_gs_at_past=sw_gs_at_past,sw_gs_at_future=sw_gs_at_future,
  le_gs_at_past=le_gs_at_past,le_gs_at_future=le_gs_at_future
)
continuous_rasters <- lapply(continuous_rasters,terra::resample,y=land_cover,method="bilinear")
list2env(continuous_rasters,envir=environment())
zone_group <- terra::resample(zone_group,land_cover,method="near")

terra::writeRaster(ta_gs_at_future-ta_gs_at_past,file.path(output_dir,"Change_TA_gs_at.tif"),overwrite=TRUE)
terra::writeRaster(p_gs_at_future-p_gs_at_past,file.path(output_dir,"Change_P_gs_at.tif"),overwrite=TRUE)
terra::writeRaster(le_gs_at_future-le_gs_at_past,file.path(output_dir,"Change_LE_gs_at.tif"),overwrite=TRUE)

# Training data and weights
trend_driver <- read.csv(file.path(data_dir,"potential_trend_driver.csv"))
sites_long <- read.csv(file.path(data_dir,"long_term_ec_sites.csv"))
if(!"category" %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(
    sites_long[,c("site_ID","IGBP","Climate.class","category","NEEpartition","source","AI")],
    by="site_ID"
  )
}

site_uncertainty <- trend_driver %>% group_by(site_ID) %>% summarise(
  NEE_at_sd_min=if(all(is.na(NEE_at_sd))) NA_real_ else min(NEE_at_sd,na.rm=TRUE),
  .groups="drop"
)
trend_driver <- trend_driver %>% left_join(site_uncertainty,by="site_ID") %>% mutate(
  NEE_at_sd_scaled=NEE_at_sd/(NEE_at_sd_min*sqrt(nyear)),
  model_weight=1/NEE_at_sd_scaled^2,TA_at=TA_gs_at,LE_at=LE_gs_at
) %>% filter(is.finite(model_weight),model_weight>0)

empty_map <- function() {
  x <- zone_group
  terra::values(x) <- NA_real_
  x
}

nee_past_maps <- setNames(lapply(c(method_ids,"Ensemble"),function(x) empty_map()),c(method_ids,"Ensemble"))
nee_future_maps <- setNames(lapply(c(method_ids,"Ensemble"),function(x) empty_map()),c(method_ids,"Ensemble"))
zone_codes <- setNames(c(2,1,3,4,5),climate_vegs)

raster_values_at <- function(x,cells) {
  if(terra::nlyr(x)!=1) stop("Expected one raster layer, found ",terra::nlyr(x),".")
  as.numeric(terra::extract(x,cells)[[1]])
}

make_projection_data <- function(cells,period=c("past","future")) {
  period <- match.arg(period)
  is_future <- period=="future"
  data.frame(
    cell_id=cells,
    CO2=raster_values_at(if(is_future) co2_future else co2_past,cells),
    TA_gs_at=raster_values_at(if(is_future) ta_gs_at_future else ta_gs_at_past,cells),
    P_gs_at=raster_values_at(if(is_future) p_gs_at_future else p_gs_at_past,cells),
    SW_gs_at=raster_values_at(if(is_future) sw_gs_at_future else sw_gs_at_past,cells),
    LE_gs_at=raster_values_at(if(is_future) le_gs_at_future else le_gs_at_past,cells),
    DIST_at=raster_values_at(if(is_future) dist_at_future else dist_at_past,cells),
    MAP=raster_values_at(map_baseline,cells),MAT=raster_values_at(mat_baseline,cells),
    NEE=raster_values_at(nee_baseline,cells),GPP_at_mean=raster_values_at(gpp_at_baseline,cells),
    IGBP_id=raster_values_at(land_cover,cells)
  ) %>% left_join(df_IGBP,by="IGBP_id") %>% mutate(TA_at=TA_gs_at,LE_at=LE_gs_at)
}

partial_results <- list()
model_diagnostics <- list()

for(i in seq_along(climate_vegs)) {
  climate_veg <- climate_vegs[i]
  message("Projecting ",climate_veg,"...")
  data_group <- get_group(trend_driver,climate_veg)
  cells <- which(terra::values(zone_group,mat=FALSE)==zone_codes[[climate_veg]])
  project_past <- make_projection_data(cells,"past")
  project_future <- make_projection_data(cells,"future")

  if(climate_veg %in% swc_groups) {
    swc_prepared <- prepare_training_data(data_group,"SWC_gs_at",swc_predictors)
    swc_rf <- list(
      model=fit_rf_component(swc_prepared$data,"SWC_gs_at",swc_prepared$predictors,1000L+i),
      train=swc_prepared$data,predictors=swc_prepared$predictors
    )
    project_past$SWC_gs_at <- predict_one_model(swc_rf,project_past,"RF")
    project_future$SWC_gs_at <- predict_one_model(swc_rf,project_future,"RF")
  }

  nee_predictors <- get_nee_predictors(key_vars[i])
  nee_models <- fit_models(data_group,"NEE_at",nee_predictors,128L+i)
  key_predictors <- all.vars(as.formula(paste("~",key_vars[i])))
  partial_results[[climate_veg]] <- calculate_pdp(
    nee_models,climate_veg,"NEE_at",key_predictors,n_grid=50L
  )

  for(method in method_ids) {
    pred_past <- predict_one_model(nee_models[[method]],project_past,method)
    pred_future <- predict_one_model(nee_models[[method]],project_future,method)
    terra::values(nee_past_maps[[method]])[cells] <- pred_past
    terra::values(nee_future_maps[[method]])[cells] <- pred_future
  }

  model_diagnostics[[climate_veg]] <- tibble(
    climate_veg=climate_veg,method="Linear CO2 + RF",
    converged=nee_models$CO2_RF$converged,
    iterations=nee_models$CO2_RF$iterations
  )
}

# Partial-dependence outputs for climate-specific key variables
df_partial <- bind_rows(partial_results)
write.csv(df_partial,file.path(output_dir,"NEE_partial_dependence_key_variables_LM_RF_CO2RF.csv"),row.names=FALSE)

for(group_i in climate_vegs) {
  plot_data <- df_partial %>% filter(climate_veg==group_i)
  if(!nrow(plot_data)) next
  p <- ggplot(plot_data,aes(value,partial_effect,color=method,linetype=method))+
    geom_hline(yintercept=0,color="grey70",linewidth=.3)+
    geom_line(linewidth=.8)+facet_wrap(~predictor,scales="free_x")+
    labs(x=NULL,y=expression("Partial effect on NEE trend (g C m"^{-2}*" yr"^{-2}*")"),
         color="Method",linetype="Method",title=group_i)+theme_classic()
  print(p)
  # ggsave(file.path(output_dir,paste0("NEE_partial_key_variables_LM_RF_CO2RF_",group_i,".pdf")),p,width=10,height=6)
}

# Equal ensemble of LM and RF
nee_past_maps$Ensemble <- (nee_past_maps$LM+nee_past_maps$RF)/2
nee_future_maps$Ensemble <- (nee_future_maps$LM+nee_future_maps$RF)/2

write_method_rasters <- function(past_maps,future_maps,prefix) {
  for(method in names(past_maps)) {
    terra::writeRaster(
      future_maps[[method]]-past_maps[[method]],
      file.path(output_dir,paste0(prefix,"_change_",method,".tif")),
      overwrite=TRUE
    )
  }
}

write_method_rasters(nee_past_maps,nee_future_maps,"NEE_trend")
# write.csv(bind_rows(model_diagnostics),file.path(output_dir,"CO2_RF_convergence_diagnostics.csv"),row.names=FALSE)

plot(nee_future_maps$Ensemble-nee_past_maps$Ensemble,main="Ensemble")
plot(nee_future_maps$LM-nee_past_maps$LM,main="LM")
plot(nee_future_maps$RF-nee_past_maps$RF,main="RF")
plot(nee_future_maps$CO2_RF-nee_past_maps$CO2_RF,main="Linear CO2 + RF")

df_direction_agreement <- purrr::map_dfr(climate_vegs,function(climate_veg) {
  cells <- which(terra::values(zone_group,mat=FALSE)==zone_codes[[climate_veg]])
  change <- cbind(
    LM=raster_values_at(nee_future_maps$LM-nee_past_maps$LM,cells),
    RF=raster_values_at(nee_future_maps$RF-nee_past_maps$RF,cells) #,
    # CO2_RF=raster_values_at(nee_future_maps$CO2_RF-nee_past_maps$CO2_RF,cells) #,
    # Ensemble=raster_values_at(nee_future_maps$Ensemble-nee_past_maps$Ensemble,cells)
  )
  valid <- rowSums(is.finite(change))==ncol(change)
  same_increase <- valid & rowSums(change>0,na.rm=TRUE)==ncol(change)
  same_decrease <- valid & rowSums(change<0,na.rm=TRUE)==ncol(change)
  same_direction <- same_increase | same_decrease
  n_valid <- sum(valid)
  tibble(
    climate_veg=climate_veg,
    n_cells=length(cells),
    n_valid=n_valid,
    percent_same_direction=ifelse(n_valid>0,100*sum(same_direction)/n_valid,NA_real_),
    percent_same_increase=ifelse(n_valid>0,100*sum(same_increase)/n_valid,NA_real_),
    percent_same_decrease=ifelse(n_valid>0,100*sum(same_decrease)/n_valid,NA_real_)
  )
})

print(df_direction_agreement)
# It is really CO2 effect determining nearly all uncertainty in future change projections. 
# LM model is consistent in projection change direction with LM_CO2_RF, and ensemble LM and RF
# LM_CO2_RF are very different from RF in projection direction in wet-forest and wet-warm ecosystems. 

# write.csv(df_direction_agreement,file.path(output_dir,"NEE_change_direction_agreement_LM_RF_CO2RF_ensemble.csv"),row.names=FALSE)

message("Finished. Outputs are in: ",output_dir)

files <- c(
  Ensemble=file.path("data2026","NEE_trend_change_Ensemble.tif"),
  CO2_RF=file.path("data2026","NEE_trend_change_CO2_RF.tif"),
  RF=file.path("data2026","NEE_trend_change_RF.tif"),
  LM=file.path("data2026","NEE_trend_change_LM.tif")
)

r <- terra::rast(unname(files))
names(r) <- names(files)

v <- terra::values(r,mat=TRUE)
v <- v[rowSums(is.finite(v))==ncol(v),,drop=FALSE]

round(cor(v),3)



