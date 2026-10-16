# Global NEE-trend projection using weighted LM, weighted RF, their equal ensemble,
# and a linear-CO2 plus RF-residual hybrid.
# Junna Wang
#
# Final outputs:
#   1. NEE_partial_dependence_RF.csv
#   2. SWC_gs_at_by_climate_vegetation.csv
#   3. NEE_trend_change_LM.tif
#   4. NEE_trend_change_RF.tif
#   5. NEE_trend_change_Ensemble.tif
#   6. NEE_trend_change_CO2_RF.tif
#   7. NEE_change_bootstrap_summary_LM_RF_Ensemble_CO2_RF.csv

rm(list=ls())
library(librarian)
shelf(tidyverse,terra,trend,ranger)
source(file.path("scripts2026", "00_functions.R"))
# ==============================================================================
# 1. Settings
# ==============================================================================

data_dir <- "data2026"
output_dir <- "data2026"
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

dist_at_past_file <- file.path(data_dir,"LAI_deficit_trends_DBF.tif")
dist_at_future_file <- dist_at_past_file

n_boot <- 1000L
bootstrap_seed <- 20260928L
checkpoint_every <- 25L
resume_from_checkpoint <- TRUE
checkpoint_schema <- 5L

climate_vegs <- c(
  "cold_evergreen_forest",
  "cold_deciduous_forest",
  "cold_non_forest",
  "dry",
  "wet-warm"
)
zone_codes <- setNames(c(2,1,3,4,5),climate_vegs)

# Climate-specific time-varying predictors used for NEE projection.
key_vars <- c(
  "CO2 + TA_gs_at",
  "CO2",
  "SWC_gs_at",
  "CO2 + TA_gs_at + SWC_gs_at",
  "CO2 + LE_at")

nee_common_predictors <- c("MAP","MAT","IGBP","GPP_at_mean","NEE")
swc_predictors <- c("CO2","TA_gs_at","SW_gs_at","P_gs_at","MAP","MAT","NEE","IGBP")

# SWC is only required by the NEE models for these two groups.
swc_groups <- c("cold_non_forest","dry")

m_try <- 3L
min_node_size <- 2L
num_trees <- 500L
sample_fraction <- 0.8
replace_samples <- TRUE
split_rule <- "extratrees"
hybrid_max_iter <- 30L
hybrid_tolerance <- 5e-3

# ==============================================================================
# 2. Model and data helper functions
# ==============================================================================

get_nee_predictors <- function(key_var) {
  unique(c(all.vars(as.formula(paste("~",key_var))),nee_common_predictors))
}

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
  out <- out[is.finite(out$model_weight)&out$model_weight>0,,drop=FALSE]
  out$IGBP <- factor(out$IGBP)
  out <- droplevels(out)
  predictors <- predictors[vapply(predictors,function(v) n_distinct(out[[v]])>1,logical(1))]
  # if(!"CO2" %in% predictors) stop("CO2 is not estimable for ",response,".")
  if(nrow(out)<10) stop("Too few complete training rows for ",response,".")
  out$.model_weight <- out$model_weight/mean(out$model_weight)
  list(data=out,predictors=predictors)
}

align_factor_levels <- function(train,newdata,predictors) {
  out <- newdata
  factor_vars <- predictors[vapply(train[predictors],is.factor,logical(1))]
  for(v in factor_vars) {
    levels_train <- levels(train[[v]])
    fallback <- names(which.max(table(train[[v]])))
    x <- as.character(out[[v]])
    x[!is.na(x)&!x %in% levels_train] <- fallback
    out[[v]] <- factor(x,levels=levels_train)
  }
  out
}

fit_rf_component <- function(data,response,predictors,seed) {
  ranger::ranger(
    reformulate(predictors,response=response),
    data=data,
    num.trees=num_trees,
    mtry=min(m_try,length(predictors)),
    min.node.size=min_node_size,
    sample.fraction=sample_fraction,
    replace=replace_samples,
    splitrule=split_rule,
    case.weights=data$.model_weight,
    importance="permutation",
    oob.error=FALSE,
    seed=seed,
    num.threads=1
  )
}

fit_rf_model <- function(data,response,predictors,seed) {
  prepared <- prepare_training_data(data,response,predictors)
  list(
    model=fit_rf_component(prepared$data,response,prepared$predictors,seed),
    train=prepared$data,
    predictors=prepared$predictors
  )
}

fit_lm_rf <- function(data,response,predictors,seed) {
  prepared <- prepare_training_data(data,response,predictors)
  train <- prepared$data
  predictors <- prepared$predictors
  list(
    LM=list(
      model=lm(reformulate(predictors,response=response),data=train,weights=.model_weight),
      train=train,
      predictors=predictors
    ),
    RF=list(
      model=fit_rf_component(train,response,predictors,seed),
      train=train,
      predictors=predictors
    )
  )
}

fit_co2_component <- function(y,data) {
  d <- data
  d$.hybrid_response <- y
  lm(.hybrid_response~CO2,data=d,weights=.model_weight)
}

# CO2 enters only when selected as a climate-group key variable. The RF fits
# the remaining factors; weighted OOB predictions backfit the CO2 component.
fit_hybrid_model <- function(
    data,response,predictors,seed,co2_selected=("CO2" %in% predictors),
    max_iter=hybrid_max_iter,tol=hybrid_tolerance) {
  hybrid_predictors <- if(co2_selected) unique(c("CO2",setdiff(predictors,"CO2")))
                       else setdiff(predictors,"CO2")
  prepared <- prepare_training_data(data,response,hybrid_predictors)
  train <- prepared$data
  if(co2_selected&&!"CO2" %in% prepared$predictors) {
    stop("CO2 is not estimable for hybrid response ",response,".")
  }
  rf_predictors <- setdiff(prepared$predictors,"CO2")
  if(length(rf_predictors)==0L) stop("Hybrid model needs a non-CO2 predictor.")

  y <- train[[response]]
  w <- train$.model_weight
  co2_model <- NULL
  co2_train <- numeric(nrow(train))
  converged <- TRUE

  if(co2_selected) {
    co2_model <- fit_co2_component(y,train)
    co2_train <- as.numeric(predict(co2_model,newdata=train))
    y_mean <- weighted.mean(y,w)
    y_scale <- max(sqrt(weighted.mean((y-y_mean)^2,w)),.Machine$double.eps)
    converged <- FALSE

    for(iteration in seq_len(max_iter)) {
      rf_data <- train
      rf_data$.hybrid_residual <- y-co2_train
      rf_model <- ranger::ranger(
        reformulate(rf_predictors,response=".hybrid_residual"),
        data=rf_data,num.trees=num_trees,
        mtry=min(m_try,length(rf_predictors)),min.node.size=min_node_size,
        sample.fraction=sample_fraction,replace=replace_samples,
        splitrule=split_rule,case.weights=rf_data$.model_weight,
        oob.error=TRUE,seed=seed+iteration,num.threads=1
      )
      rf_oob_raw <- rf_model$predictions
      if(anyNA(rf_oob_raw)) stop("Missing OOB predictions in hybrid model.")
      rf_oob <- rf_oob_raw-weighted.mean(rf_oob_raw,w)
      new_co2_model <- fit_co2_component(y-rf_oob,train)
      new_co2_train <- as.numeric(predict(new_co2_model,newdata=train))
      change <- sqrt(weighted.mean((new_co2_train-co2_train)^2,w))/y_scale
      co2_model <- new_co2_model
      co2_train <- new_co2_train
      if(change<tol) {
        converged <- TRUE
        break
      }
    }
    if(!converged) warning("Hybrid LM-RF did not converge after ",max_iter," iterations.")
  }

  rf_data <- train
  rf_data$.hybrid_residual <- y-co2_train
  rf_model <- ranger::ranger(
    reformulate(rf_predictors,response=".hybrid_residual"),
    data=rf_data,num.trees=num_trees,
    mtry=min(m_try,length(rf_predictors)),min.node.size=min_node_size,
    sample.fraction=sample_fraction,replace=replace_samples,
    splitrule=split_rule,case.weights=rf_data$.model_weight,
    oob.error=TRUE,seed=seed+max_iter+1L,num.threads=1
  )

  list(
    has_co2=co2_selected,co2_model=co2_model,rf_model=rf_model,
    rf_center=weighted.mean(rf_model$predictions,w),train=train,
    rf_predictors=rf_predictors,converged=converged
  )
}

predict_hybrid_model <- function(object,newdata) {
  out <- rep(NA_real_,nrow(newdata))
  predictors <- c(if(object$has_co2) "CO2",object$rf_predictors)
  d <- align_factor_levels(object$train,newdata,predictors)
  ok <- complete.cases(d[,predictors,drop=FALSE])
  if(!any(ok)) return(out)
  co2_prediction <- if(object$has_co2)
    as.numeric(predict(object$co2_model,newdata=d[ok,,drop=FALSE]))
  else numeric(sum(ok))
  out[ok] <- co2_prediction+
    as.numeric(predict(object$rf_model,data=d[ok,,drop=FALSE])$predictions)-
    object$rf_center
  out
}

predict_one_model <- function(object,newdata,method=c("LM","RF")) {
  method <- match.arg(method)
  out <- rep(NA_real_,nrow(newdata))
  d <- align_factor_levels(object$train,newdata,object$predictors)
  ok <- complete.cases(d[,object$predictors,drop=FALSE])
  if(!any(ok)) return(out)
  if(method=="LM") {
    out[ok] <- as.numeric(predict(object$model,newdata=d[ok,,drop=FALSE]))
  } else {
    out[ok] <- as.numeric(predict(object$model,data=d[ok,,drop=FALSE])$predictions)
  }
  out
}

# Assign each site to a fixed IGBP x observation-count stratum before bootstrap.
# The data were not thinned, so n_observations is the number of rows per site.
make_bootstrap_registry <- function(d,climate_veg) {
  site_key <- d %>%
    filter(!is.na(site_ID),!is.na(IGBP)) %>%
    count(site_ID,IGBP,name="n_observations")

  if(nrow(site_key)==0L) stop("No sites available for bootstrap in ",climate_veg,".")
  if(any(duplicated(site_key$site_ID))) {
    stop("At least one site has multiple IGBP values in ",climate_veg,".")
  }

  site_key %>%
    mutate(
      observation_class=case_when(
        n_observations<=2L ~ "n<=2",
        n_observations<=6L ~ "3<=n<=6",
        TRUE ~ "n>6"
      ),
      observation_class=factor(
        observation_class,
        levels=c("n<=2","3<=n<=6","n>6")
      )
    )
}

# Site-cluster bootstrap, stratified by IGBP and observation-count class.
# A sampled site contributes all of its unthinned rows.
bootstrap_sites <- function(d,registry,seed) {
  set.seed(seed)
  draws <- registry %>%
    group_by(IGBP,observation_class) %>%
    group_modify(~tibble(
      source_site_ID=sample(.x$site_ID,size=nrow(.x),replace=TRUE),
      draw_within_stratum=seq_len(nrow(.x))
    )) %>%
    ungroup() %>%
    mutate(bootstrap_draw=row_number())

  bind_rows(lapply(seq_len(nrow(draws)),function(j) {
    source_id <- draws$source_site_ID[j]
    z <- d[d$site_ID==source_id,,drop=FALSE]
    z$.bootstrap_site <- paste(
      draws$IGBP[j],
      as.character(draws$observation_class[j]),
      source_id,
      draws$bootstrap_draw[j],
      sep="__"
    )
    z
  }))
}

# Weighted RF partial dependence over the observed 2nd-98th percentile range.
calculate_rf_pdp <- function(object,climate_veg,plot_predictors,n_grid=50L) {
  train <- object$train
  weight <- as.numeric(train$.model_weight)
  predictors <- intersect(plot_predictors,object$predictors)
  baseline <- weighted.mean(predict_one_model(object,train,"RF"),weight,na.rm=TRUE)

  map_dfr(predictors,function(variable) {
    x <- as.numeric(train[[variable]])
    limits <- as.numeric(quantile(x,c(.02,.98),na.rm=TRUE,names=FALSE))
    if(!all(is.finite(limits))||limits[1]==limits[2]) return(tibble())

    map_dfr(seq(limits[1],limits[2],length.out=n_grid),function(value) {
      newdata <- train
      newdata[[variable]] <- value
      pdp_prediction <- weighted.mean(
        predict_one_model(object,newdata,"RF"),
        weight,
        na.rm=TRUE
      )
      tibble(
        climate_veg=climate_veg,
        method="RF",
        predictor=variable,
        value=value,
        prediction=pdp_prediction,
        partial_effect=pdp_prediction-baseline
      )
    })
  })
}

sen_slope_raster <- function(x) {
  terra::app(x,function(z) {
    if(anyNA(z)) NA_real_ else as.numeric(trend::sens.slope(z)$estimate)
  })
}

raster_values_at <- function(x,cells) {
  if(terra::nlyr(x)!=1) stop("Expected one raster layer, found ",terra::nlyr(x),".")
  as.numeric(terra::values(x,mat=FALSE)[cells])
}

empty_map <- function(template) {
  x <- template
  terra::values(x) <- NA_real_
  x
}

area_weighted_mean <- function(values,area) {
  ok <- is.finite(values)&is.finite(area)&area>0
  if(!any(ok)) return(NA_real_)
  weighted.mean(values[ok],area[ok])
}

# ==============================================================================
# 3. Read rasters and calculate environmental trends once
# ==============================================================================

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

message("Calculating environmental trends...")
ta_gs_at_past <- sen_slope_raster(ta_gs_past)
ta_gs_at_future <- sen_slope_raster(ta_gs_future)
p_gs_at_past <- sen_slope_raster(tp_gs_past)
p_gs_at_future <- sen_slope_raster(tp_gs_future)
sw_gs_at_past <- sen_slope_raster(sw_gs_past)
sw_gs_at_future <- sen_slope_raster(sw_gs_future)
le_gs_at_past <- sen_slope_raster(et_gs_past/3600*2.45e6)
le_gs_at_future <- le_gs_at_past+et_change/3600*2.45e6

# Put all continuous predictors on the land-cover grid. Categorical rasters use nearest neighbor.
continuous_rasters <- list(
  nee_baseline=nee_baseline,
  mat_baseline=mat_baseline,
  map_baseline=map_baseline,
  gpp_at_baseline=gpp_at_baseline,
  co2_past=co2_past,
  co2_future=co2_future,
  dist_at_past=dist_at_past,
  dist_at_future=dist_at_future,
  ta_gs_at_past=ta_gs_at_past,
  ta_gs_at_future=ta_gs_at_future,
  p_gs_at_past=p_gs_at_past,
  p_gs_at_future=p_gs_at_future,
  sw_gs_at_past=sw_gs_at_past,
  sw_gs_at_future=sw_gs_at_future,
  le_gs_at_past=le_gs_at_past,
  le_gs_at_future=le_gs_at_future
)
continuous_rasters <- lapply(continuous_rasters,terra::resample,y=land_cover,method="bilinear")
list2env(continuous_rasters,envir=environment())
zone_group <- terra::resample(zone_group,land_cover,method="near")
cell_area_values <- terra::values(terra::cellSize(zone_group,unit="km"),mat=FALSE)

# ==============================================================================
# 4. Read training data and calculate fixed model weights
# ==============================================================================

trend_driver <- read.csv(file.path(data_dir,"potential_trend_driver.csv"))
sites_long <- read.csv(file.path(data_dir,"long_term_ec_sites.csv"))
meta_cols <- c("IGBP","Climate.class","category","NEEpartition","source","AI")
missing_meta <- setdiff(meta_cols,names(trend_driver))
if(length(missing_meta)>0) {
  trend_driver <- trend_driver %>%
    left_join(sites_long[,c("site_ID",missing_meta),drop=FALSE],by="site_ID")
}

site_uncertainty <- trend_driver %>%
  group_by(site_ID) %>%
  summarise(
    NEE_at_sd_min=if(all(is.na(NEE_at_sd))) NA_real_ else min(NEE_at_sd,na.rm=TRUE),
    .groups="drop"
  )

trend_driver <- trend_driver %>%
  left_join(site_uncertainty,by="site_ID") %>%
  mutate(
    NEE_at_sd_scaled=NEE_at_sd/(NEE_at_sd_min*sqrt(nyear)),
    model_weight=1/NEE_at_sd_scaled^2,
    TA_at=TA_gs_at,
    LE_at=LE_gs_at
  ) %>%
  filter(is.finite(model_weight),model_weight>0)

# ==============================================================================
# 5. Construct reusable past and future projection data
# ==============================================================================

make_projection_data <- function(cells,period=c("past","future")) {
  future <- match.arg(period)=="future"
  data.frame(
    cell_id=cells,
    CO2=raster_values_at(if(future) co2_future else co2_past,cells),
    TA_gs_at=raster_values_at(if(future) ta_gs_at_future else ta_gs_at_past,cells),
    P_gs_at=raster_values_at(if(future) p_gs_at_future else p_gs_at_past,cells),
    SW_gs_at=raster_values_at(if(future) sw_gs_at_future else sw_gs_at_past,cells),
    LE_gs_at=raster_values_at(if(future) le_gs_at_future else le_gs_at_past,cells),
    DIST_at=raster_values_at(if(future) dist_at_future else dist_at_past,cells),
    MAP=raster_values_at(map_baseline,cells),
    MAT=raster_values_at(mat_baseline,cells),
    NEE=raster_values_at(nee_baseline,cells),
    GPP_at_mean=raster_values_at(gpp_at_baseline,cells),
    IGBP_id=raster_values_at(land_cover,cells)
  ) %>%
    left_join(df_IGBP,by="IGBP_id") %>%
    mutate(TA_at=TA_gs_at,LE_at=LE_gs_at)
}

group_data <- setNames(
  lapply(climate_vegs,function(g) get_group(trend_driver,g)),
  climate_vegs
)

bootstrap_registry <- setNames(
  lapply(climate_vegs,function(g) make_bootstrap_registry(group_data[[g]],g)),
  climate_vegs
)

projection_frames <- setNames(
  lapply(climate_vegs,function(g) {
    cells <- which(terra::values(zone_group,mat=FALSE)==zone_codes[[g]])
    list(
      cells=cells,
      past=make_projection_data(cells,"past"),
      future=make_projection_data(cells,"future"),
      area=cell_area_values[cells]
    )
  }),
  climate_vegs
)

# ==============================================================================
# 6. Full-data RF partial dependence and LM/RF/ensemble/hybrid projections
# ==============================================================================

nee_change_maps <- list(
  LM=empty_map(zone_group),
  RF=empty_map(zone_group),
  Ensemble=empty_map(zone_group),
  CO2_RF=empty_map(zone_group)
)
partial_results <- vector("list",length(climate_vegs))
names(partial_results) <- climate_vegs
swc_results_past <- vector("list",length(climate_vegs))
swc_results_future <- vector("list",length(climate_vegs))
names(swc_results_past) <- climate_vegs
names(swc_results_future) <- climate_vegs

for(i in seq_along(climate_vegs)) {
  climate_veg <- climate_vegs[i]
  message("Full-data projection: ",climate_veg)

  data_group <- group_data[[climate_veg]]
  frame <- projection_frames[[climate_veg]]
  project_past <- frame$past
  project_future <- frame$future

  # Calculate SWC for every climate-vegetation group. The NEE models only use
  # these predictions where SWC_gs_at is present in the group's key variables.
  swc_rf <- fit_rf_model(
    data_group,
    response="SWC_gs_at",
    predictors=swc_predictors,
    seed=1000L+i
  )
  project_past$SWC_gs_at <- predict_one_model(swc_rf,project_past,"RF")
  project_future$SWC_gs_at <- predict_one_model(swc_rf,project_future,"RF")

  swc_results_past[[climate_veg]] <- tibble(
    `climate-veg`=climate_veg,
    `cell-id`=project_past$cell_id,
    SWC_gs_at=project_past$SWC_gs_at
  )
  swc_results_future[[climate_veg]] <- tibble(
    `climate-veg`=climate_veg,
    `cell-id`=project_future$cell_id,
    SWC_gs_at=project_future$SWC_gs_at
  )

  nee_models <- fit_lm_rf(
    data_group,
    response="NEE_at",
    predictors=get_nee_predictors(key_vars[i]),
    seed=2000L+i
  )
  hybrid_model <- fit_hybrid_model(
    data_group,
    response="NEE_at",
    predictors=get_nee_predictors(key_vars[i]),
    seed=3000L+i,
    co2_selected="CO2" %in% all.vars(as.formula(paste("~",key_vars[i])))
  )

  # RF partial dependence is calculated only for the climate-specific key variables.
  partial_results[[climate_veg]] <- calculate_rf_pdp(
    nee_models$RF,
    climate_veg=climate_veg,
    plot_predictors=all.vars(as.formula(paste("~",key_vars[i]))),
    n_grid=50L
  )

  change_lm <- predict_one_model(nee_models$LM,project_future,"LM")-
    predict_one_model(nee_models$LM,project_past,"LM")
  change_rf <- predict_one_model(nee_models$RF,project_future,"RF")-
    predict_one_model(nee_models$RF,project_past,"RF")
  change_ensemble <- (change_lm+change_rf)/2
  change_hybrid <- predict_hybrid_model(hybrid_model,project_future)-
    predict_hybrid_model(hybrid_model,project_past)

  terra::values(nee_change_maps$LM)[frame$cells] <- change_lm
  terra::values(nee_change_maps$RF)[frame$cells] <- change_rf
  terra::values(nee_change_maps$Ensemble)[frame$cells] <- change_ensemble
  terra::values(nee_change_maps$CO2_RF)[frame$cells] <- change_hybrid
}

# Store past, future, and future-minus-past SWC projections in one row per cell.
df_swc_gs_at <- bind_rows(swc_results_past) %>%
  rename(SWC_gs_at_past=SWC_gs_at) %>%
  full_join(
    bind_rows(swc_results_future) %>%
      rename(SWC_gs_at_fut=SWC_gs_at),
    by=c("climate-veg","cell-id")
  ) %>%
  mutate(SWC_gs_at=SWC_gs_at_fut-SWC_gs_at_past) %>%
  select(
    `climate-veg`,
    `cell-id`,
    SWC_gs_at_fut,
    SWC_gs_at_past,
    SWC_gs_at
  )

write.csv(
  df_swc_gs_at,
  file.path(output_dir,"SWC_gs_at_change_by_climate_vegetation.csv"),
  row.names=FALSE
)

# Save the RF partial-dependence values used by the partial plots.
df_partial_RF <- bind_rows(partial_results)
write.csv(
  df_partial_RF,
  file.path(output_dir,"NEE_partial_dependence_RF.csv"),
  row.names=FALSE
)

# Display one RF partial-dependence plot for each climate-vegetation group.
partial_plots_RF <- setNames(lapply(climate_vegs,function(g) {
  p <- df_partial_RF %>%
    filter(climate_veg==g) %>%
    ggplot(aes(value,partial_effect))+
    geom_hline(yintercept=0,color="grey70",linewidth=.3)+
    geom_line(color="#00A087",linewidth=.8)+
    facet_wrap(~predictor,scales="free_x")+
    labs(
      x=NULL,
      y=expression("Partial effect on NEE trend (g C m"^{-2}*" yr"^{-2}*")"),
      title=g
    )+
    theme_classic()
  print(p)
  p
}),climate_vegs)

# Save all future-minus-past NEE-trend projection maps.
for(method in names(nee_change_maps)) {
  names(nee_change_maps[[method]]) <- paste0("NEE_trend_change_",method)
  terra::writeRaster(
    nee_change_maps[[method]],
    file.path(output_dir,paste0("NEE_trend_change_",method,".tif")),
    overwrite=TRUE
  )
}

# ==============================================================================
# 7. Stratified site-cluster bootstrap of LM, RF, ensemble, and CO2-RF
# ==============================================================================

# Refit all models in every bootstrap and retain area-weighted projections.
bootstrap_methods <- c("LM","RF","Ensemble","CO2_RF")
global_boot <- matrix(
  NA_real_,
  nrow=n_boot,
  ncol=length(bootstrap_methods),
  dimnames=list(bootstrap=seq_len(n_boot),method=bootstrap_methods)
)
climate_boot <- array(
  NA_real_,
  dim=c(n_boot,length(climate_vegs),length(bootstrap_methods)),
  dimnames=list(
    bootstrap=seq_len(n_boot),
    climate_veg=climate_vegs,
    method=bootstrap_methods
  )
)

checkpoint_file <- file.path(output_dir,"LM_RF_ensemble_bootstrap_checkpoint.rds")
start_boot <- 1L

if(resume_from_checkpoint&&file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)
  compatible <- identical(checkpoint$schema,checkpoint_schema)&&
    identical(checkpoint$n_boot,n_boot)&&
    identical(checkpoint$climate_vegs,climate_vegs)&&
    identical(checkpoint$bootstrap_methods,bootstrap_methods)

  if(compatible) {
    global_boot <- checkpoint$global_boot
    climate_boot <- checkpoint$climate_boot
    start_boot <- checkpoint$last_boot+1L
    message("Resuming at bootstrap ",start_boot,".")
  } else {
    warning("Existing checkpoint is incompatible; restarting at bootstrap 1.")
  }
}

if(start_boot<=n_boot) {
  for(b in seq.int(start_boot,n_boot)) {
    message("LM/RF/ensemble/CO2-RF bootstrap ",b," / ",n_boot)
    bootstrap_ok <- TRUE
    global_numerator <- setNames(numeric(length(bootstrap_methods)),bootstrap_methods)
    global_denominator <- setNames(numeric(length(bootstrap_methods)),bootstrap_methods)

    for(i in seq_along(climate_vegs)) {
      climate_veg <- climate_vegs[i]
      frame <- projection_frames[[climate_veg]]

      result <- tryCatch({
        boot_data <- bootstrap_sites(
          group_data[[climate_veg]],
          registry=bootstrap_registry[[climate_veg]],
          seed=bootstrap_seed+b*1000L+i
        )
        project_past <- frame$past
        project_future <- frame$future

        if(climate_veg %in% swc_groups) {
          swc_rf <- fit_rf_model(
            boot_data,
            response="SWC_gs_at",
            predictors=swc_predictors,
            seed=bootstrap_seed+b*10000L+i
          )
          project_past$SWC_gs_at <- predict_one_model(swc_rf,project_past,"RF")
          project_future$SWC_gs_at <- predict_one_model(swc_rf,project_future,"RF")
        }

        nee_models <- fit_lm_rf(
          boot_data,
          response="NEE_at",
          predictors=get_nee_predictors(key_vars[i]),
          seed=bootstrap_seed+b*10000L+100L+i
        )
        hybrid_model <- fit_hybrid_model(
          boot_data,
          response="NEE_at",
          predictors=get_nee_predictors(key_vars[i]),
          seed=bootstrap_seed+b*10000L+200L+i,
          co2_selected="CO2" %in% all.vars(as.formula(paste("~",key_vars[i])))
        )
        change_lm <- predict_one_model(nee_models$LM,project_future,"LM")-
          predict_one_model(nee_models$LM,project_past,"LM")
        change_rf <- predict_one_model(nee_models$RF,project_future,"RF")-
          predict_one_model(nee_models$RF,project_past,"RF")
        change_hybrid <- predict_hybrid_model(hybrid_model,project_future)-
          predict_hybrid_model(hybrid_model,project_past)
        list(
          LM=change_lm,
          RF=change_rf,
          Ensemble=(change_lm+change_rf)/2,
          CO2_RF=change_hybrid
        )
      },error=function(e) {
        warning(
          "Bootstrap ",b,", group ",climate_veg," failed: ",
          conditionMessage(e)
        )
        NULL
      })

      if(is.null(result)) {
        bootstrap_ok <- FALSE
        break
      }

      for(method in bootstrap_methods) {
        values <- result[[method]]
        valid <- is.finite(values)&is.finite(frame$area)&frame$area>0
        if(!any(valid)) {
          bootstrap_ok <- FALSE
          warning(
            "Bootstrap ",b,", group ",climate_veg,", method ",method,
            " has no valid projected cells."
          )
          break
        }

        climate_boot[b,climate_veg,method] <- weighted.mean(
          values[valid],
          frame$area[valid]
        )
        global_numerator[method] <- global_numerator[method]+
          sum(values[valid]*frame$area[valid])
        global_denominator[method] <- global_denominator[method]+
          sum(frame$area[valid])
      }
      if(!bootstrap_ok) break
    }

    if(bootstrap_ok&&all(global_denominator>0)) {
      global_boot[b,] <- global_numerator/global_denominator
    } else {
      climate_boot[b,,] <- NA_real_
      global_boot[b,] <- NA_real_
    }

    if(b%%checkpoint_every==0L||b==n_boot) {
      saveRDS(
        list(
          schema=checkpoint_schema,
          n_boot=n_boot,
          climate_vegs=climate_vegs,
          bootstrap_methods=bootstrap_methods,
          last_boot=b,
          global_boot=global_boot,
          climate_boot=climate_boot
        ),
        checkpoint_file,
        compress=TRUE
      )
    }
  }
}

# ==============================================================================
# 8. Summarise area-averaged bootstrap results for all methods
# ==============================================================================

summarise_bootstrap <- function(values,climate_veg,method) {
  values <- values[is.finite(values)]
  tibble(
    method=method,
    climate_veg=climate_veg,
    n_boot=length(values),
    mean=if(length(values)>0) mean(values) else NA_real_,
    q025=if(length(values)>0) quantile(values,.025,names=FALSE) else NA_real_,
    q25=if(length(values)>0) quantile(values,.25,names=FALSE) else NA_real_,
    q50=if(length(values)>0) quantile(values,.50,names=FALSE) else NA_real_,
    q75=if(length(values)>0) quantile(values,.75,names=FALSE) else NA_real_,
    q975=if(length(values)>0) quantile(values,.975,names=FALSE) else NA_real_
  )
}

df_all_area_mean_NEE_change <- map_dfr(bootstrap_methods,function(method) {
  bind_rows(
    summarise_bootstrap(global_boot[,method],"Global",method),
    map_dfr(climate_vegs,function(g) summarise_bootstrap(
      climate_boot[,g,method],g,method
    ))
  )
}) %>%
  mutate(
    method=factor(method,levels=bootstrap_methods),
    climate_veg=factor(climate_veg,levels=c("Global",climate_vegs))
  ) %>%
  arrange(method,climate_veg) %>%
  mutate(method=as.character(method),climate_veg=as.character(climate_veg))

print(df_all_area_mean_NEE_change)
write.csv(
  df_all_area_mean_NEE_change,
  file.path(output_dir,"NEE_change_bootstrap_summary_LM_RF_Ensemble_CO2_RF.csv"),
  row.names=FALSE
)

if(file.exists(checkpoint_file)) file.remove(checkpoint_file)
message("Finished. Outputs are in: ",output_dir)
