# Export per-group future projection frames for applicability-domain checks.
# This reconstructs the df_project_fut objects used in 03_04 with plain vector
# raster extraction, then writes data2026/proj_fut_<group>.csv.

library(tidyverse)
library(terra)
library(ranger)

source(file.path("scripts2026", "review_response", "00_review_helpers.R"))
source(file.path("scripts2026", "00_functions.R"))

sens_slope_est <- function(x) {
  ok <- is.finite(x)
  y <- x[ok]
  time <- seq_along(x)[ok]
  if (length(y) < 2) return(NA_real_)
  ij <- utils::combn(seq_along(y), 2)
  dt <- time[ij[2, ]] - time[ij[1, ]]
  dy <- y[ij[2, ]] - y[ij[1, ]]
  stats::median(dy[dt != 0] / dt[dt != 0], na.rm = TRUE)
}

cell_values <- function(r, cells) {
  v <- terra::values(r)
  if (is.matrix(v) || is.data.frame(v)) {
    return(as.vector(v[cells, 1]))
  }
  as.vector(v[cells])
}

trend_raster <- function(path, transform = identity) {
  r <- transform(terra::rast(path))
  terra::app(r, function(x) ifelse(any(is.na(x)), NA, sens_slope_est(x)))
}

flatten_frame <- function(df) {
  df[] <- lapply(df, function(x) if (!is.null(dim(x))) as.vector(x) else x)
  as.data.frame(df)
}

message("Reading rasters and calculating future/past trend rasters...")

nee_baseline <- terra::rast(file.path(data_dir, "baseline_nee_2005-2021.tif"))
gpp_baseline <- terra::rast(file.path(data_dir, "baseline_gpp_2005-2021.tif"))
et_baseline <- terra::rast(file.path(data_dir, "baseline_et_2005-2021.tif"))
le_baseline <- et_baseline / 3600 * 2.45 * 1000000
mat_baseline <- terra::rast(file.path(data_dir, "baseline_mat.tif"))
map_baseline <- terra::rast(file.path(data_dir, "baseline_map.tif"))
gpp_at_baseline <- terra::rast(file.path(data_dir, "CEDAR_GPP_Sen_slope_ML.tif"))

land_cover <- terra::rast(file.path(data_dir, "968_Land_Cover_Class_0.25degree_layer2.tif"))
koppen_climate <- terra::rast(file.path(data_dir, "koppen_geiger_0p00833333.tif"))
aridity <- terra::rast(file.path(data_dir, "Global-AI.tif"))
zone_group <- terra::rast(file.path(data_dir, "zones_climate_vegetation_group.tif"))

koppen_climate <- terra::resample(koppen_climate, land_cover, method = "near")
aridity <- terra::resample(aridity, land_cover) / 10000

co2_past <- terra::rast(file.path(data_dir, "co2_mean_2005-2023.tif"))
co2_fut <- terra::rast(file.path(data_dir, "co2_mean_2025-2045.tif"))

ta_gs_at_past <- trend_raster(file.path(data_dir, "t2m_gs_2005_2023.tif"), function(r) r - 273.15)
ta_ngs_at_past <- trend_raster(file.path(data_dir, "t2m_ngs_2005_2023.tif"), function(r) r - 273.15)
ta_gs_at_fut <- trend_raster(file.path(data_dir, "tas_gs_2026_2050.tif"), function(r) r - 273.15)
ta_ngs_at_fut <- trend_raster(file.path(data_dir, "tas_ngs_2026_2050.tif"), function(r) r - 273.15)

tp_gs_at_past <- trend_raster(file.path(data_dir, "tp_gs_2005_2023.tif"), function(r) r * 1000 * 24)
tp_ngs_at_past <- trend_raster(file.path(data_dir, "tp_ngs_2005_2023.tif"), function(r) r * 1000 * 24)
tp_gs_at_fut <- trend_raster(file.path(data_dir, "pr_gs_2026_2050.tif"), function(r) r * 86400)
tp_ngs_at_fut <- trend_raster(file.path(data_dir, "pr_ngs_2026_2050.tif"), function(r) r * 86400)

sw_gs_at_past <- trend_raster(file.path(data_dir, "ssrd_gs_2005_2023.tif"), function(r) r / 3600)
sw_gs_at_fut <- trend_raster(file.path(data_dir, "rsds_gs_2026_2050.tif"))

et_gs_past <- terra::rast(file.path(data_dir, "et_gs_2005_2023.tif"))
et_at_past_fut_change <- terra::rast(file.path(data_dir, "et_at_past_fut_change.tif"))
le_gs_past <- et_gs_past / 3600 * 2.45 * 1000000
le_gs_at_past <- trend_raster(file.path(data_dir, "et_gs_2005_2023.tif"),
                              function(r) r / 3600 * 2.45 * 1000000)
le_gs_at_fut <- le_gs_at_past + et_at_past_fut_change / 3600 * 2.45 * 1000000

DIST_at_wet_past <- terra::rast(file.path(data_dir, "DIST_at_wet_past.tif"))
DIST_at_wet_fut <- DIST_at_wet_past

trend_layers <- list(
  ta_gs_at_past = ta_gs_at_past, ta_gs_at_fut = ta_gs_at_fut,
  ta_ngs_at_past = ta_ngs_at_past, ta_ngs_at_fut = ta_ngs_at_fut,
  tp_gs_at_past = tp_gs_at_past, tp_gs_at_fut = tp_gs_at_fut,
  tp_ngs_at_past = tp_ngs_at_past, tp_ngs_at_fut = tp_ngs_at_fut,
  sw_gs_at_past = sw_gs_at_past, sw_gs_at_fut = sw_gs_at_fut,
  le_gs_at_past = le_gs_at_past, le_gs_at_fut = le_gs_at_fut,
  DIST_at_wet_past = DIST_at_wet_past, DIST_at_wet_fut = DIST_at_wet_fut
)
trend_layers <- lapply(trend_layers, terra::resample, y = land_cover)
list2env(trend_layers, envir = environment())

trend_driver <- read.csv(file.path(data_dir, "potential_trend_driver.csv"))
sites_long <- read.csv(file.path(data_dir, "long_term_ec_sites.csv"))
sites_long <- sites_long %>% mutate(source = ifelse(source == "AmeriFlux_BASE", "Eddyproc", "Oneflux"))
if (!"category" %in% names(trend_driver)) {
  trend_driver <- trend_driver %>%
    left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category",
                             "NEEpartition", "source", "AI")],
              by = "site_ID")
}

climate_vegs <- c("cold_ENF", "cold_DBF", "cold_nonforest", "dry", "wet")
key_vars <- c(
  "CO2 + TA_gs_at",
  "CO2 + TA_ngs_at",
  "SWC_gs_at",
  "CO2 + TA_gs_at + P_gs_at + SWC_gs_at",
  "CO2 + LE_at + DIST_at + TA_at"
)
projection_frame_names <- c(
  cold_ENF = "cold_evergreen_forest",
  cold_DBF = "cold_deciduous_forest",
  cold_nonforest = "cold_non_forest",
  dry = "dry",
  wet = "wet-warm"
)
unique_climate_veg_groups <- data.frame(value = c(2, 1, 3:5), label = climate_vegs)

m_try <- 2
min_node_size <- 5

for (iclimate_veg in seq_along(climate_vegs)) {
  climate_veg <- climate_vegs[iclimate_veg]
  key_var <- key_vars[iclimate_veg]
  message("Exporting ", climate_veg)

  if (climate_veg == "cold_ENF") {
    data_group <- trend_driver %>% filter(IGBP %in% c("ENF", "EBF") & category == "cold")
  } else if (climate_veg == "cold_DBF") {
    data_group <- trend_driver %>% filter(IGBP %in% c("DBF", "DNF", "MF") & category == "cold")
  } else if (climate_veg == "cold_nonforest") {
    data_group <- trend_driver %>% filter(!IGBP %in% c("ENF", "EBF", "DBF", "DNF", "MF") & category == "cold")
  } else if (climate_veg == "dry") {
    data_group <- trend_driver %>% filter(category == "dry")
  } else if (climate_veg == "wet") {
    data_group <- trend_driver %>% filter(category == "wet")
  }

  cells <- which(terra::values(zone_group) == unique_climate_veg_groups$value[unique_climate_veg_groups$label == climate_veg])

  df_project_past <- data.frame(
    cell_id = cells,
    CO2 = cell_values(co2_past, cells),
    TA_gs_at = cell_values(ta_gs_at_past, cells),
    TA_ngs_at = cell_values(ta_ngs_at_past, cells),
    P_gs_at = cell_values(tp_gs_at_past, cells),
    P_ngs_at = cell_values(tp_ngs_at_past, cells),
    LE_gs_at = cell_values(le_gs_at_past, cells),
    SW_gs_at = cell_values(sw_gs_at_past, cells),
    DIST_at = cell_values(DIST_at_wet_past, cells),
    MAP = cell_values(map_baseline, cells),
    MAT = cell_values(mat_baseline, cells),
    NEE = cell_values(nee_baseline, cells),
    GPP = cell_values(gpp_baseline, cells),
    LE = cell_values(le_baseline, cells),
    AI = cell_values(aridity, cells),
    GPP_at_mean = cell_values(gpp_at_baseline, cells),
    climate_id = cell_values(koppen_climate, cells),
    IGBP_id = cell_values(land_cover, cells)
  ) %>%
    left_join(df_koppen_climate, by = "climate_id") %>%
    left_join(df_IGBP, by = "IGBP_id") %>%
    flatten_frame()

  new_climate <- setdiff(unique(df_project_past$Climate.class), unique(data_group$Climate.class))
  new_IGBP <- setdiff(unique(df_project_past$IGBP), unique(data_group$IGBP))
  if (length(new_climate) > 0) {
    df_project_past$Climate.class[df_project_past$Climate.class %in% new_climate] <- names(which.max(table(df_project_past$Climate.class)))
  }
  if (length(new_IGBP) > 0) {
    df_project_past$IGBP[df_project_past$IGBP %in% new_IGBP] <- names(which.max(table(df_project_past$IGBP)))
  }

  df_project_fut <- df_project_past %>%
    mutate(
      CO2 = cell_values(co2_fut, cells),
      TA_gs_at = cell_values(ta_gs_at_fut, cells),
      TA_ngs_at = cell_values(ta_ngs_at_fut, cells),
      P_gs_at = cell_values(tp_gs_at_fut, cells),
      P_ngs_at = cell_values(tp_ngs_at_fut, cells),
      LE_gs_at = cell_values(le_gs_at_fut, cells),
      SW_gs_at = cell_values(sw_gs_at_fut, cells),
      DIST_at = cell_values(DIST_at_wet_fut, cells)
    ) %>%
    flatten_frame()

  if (climate_veg != "wet") {
    df_project_past <- df_project_past %>% filter(if_all(-c(DIST_at), ~ !is.na(.))) %>% flatten_frame()
    df_project_fut <- df_project_fut %>% filter(if_all(-c(DIST_at), ~ !is.na(.))) %>% flatten_frame()
  } else {
    df_project_past <- df_project_past %>%
      filter(if_all(-c(TA_ngs_at, P_ngs_at), ~ !is.na(.))) %>%
      rename(TA_at = TA_gs_at, SW_at = SW_gs_at, P_at = P_gs_at, LE_at = LE_gs_at) %>%
      flatten_frame()
    df_project_fut <- df_project_fut %>%
      filter(if_all(-c(TA_ngs_at, P_ngs_at), ~ !is.na(.))) %>%
      rename(TA_at = TA_gs_at, SW_at = SW_gs_at, P_at = P_gs_at, LE_at = LE_gs_at) %>%
      flatten_frame()
  }

  if (climate_veg %in% c("cold_nonforest", "dry")) {
    rf_model <- ranger(
      formula = paste0("SWC_gs_at ~ TA_gs_at + SW_gs_at + P_gs_at + CO2",
                       " + MAP + MAT + IGBP + Climate.class + AI"),
      data = data_group,
      mtry = m_try,
      min.node.size = min_node_size,
      importance = "permutation"
    )
    df_project_fut$SWC_gs_at <- predict(rf_model, data = df_project_fut)$predictions
  }

  rf_model <- ranger(
    formula = paste0("NEE_at ~ ", key_var, " + MAP + MAT + AI + IGBP + Climate.class + GPP_at_mean"),
    data = data_group,
    mtry = m_try,
    min.node.size = min_node_size
  )
  df_project_fut$NEE_at <- predict(rf_model, data = df_project_fut)$predictions

  out_path <- file.path(data_dir, paste0("proj_fut_", projection_frame_names[[climate_veg]], ".csv"))
  write.csv(df_project_fut, out_path, row.names = FALSE)
  message("  wrote ", out_path, " (", nrow(df_project_fut), " rows)")
}

cat("\nDone exporting applicability-domain projection frames.\n")
