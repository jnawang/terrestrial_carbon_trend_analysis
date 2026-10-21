# Bootstrap sensitivity analysis for thinning choices
#
# Each bootstrap iteration retains every site and randomly selects one valid
# year_mid thinning option for that site. Only total effects on NEE_at are fit.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(dagitty)
  library(nlme)
  library(lme4)
  library(MuMIn)
  library(lubridate)
})


# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
dir_rawdata <- '/Volumes/WZZ_disk/Carbon_trend/SiteData'

project_dir  <- path.expand("~/YaleLab/carbon_trend")
data_file    <- file.path(project_dir, "data2026", "potential_trend_driver.csv")
site_file    <- file.path(project_dir, "data2026", "long_term_ec_sites.csv")
output_dir   <- file.path(project_dir, "sensitivity_analysis")  # , "SWC_multilayer_sensitivity"

trend_driver <- read.csv(data_file, check.names = FALSE)
sites_long   <- read.csv(site_file, check.names = FALSE)

metadata_vars <- c(
  "site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source"
)

missing_metadata <- setdiff(metadata_vars[-1], names(trend_driver))

if (length(missing_metadata) > 0) {
  site_metadata <- sites_long %>%
    select(any_of(c("site_ID", missing_metadata))) %>%
    distinct(site_ID, .keep_all = TRUE)
  
  trend_driver <- trend_driver %>%
    left_join(site_metadata, by = "site_ID")
}

required_columns <- c("site_ID", "year_mid", "NEE_at")
missing_columns <- setdiff(required_columns, names(trend_driver))

if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

#-----------------------

data_grow_season <- read.csv(file.path(project_dir, 'data2026', 'growing_season_yearly.csv'))
grow_season_site <- data_grow_season %>%
  group_by(site_ID) %>% 
  summarise(sos=floor((mean(sos_DER) + mean(sos_TRS)) / 2), 
            eos=ceiling((mean(eos_DER) + mean(eos_TRS)) / 2))
# get growing seasons for all sites. 
source(file.path(project_dir, "scripts2026", "00_functions.R"))
grow_season_site <- bind_rows(grow_season_site, data.frame(site_ID = sites_grow_yearlong, sos = 1, eos = 365))

df_file_structure <- get_rawflux_file_structure(dir_rawdata)

txt <- "
site_ID, goodyear_start, goodyear_end
FI-Sod, 2000, 2016
AU-ASM, 2010, 2024
AU-Stp, 2008, 2022
CA-Mer, 2010, 2022
JP-Spp, 1999, 2014
JP-Ta2, 2007, 2017
US-NC4, 2009, 2023
US-Syv, 2012, 2015
US-BZB, 2018, 2026
RU-Fyo, 2018, 2025
US-EML, 2015, 2025
US-Wkg, 2005, 2018
BE-Dor, 2030, 2030
"
# we should remove BE-Dor, that is why its good years are unrealistic

sites_SWC_trunck <- read.csv(text = txt)

# The selected FLUXNET_shuttle CSVs can contain only layer 1 even when the
# corresponding full FLUXMET file has additional layers.
find_full_swc_file <- function(site_ID, source_name, selected_file, raw_root) {
  if (source_name == "FLUXNET_shuttle") {
    folder <- file.path(raw_root, "Fluxnet_shuttle", "unzip")
    pattern <- paste0("_", site_ID, "_FLUXNET_FLUXMET_DD_")
  } else if (source_name == "FLUXNET2015") {
    folder <- file.path(raw_root, "FLUXNET2015", "unzip")
    pattern <- paste0("_", site_ID, "_FLUXNET2015_FULLSET_DD_")
  } else {
    return(selected_file)
  }
  candidates <- list.files(folder, full.names = TRUE)
  candidates <- candidates[grepl(pattern, basename(candidates), fixed = TRUE) &
                             grepl("\\.csv$", candidates)]
  if (length(candidates) == 1L) candidates else selected_file
}

read_swc_depth_metadata <- function(bif_file, available_layers) {
  empty <- data.frame(layer = character(), depth_m = numeric())
  if (length(bif_file) != 1L || !file.exists(bif_file)) return(empty)
  bif <- read.csv(bif_file, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("GROUP_ID", "VARIABLE", "DATAVALUE")
  if (!all(required %in% names(bif))) return(empty)
  var_rows <- bif[bif$VARIABLE == "VAR_INFO_VARNAME" &
                    bif$DATAVALUE %in% available_layers,
                  c("GROUP_ID", "DATAVALUE")]
  height_rows <- bif[bif$VARIABLE == "VAR_INFO_HEIGHT",
                     c("GROUP_ID", "DATAVALUE")]
  names(var_rows)[2] <- "layer"
  names(height_rows)[2] <- "height_m"
  out <- merge(var_rows, height_rows, by = "GROUP_ID")
  out$height_m <- suppressWarnings(as.numeric(out$height_m))
  out <- out[is.finite(out$height_m) & out$height_m <= 0, ]
  if (nrow(out) == 0L) return(empty)
  out$depth_m <- -out$height_m
  # A layer index can be associated with multiple sensor depths over time.
  # Use its deepest documented position to avoid admitting a partly deep layer.
  out <- aggregate(depth_m ~ layer, out, max)
  out <- out[order(as.integer(sub("SWC_F_MDS_", "", out$layer))), ]
  out[, c("layer", "depth_m")]
}

swc_excluded_years <- function(site_ID) {
  if (site_ID == "US-CMW") c(2008L, 2016L:2021L) else integer()
}

# Return one averaged slope per trend window, plus the contributing layer slopes.
# A single-layer site retains its original layer-1 behavior. A multi-layer
# site with no depth metadata is omitted rather than guessing which layers
# are shallower than 0.1 m.
compute_swc_window_trends <- function(data, windows, layer_depths, sos, eos,
                                      excluded_years = integer(), min_qc = 0.8,
                                      min_coverage = 0.95, min_years = 10L) {
  available <- grep("^SWC_F_MDS_[0-9]+$", names(data), value = TRUE)
  result <- windows
  result$SWC_at_obs <- NA_real_
  result$n_layers_used <- 0L
  details <- list()
  if (length(available) == 0L || nrow(windows) == 0L ||
      length(sos) != 1L || length(eos) != 1L ||
      !is.finite(sos) || !is.finite(eos)) {
    return(list(summary = result, details = data.frame()))
  }

  if (length(available) == 1L) {
    layers <- available
  } else if (nrow(layer_depths) > 0L) {
    layers <- intersect(available, layer_depths$layer[
      is.finite(layer_depths$depth_m) & layer_depths$depth_m < 0.5])
  } else {
    layers <- character()
  }

  #
  if (length(layers) == 0L) {
    return(list(summary = result, details = data.frame()))
  }
  #
  for (layer in layers) {
    values <- suppressWarnings(as.numeric(data[[layer]]))
    values[!is.finite(values) | values == -9999] <- NA_real_
    qc_column <- paste0(layer, "_QC")
    if (qc_column %in% names(data)) {
      qc <- suppressWarnings(as.numeric(data[[qc_column]]))
      values[!is.finite(qc) | qc < min_qc] <- NA_real_
    }
    annual <- do.call(rbind, lapply(sort(unique(data$year)), function(y) {
      if (is.na(y) || y %in% excluded_years) return(NULL)
      in_season <- data$year == y & data$doy >= sos & data$doy <= eos
      in_season[is.na(in_season)] <- FALSE
      if (!any(in_season)) return(NULL)
      n_valid <- sum(is.finite(values[in_season]))
      mean_swc <- if (n_valid > 0L) mean(values[in_season], na.rm = TRUE) else NA_real_
      data.frame(year = y, SWC_obs = mean_swc,
                 coverage = n_valid / (eos - sos + 1))
    }))
    if (is.null(annual)) next
    annual <- annual[is.finite(annual$SWC_obs) &
                       annual$coverage > min_coverage, ]
    if (nrow(annual) < min_years) next
    for (j in seq_len(nrow(windows))) {
      sub <- annual[annual$year >= windows$year_start[j] &
                      annual$year <= windows$year_end[j], ]
      if (nrow(sub) < min_years) next
      details[[length(details) + 1L]] <- data.frame(
        window_row = j, layer = layer,
        depth_m = {
          match_row <- match(layer, layer_depths$layer)
          if (is.na(match_row)) NA_real_ else layer_depths$depth_m[match_row]
        },
        n_years = nrow(sub), first_year = min(sub$year),
        last_year = max(sub$year),
        slope = cal_sens_slope(sub$year, sub$SWC_obs)
      )
    }
  }
  if (length(details) > 0L) {
    details <- do.call(rbind, details)
    for (j in unique(details$window_row)) {
      slopes <- details$slope[details$window_row == j]
      slopes <- slopes[is.finite(slopes)]
      if (length(slopes) > 0L) {
        result$SWC_at_obs[j] <- mean(slopes)
        result$n_layers_used[j] <- length(slopes)
      }
    }
  } else {
    details <- data.frame()
  }
  list(summary = result, details = details)
}

find_swc_bif_file <- function(site_ID, raw_root) {
  folder <- file.path(raw_root, "Fluxnet_shuttle", "unzip")
  candidates <- list.files(folder, full.names = TRUE)
  candidates <- candidates[
    grepl(paste0("_", site_ID, "_FLUXNET_BIFVARINFO_DD_"),
          basename(candidates), fixed = TRUE) & grepl("\\.csv$", candidates)]
  if (length(candidates) == 1L) candidates else NA_character_
}

trend_driver_SWC <- data.frame()
swc_layer_diagnostics <- list()
swc_source_diagnostics <- list()

for (i in 1:nrow(sites_long)) {
  print(i)
  site_ID <- sites_long$site_ID[i]
  #
  print(site_ID)
  folder <- unlist(df_file_structure$folders[df_file_structure$sources == sites_long$source[i]])
  DDfile <- df_file_structure$DDfiles[df_file_structure$sources == sites_long$source[i]]
  #
  selected <- folder[grepl(pattern = paste0(site_ID, DDfile), folder)]
  if (length(selected) != 1L) {
    warning("Expected one daily source for ", site_ID, "; found ", length(selected))
    next
  }
  swc_file <- find_full_swc_file(site_ID, sites_long$source[i], selected, dir_rawdata)
  a_DD <- read.csv(swc_file, check.names = FALSE)
  a_DD[a_DD==-9999] <- NA
  #
  a_DD$TIMESTAMP <- ymd(a_DD$TIMESTAMP)
  a_DD <- a_DD %>% mutate(year=year(a_DD$TIMESTAMP), doy=yday(a_DD$TIMESTAMP))
  
  # use the south of -10 degree to be consistent with growing season detection.
  if (sites_long$Lat[i] < -10) {
    a_DD$doy <- southern_hemisphere_doy(a_DD$TIMESTAMP)
  }
  available_layers <- grep("^SWC_F_MDS_[0-9]+$", names(a_DD), value = TRUE)
  if (length(available_layers) == 0L) next
  bif_file <- find_swc_bif_file(site_ID, dir_rawdata)
  layer_depths <- read_swc_depth_metadata(bif_file, available_layers)
  swc_source_diagnostics[[length(swc_source_diagnostics) + 1L]] <- data.frame(
    site_ID = site_ID, source_file = swc_file,
    n_available_layers = length(available_layers),
    n_layers_with_depth = nrow(layer_depths),
    shallow_layers = paste(layer_depths$layer[layer_depths$depth_m < 0.1],
                           collapse = ";"),
    # skipped_missing_depth = length(available_layers) > 1L && nrow(layer_depths) == 0L
    skipped_missing_depth = nrow(layer_depths) == 0L
  )

  if (site_ID %in% sites_SWC_trunck$site_ID) {
    goodyear_start <- sites_SWC_trunck$goodyear_start[sites_SWC_trunck$site_ID == site_ID]
    goodyear_end <- sites_SWC_trunck$goodyear_end[sites_SWC_trunck$site_ID == site_ID]
    a_DD <- a_DD %>% filter(between(year, goodyear_start, goodyear_end))
  }
  #
  sos <- grow_season_site$sos[grow_season_site$site_ID == site_ID]
  eos <- grow_season_site$eos[grow_season_site$site_ID == site_ID]
  windows <- trend_driver[trend_driver$site_ID == site_ID,
                          c("site_ID", "year_mid", "nyear", "year_start", "year_end")]
  swc_result <- compute_swc_window_trends(
    a_DD, windows, layer_depths, sos, eos,
    excluded_years = swc_excluded_years(site_ID)
  )
  trend_driver_SWC <- bind_rows(
    trend_driver_SWC,
    swc_result$summary %>% filter(is.finite(SWC_at_obs))
  )
  if (nrow(swc_result$details) > 0L) {
    swc_layer_diagnostics[[length(swc_layer_diagnostics) + 1L]] <-
      swc_result$details %>%
      mutate(site_ID = .env$site_ID,
             year_mid = windows$year_mid[window_row],
             year_start = windows$year_start[window_row],
             year_end = windows$year_end[window_row])
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (length(swc_source_diagnostics) > 0L) {
  write.csv(bind_rows(swc_source_diagnostics),
            file.path(output_dir, "SWC_layer_source_diagnostics.csv"), row.names = FALSE)
}
if (length(swc_layer_diagnostics) > 0L) {
  write.csv(bind_rows(swc_layer_diagnostics),
            file.path(output_dir, "SWC_layer_window_trends.csv"), row.names = FALSE)
}
write.csv(trend_driver_SWC,
          file.path(output_dir, "SWC_window_mean_layer_trends.csv"), row.names = FALSE)
length(unique(trend_driver_SWC$site_ID))
plot(trend_driver_SWC$SWC_at_obs)


trend_driver <- trend_driver %>%
  select(-any_of(c("SWC_at_obs", "n_layers_used"))) %>%
  left_join(trend_driver_SWC %>%
              select(site_ID, year_mid, nyear, year_start, year_end, SWC_at_obs),
            by = c("site_ID", "year_mid", "nyear", "year_start", "year_end"))

mod <- lmer(data = trend_driver,  SWC_gs_at ~ SWC_at_obs + (1|site_ID))
summary(mod)
MuMIn::r.squaredGLMM(mod)
# R2m       R2c
# [1,] 0.1180027 0.5788091
# SWC_at_obs    0.10688    0.01801 307.89503   5.935 7.92e-09 ***


# data_group <- trend_driver %>% filter(category == "cold" & IGBP %in% c("DBF", "DNF", "MF"))
data_group <- trend_driver %>% filter(category == "dry")
plot(data_group$SWC_at, data_group$NEE_at)

plot(data_group$SWC_at_obs, data_group$SWC_gs_at)
mod <- lmer(data = data_group,  SWC_gs_at ~ SWC_at_obs + (1|site_ID))
summary(mod)
MuMIn::r.squaredGLMM(mod)



trend_driver <- trend_driver %>% filter(!is.na(SWC_at_obs)) # %>% mutate(SWC_at = SWC_at_obs, SWC_gs_at = SWC_at_obs)
##------------------calculation below

thin_interval   <- 2
n_boot          <- 1000       # use 10-20 first to test the complete workflow
random_seed     <- 20260918
checkpoint_every <- 25

# Fixed inverse-variance weights. These are calculated once from the complete,
# unthinned dataset and are not recalculated inside the bootstrap loop.
weight_sd_source <- "NEE_at_sd"
weight_nyear_var <- "nyear"
model_weight_var <- ".fixed_weight_sd"

use_weights <- TRUE
  
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read and prepare the unthinned data
# -----------------------------------------------------------------------------
trend_driver <- trend_driver %>%
  mutate(source = if_else(source == "AmeriFlux_BASE", "Eddyproc", "Oneflux"))


trend_driver <- trend_driver %>%
  mutate(
    site_ID = as.character(site_ID),
    year_mid = as.numeric(year_mid),
    .row_id = row_number()
  ) %>%
  arrange(site_ID, year_mid, .row_id)

duplicate_site_year <- trend_driver %>%
  filter(!is.na(site_ID), !is.na(year_mid)) %>%
  count(site_ID, year_mid, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_site_year) > 0) {
  stop(
    "site_ID-year_mid must be unique. Duplicate combinations were found; ",
    "inspect object duplicate_site_year."
  )
}

# Calculate the original combined data-length/slope-SD weights before thinning:
#
# fixed SD = NEE_at_sd / (site minimum NEE_at_sd * sqrt(nyear))
# effective precision weight is proportional to 1 / fixed_SD^2.
weight_columns <- c(weight_sd_source, weight_nyear_var)
missing_weight_columns <- setdiff(weight_columns, names(trend_driver))

if (length(missing_weight_columns) > 0) {
  stop(
    "Missing weighting columns: ",
    paste(missing_weight_columns, collapse = ", ")
  )
}

trend_driver[[weight_sd_source]] <- suppressWarnings(
  as.numeric(as.character(trend_driver[[weight_sd_source]]))
)
trend_driver[[weight_nyear_var]] <- suppressWarnings(
  as.numeric(as.character(trend_driver[[weight_nyear_var]]))
)

trend_driver <- trend_driver %>%
  group_by(site_ID) %>%
  mutate(
    .site_sd_min_full = {
      valid_sd <- .data[[weight_sd_source]][
        is.finite(.data[[weight_sd_source]]) & .data[[weight_sd_source]] > 0
      ]
      if (length(valid_sd) > 0) min(valid_sd) else NA_real_
    },
    "{model_weight_var}" := .data[[weight_sd_source]] /
      (.site_sd_min_full * sqrt(.data[[weight_nyear_var]]))
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Step 1: generate the unique thinning options for each site
# -----------------------------------------------------------------------------

# Make a directional partition. Within a group, an exact min_gap is selected
# whenever available; otherwise the nearest year satisfying a larger gap is used.
# After a group is produced, it is extended toward earlier years using the full
# site record. This reproduces the forward/backward option rules discussed.
directional_groups <- function(years, min_gap = 2,
                               direction = c("forward", "backward")) {
  direction <- match.arg(direction)
  years <- sort(unique(years[is.finite(years)]))

  if (length(years) == 0) return(list())

  remaining <- years
  groups <- list()

  while (length(remaining) > 0) {
    current <- if (direction == "forward") min(remaining) else max(remaining)
    selected <- current

    repeat {
      if (direction == "forward") {
        candidates <- remaining[remaining >= current + min_gap]
        if (length(candidates) == 0) break

        exact <- candidates[candidates == current + min_gap]
        next_year <- if (length(exact) > 0) exact[1] else min(candidates)
      } else {
        candidates <- remaining[remaining <= current - min_gap]
        if (length(candidates) == 0) break

        exact <- candidates[candidates == current - min_gap]
        next_year <- if (length(exact) > 0) exact[1] else max(candidates)
      }

      selected <- c(selected, next_year)
      current <- next_year
    }

    selected <- sort(unique(selected))
    remaining <- setdiff(remaining, selected)

    # Make the option maximal toward the beginning of the site's time series.
    first_year <- min(selected)
    repeat {
      candidates <- years[years <= first_year - min_gap]
      if (length(candidates) == 0) break

      exact <- candidates[candidates == first_year - min_gap]
      previous_year <- if (length(exact) > 0) exact[1] else max(candidates)

      if (previous_year %in% selected) break
      selected <- sort(c(previous_year, selected))
      first_year <- previous_year
    }

    groups[[length(groups) + 1L]] <- selected
  }

  groups
}

make_site_options <- function(site_data, min_gap = 2) {
  site_id <- unique(site_data$site_ID)
  years <- sort(unique(site_data$year_mid[is.finite(site_data$year_mid)]))

  if (length(site_id) != 1) stop("make_site_options received more than one site.")
  if (length(years) == 0) stop("Site ", site_id, " has no finite year_mid.")

  groups <- c(
    directional_groups(years, min_gap, "forward"),
    directional_groups(years, min_gap, "backward")
  )

  # Remove identical forward/backward results while retaining distinct options.
  keys <- vapply(groups, function(x) paste(x, collapse = "|"), character(1))
  groups <- groups[!duplicated(keys)]

  tibble(
    site_ID = site_id,
    option_id = seq_along(groups),
    year_mid_option = groups,
    row_id_option = map(
      groups,
      ~ site_data$.row_id[match(.x, site_data$year_mid)]
    )
  )
}

site_options <- trend_driver %>%
  filter(!is.na(site_ID), is.finite(year_mid)) %>%
  group_split(site_ID, .keep = TRUE) %>%
  map_dfr(make_site_options, min_gap = thin_interval)

# Validate all stored options.
invalid_options <- site_options %>%
  mutate(
    valid_gap = map_lgl(
      year_mid_option,
      ~ length(.x) <= 1 || all(diff(.x) >= thin_interval)
    ),
    missing_row_id = map_lgl(row_id_option, ~ anyNA(.x))
  ) %>%
  filter(!valid_gap | missing_row_id)

if (nrow(invalid_options) > 0) {
  stop("Invalid thinning options were generated; inspect invalid_options.")
}

sites_in_data <- sort(unique(trend_driver$site_ID[!is.na(trend_driver$site_ID)]))
sites_in_options <- sort(unique(site_options$site_ID))

if (!identical(sites_in_data, sites_in_options)) {
  stop("Not every site is represented in site_options.")
}

option_summary <- site_options %>%
  mutate(
    years = map_chr(year_mid_option, ~ paste(.x, collapse = ", ")),
    n_years = map_int(year_mid_option, length)
  ) %>%
  select(site_ID, option_id, n_years, years)

write.csv(
  option_summary,
  file.path(output_dir, "site_year_mid_thinning_options.csv"),
  row.names = FALSE
)

site_option_counts <- option_summary %>%
  count(site_ID, name = "n_options")

write.csv(
  site_option_counts,
  file.path(output_dir, "site_thinning_option_counts.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# DAGs used by the original analysis
# -----------------------------------------------------------------------------

# The three cold-climate vegetation analyses used the same graph structure.
DAG_cold <- dagitty("dag {
  CO2 -> TA_gs_at
  CO2 -> SW_gs_at
  CO2 -> P_gs_at
  CO2 -> LE_gs_at
  CO2 -> SWC_gs_at
  CO2 -> Lgs_at
  CO2 -> NEE_at
  TA_gs_at -> SWC_gs_at
  SW_gs_at -> SWC_gs_at
  P_gs_at -> SWC_gs_at
  LE_gs_at -> SWC_gs_at
  SWC_gs_at -> NEE_at
  P_gs_at -> SW_gs_at
  SW_gs_at -> TA_gs_at
  P_gs_at -> TA_gs_at
  P_gs_at -> NEE_at
  SW_gs_at -> NEE_at
  TA_gs_at -> LE_gs_at <- P_gs_at
  SW_gs_at -> LE_gs_at
  TA_gs_at -> NEE_at
  LE_gs_at -> NEE_at
  TA_gs_at -> Lgs_at -> NEE_at
  SWC_gs_at -> Lgs_at
  P_gs_at -> Lgs_at
  SW_gs_at -> Lgs_at
  LE_gs_at -> Lgs_at
  TA_ngs_at -> NEE_at
  CO2 -> TA_ngs_at
  TA_gs_at -> TA_ngs_at
  P_gs_at -> TA_ngs_at
  SW_gs_at -> TA_ngs_at
  LE_gs_at -> TA_ngs_at
  SWC_gs_at -> TA_ngs_at
  TA_ngs_at -> Lgs_at
  DIST_at -> NEE_at
  CO2 -> DIST_at
  TA_gs_at -> DIST_at
  SW_gs_at -> DIST_at
  P_gs_at -> DIST_at
  DIST_at -> LE_gs_at
  DIST_at -> Lgs_at
  DIST_at -> SWC_gs_at
  DIST_at -> TA_ngs_at
}")

DAG_dry <- dagitty("dag {
  CO2 -> TA_gs_at
  CO2 -> SW_gs_at
  CO2 -> P_gs_at
  CO2 -> LE_gs_at
  CO2 -> SWC_gs_at
  CO2 -> NEE_at
  TA_gs_at -> SWC_gs_at
  SW_gs_at -> SWC_gs_at
  P_gs_at -> SWC_gs_at
  LE_gs_at -> SWC_gs_at
  SWC_gs_at -> NEE_at
  P_gs_at -> SW_gs_at
  SW_gs_at -> TA_gs_at
  P_gs_at -> TA_gs_at
  P_gs_at -> NEE_at
  SW_gs_at -> NEE_at
  TA_gs_at -> LE_gs_at <- P_gs_at
  SW_gs_at -> LE_gs_at
  TA_gs_at -> NEE_at
  LE_gs_at -> NEE_at
  TA_ngs_at -> NEE_at
  CO2 -> TA_ngs_at
  TA_gs_at -> TA_ngs_at
  P_gs_at -> TA_ngs_at
  SW_gs_at -> TA_ngs_at
  P_gs_at -> P_ngs_at
  SWC_gs_at -> TA_ngs_at
  CO2 -> P_ngs_at
  P_ngs_at -> NEE_at
  P_ngs_at -> TA_ngs_at
  P_ngs_at -> SWC_gs_at
  P_ngs_at -> TA_gs_at
  P_ngs_at -> SW_gs_at
  DIST_at -> NEE_at
  CO2 -> DIST_at
  TA_gs_at -> DIST_at
  SW_gs_at -> DIST_at
  P_gs_at -> DIST_at
  P_ngs_at -> DIST_at
  DIST_at -> LE_gs_at
  DIST_at -> SWC_gs_at
  DIST_at -> TA_ngs_at
}")

DAG_wet <- dagitty("dag {
  CO2 -> TA_at
  CO2 -> SW_at
  CO2 -> P_at
  CO2 -> LE_at
  CO2 -> SWC_at
  CO2 -> NEE_at
  TA_at -> SWC_at
  SW_at -> SWC_at
  P_at -> SWC_at
  LE_at -> SWC_at
  SWC_at -> NEE_at
  TA_at -> LE_at <- P_at
  TA_at -> NEE_at
  P_at -> NEE_at
  LE_at -> NEE_at
  SW_at -> NEE_at
  SW_at -> TA_at
  SW_at -> LE_at
  P_at -> SW_at
  P_at -> TA_at
  DIST_at -> NEE_at
  CO2 -> DIST_at
  TA_at -> DIST_at
  P_at -> DIST_at
  SW_at -> DIST_at
  DIST_at -> LE_at
  DIST_at -> SWC_at
}")

if (!all(vapply(list(DAG_cold, DAG_dry, DAG_wet), isAcyclic, logical(1)))) {
  stop("At least one DAG is cyclic.")
}

# -----------------------------------------------------------------------------
# Analysis specifications
# -----------------------------------------------------------------------------

cold_exposures <- c(
  "CO2", "TA_gs_at", "P_gs_at", "SW_gs_at", "LE_gs_at",
  "SWC_gs_at", "Lgs_at", "TA_ngs_at", "DIST_at"
)

analysis_specs <- list(
  cold_evergreen_forest = list(
    dag = DAG_cold,
    exposures = cold_exposures,
    subset = function(x) {
      x %>% filter(IGBP %in% c("ENF", "EBF"), category == "cold")
    }
  ),
  cold_deciduous_forest = list(
    dag = DAG_cold,
    exposures = cold_exposures,
    subset = function(x) {
      x %>%
        filter(
          IGBP %in% c("DBF", "DNF", "MF"),
          category == "cold",
          !is.na(DIST_at),
          !is.na(Lgs_at)
        )
    }
  ),
  cold_non_forest = list(
    dag = DAG_cold,
    exposures = cold_exposures,
    subset = function(x) {
      x %>%
        filter(
          !IGBP %in% c("ENF", "EBF", "DBF", "DNF", "MF"),
          category == "cold"
        )
    }
  ),
  dry = list(
    dag = DAG_dry,
    exposures = c(
      "CO2", "TA_gs_at", "SWC_gs_at", "SW_gs_at", "P_gs_at",
      "LE_gs_at", "TA_ngs_at", "P_ngs_at", "DIST_at"
    ),
    subset = function(x) {
      x %>% filter(category == "dry", !is.na(LE_gs_at), !is.na(P_ngs_at))
    }
  ),
  `wet-warm` = list(
    dag = DAG_wet,
    exposures = c("CO2", "SWC_at", "SW_at", "P_at", "LE_at", "TA_at", "DIST_at"),
    subset = function(x) {
      # Preserves the filter in the original script.
      x %>% filter(category == "wet", !is.na(LE_gs_at))
    }
  )
)

# Cache the first minimal total-effect adjustment set outside the bootstrap loop.
model_specs <- imap_dfr(analysis_specs, function(spec, climate_veg) {
  map_dfr(spec$exposures, function(exposure) {
    sets <- adjustmentSets(
      spec$dag,
      exposure = exposure,
      outcome = "NEE_at",
      type = "minimal",
      effect = "total"
    )

    adjust_vars <- if (length(sets) == 0) character(0) else unlist(sets[[1]])

    tibble(
      climate_veg = climate_veg,
      exposure = exposure,
      outcome = "NEE_at",
      adjust_vars = list(adjust_vars)
    )
  })
})

# -----------------------------------------------------------------------------
# Model helpers
# -----------------------------------------------------------------------------

scale_for_analysis <- function(data) {
  never_scale <- c(
    ".row_id", "site_ID", "IGBP", "Climate.class", "category",
    "NEEpartition", "source", "year_mid", "year_start", "year_end",
    "nyear", "nyear_sd", "NEE_at_sd", ".site_sd_min_full",
    ".fixed_weight_sd"
  )

  variables <- names(data)[
    vapply(data, is.numeric, logical(1)) & !names(data) %in% never_scale
  ]

  for (variable in variables) {
    x <- data[[variable]]
    x_sd <- sd(x, na.rm = TRUE)

    if (is.finite(x_sd) && x_sd > 0) {
      data[[variable]] <- as.numeric(scale(x))
    }
  }

  data
}

retain_usable_optional_terms <- function(data, variables) {
  variables[
    vapply(variables, function(variable) {
      variable %in% names(data) &&
        n_distinct(data[[variable]][!is.na(data[[variable]])]) > 1
    }, logical(1))
  ]
}

fit_total_effect <- function(
    data,
    exposure,
    outcome,
    adjust_vars,
    optional_fixed = c("IGBP", "source"),
    group_var = "site_ID",
    time_var = "year_mid",
    weight_sd_var = ".fixed_weight_sd",
    cor_type = "CAR1",
    method = "REML"
) {
  core_terms <- unique(c(exposure, adjust_vars))
  use_weights <- !is.null(weight_sd_var)
  core_vars <- unique(c(
    outcome,
    core_terms,
    group_var,
    time_var,
    if (use_weights) weight_sd_var
  ))
  missing_vars <- setdiff(core_vars, names(data))

  if (length(missing_vars) > 0) {
    stop("Missing model columns: ", paste(missing_vars, collapse = ", "))
  }

  # First retain rows complete for the scientifically required model variables.
  model_data <- as.data.frame(data[complete.cases(data[, core_vars, drop = FALSE]), ])

  if (nrow(model_data) == 0) stop("No complete rows for the core model variables.")

  if (use_weights) {
    model_data[[weight_sd_var]] <- suppressWarnings(
      as.numeric(as.character(model_data[[weight_sd_var]]))
    )

    model_data <- model_data[
      is.finite(model_data[[weight_sd_var]]) & model_data[[weight_sd_var]] > 0,
      ,
      drop = FALSE
    ]

    if (nrow(model_data) == 0) {
      stop("No rows remain after filtering invalid fixed weights.")
    }
  }

  optional_fixed <- retain_usable_optional_terms(model_data, optional_fixed)

  # Include optional terms only when they remain informative after NA filtering.
  if (length(optional_fixed) > 0) {
    model_data <- model_data[
      complete.cases(model_data[, optional_fixed, drop = FALSE]),
      ,
      drop = FALSE
    ]
    optional_fixed <- retain_usable_optional_terms(model_data, optional_fixed)
  }

  rhs_terms <- unique(c(core_terms, optional_fixed))

  if (nrow(model_data) <= length(rhs_terms) + 2) {
    stop("Too few complete observations for this model.")
  }
  if (n_distinct(model_data[[group_var]]) < 2) {
    stop("Fewer than two sites remain for this model.")
  }
  if (n_distinct(model_data[[exposure]]) < 2) {
    stop("Exposure has no variation in retained rows.")
  }
  if (anyDuplicated(model_data[c(group_var, time_var)]) > 0) {
    stop("The correlation structure requires unique time values within each site.")
  }

  model_data <- model_data[order(model_data[[group_var]], model_data[[time_var]]), ]
  rownames(model_data) <- NULL

  fixed_formula <- reformulate(rhs_terms, response = outcome)
  random_formula <- as.formula(paste0("~1|", group_var))
  cor_formula <- as.formula(paste0("~", time_var, "|", group_var))
  cor_struct <- if (cor_type == "CAR1") {
    corCAR1(form = cor_formula)
  } else {
    corAR1(form = cor_formula)
  }

  weight_struct <- NULL
  if (use_weights) {
    model_data$.sampling_variance <- model_data[[weight_sd_var]]^2

    # Constant weights are equivalent to no weighting. This fallback avoids
    # numerical problems when all retained rows happen to have equal weights.
    if (n_distinct(model_data$.sampling_variance) > 1) {
      weight_struct <- varFixed(~ .sampling_variance)
    }
  }

  fit <- lme(
    fixed = fixed_formula,
    random = random_formula,
    correlation = cor_struct,      
    # correlation = NULL,              # testing line
    weights = weight_struct,
    data = model_data,
    method = method,
    na.action = na.fail,
    control = lmeControl(
      opt = "optim",
      maxIter = 200,
      msMaxIter = 200,
      niterEM = 50,
      returnObject = TRUE
    )
  )

  coefficient_table <- summary(fit)$tTable
  if (!(exposure %in% rownames(coefficient_table))) {
    stop("Coefficient not found: ", exposure)
  }

  unname(coefficient_table[exposure, "Value"])
}

sample_one_option_per_site <- function(data, options) {
  selected_options <- options %>%
    group_by(site_ID) %>%
    slice_sample(n = 1) %>%
    ungroup()

  selected_row_ids <- unlist(selected_options$row_id_option, use.names = FALSE)
  selected_data <- data[match(selected_row_ids, data$.row_id), , drop = FALSE]

  if (anyNA(selected_data$.row_id)) stop("A sampled row ID was not found.")
  if (n_distinct(selected_data$site_ID) != n_distinct(data$site_ID)) {
    stop("A bootstrap iteration failed to retain every site.")
  }

  list(data = selected_data, choices = selected_options)
}

fit_one_bootstrap <- function(iteration, data, options, specs, model_specs) {
  sampled <- sample_one_option_per_site(data, options)

  map_dfr(names(specs), function(climate_veg) {
    analysis_data <- specs[[climate_veg]]$subset(sampled$data) %>%
      scale_for_analysis()

    models <- model_specs %>%
      filter(.data$climate_veg == .env$climate_veg)

    map_dfr(seq_len(nrow(models)), function(i) {
      exposure_i <- models$exposure[i]
      outcome_i <- models$outcome[i]
      adjust_i <- models$adjust_vars[[i]]

      result <- tryCatch(
        {
          estimate_i <- fit_total_effect(
            data = analysis_data,
            exposure = exposure_i,
            outcome = outcome_i,
            adjust_vars = adjust_i,
            
# ON and OFF depending on whether I am testing SWC effects
            optional_fixed = if (climate_veg == "wet-warm") {
              "IGBP"
            } else {
              c("IGBP", "source")
            },

            weight_sd_var = if (use_weights) model_weight_var else NULL
          )
          list(estimate = estimate_i, error = NA_character_)
        },
        error = function(e) {
          list(estimate = NA_real_, error = conditionMessage(e))
        }
      )

      tibble(
        bootstrap_id = iteration,
        climate_veg = climate_veg,
        exposure = exposure_i,
        outcome = outcome_i,
        estimate = result$estimate,
        error = result$error
      )
    })
  })
}

# -----------------------------------------------------------------------------
# Step 2: bootstrap loop
# -----------------------------------------------------------------------------

set.seed(random_seed)
bootstrap_results <- vector("list", n_boot)

for (bootstrap_id in seq_len(n_boot)) {
  message("Bootstrap iteration ", bootstrap_id, " of ", n_boot)

  bootstrap_results[[bootstrap_id]] <- fit_one_bootstrap(
    iteration = bootstrap_id,
    data = trend_driver,
    options = site_options,
    specs = analysis_specs,
    model_specs = model_specs
  )

  if (bootstrap_id %% checkpoint_every == 0 || bootstrap_id == n_boot) {
    checkpoint <- bind_rows(bootstrap_results[seq_len(bootstrap_id)])
    # saveRDS(
    #   checkpoint,
    #   file.path(output_dir, "bootstrap_total_effect_checkpoint.rds")
    # )
  }
}

bootstrap_results <- bind_rows(bootstrap_results)

# The requested compact result: one estimate for every fitted effect/iteration.
bootstrap_estimates <- bootstrap_results %>%
  select(bootstrap_id, climate_veg, exposure, outcome, estimate)

bootstrap_failures <- bootstrap_results %>%
  filter(!is.na(error)) %>%
  select(bootstrap_id, climate_veg, exposure, outcome, error)

# write.csv(
#   bootstrap_estimates,
#   file.path(output_dir, "total_effect_bootstrap_estimates.csv"),
#   row.names = FALSE
# )
# 
# write.csv(
#   bootstrap_failures,
#   file.path(output_dir, "total_effect_bootstrap_failures.csv"),
#   row.names = FALSE
# )

# -----------------------------------------------------------------------------
# Percentile 95% bootstrap confidence intervals
# -----------------------------------------------------------------------------

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  unname(quantile(x, probability, names = FALSE, type = 7))
}

bootstrap_summary <- bootstrap_estimates %>%
  group_by(climate_veg, exposure, outcome) %>%
  summarise(
    n_boot = n(),
    n_success = sum(is.finite(estimate)),
    failure_rate = 1 - n_success / n_boot,
    estimate_median = if_else(
      n_success > 0,
      median(estimate, na.rm = TRUE),
      NA_real_
    ),
    estimate_mean = if_else(
      n_success > 0,
      mean(estimate, na.rm = TRUE),
      NA_real_
    ),
    bootstrap_sd = if_else(
      n_success > 1,
      sd(estimate, na.rm = TRUE),
      NA_real_
    ),
    ci_lower_95 = safe_quantile(estimate, 0.025),
    ci_upper_95 = safe_quantile(estimate, 0.975),
    ci_lower_90 = safe_quantile(estimate, 0.05),
    ci_upper_90 = safe_quantile(estimate, 0.95),
    .groups = "drop"
  )

write.csv(
  bootstrap_summary,
  file.path(output_dir, "total_effect_bootstrap_95CI_SWC_modeled.csv"),
  row.names = FALSE
)

# saveRDS(
#   list(
#     settings = list(
#       thin_interval = thin_interval,
#       n_boot = n_boot,
#       random_seed = random_seed,
#       weighting = paste0(
#         "fixed before thinning: [", weight_sd_source,
#         " / (site minimum ", weight_sd_source,
#         " * sqrt(", weight_nyear_var, "))]^2"
#       ),
#       bootstrap_design = "all sites retained; one thinning option sampled per site"
#     ),
#     site_options = site_options,
#     estimates = bootstrap_estimates,
#     failures = bootstrap_failures,
#     summary = bootstrap_summary
#   ),
#   file.path(output_dir, "total_effect_bootstrap_complete.rds")
# )

message("Finished. Results were written to: ", output_dir)
