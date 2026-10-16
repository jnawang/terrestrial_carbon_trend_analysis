# Stratified site-level bootstrap with two-year thinning
#
# Each bootstrap iteration treats every climate-vegetation class separately,
# samples sites with replacement within each IGBP x observation-count stratum,
# and preserves the original number of sites in every climate-vegetation x IGBP
# x observation-class combination. Observation classes are 1, 2-3, and 4+
# retained observations after thinning.
# Every sampled site copy receives a unique cluster ID and one randomly selected
# valid year_mid thinning option. Only total effects on NEE_at are fit.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(dagitty)
  library(nlme)
})

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------

project_dir  <- path.expand("~/YaleLab/carbon_trend")
data_file    <- file.path(project_dir, "data2026", "potential_trend_driver.csv")
site_file    <- file.path(project_dir, "data2026", "long_term_ec_sites.csv")
output_dir   <- file.path(project_dir, "sensitivity_analysis")

thin_interval   <- 2
n_boot          <- 2000       # use 10-20 first to test the complete workflow
random_seed     <- 20260918
checkpoint_every <- 25

# Fixed inverse-variance weights. These are calculated once from the complete,
# unthinned dataset and are not recalculated inside the bootstrap loop.
weight_sd_source <- "NEE_at_sd"
weight_nyear_var <- "nyear"
model_weight_var <- ".fixed_weight_sd"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Read and prepare the unthinned data
# -----------------------------------------------------------------------------

trend_driver <- read.csv(data_file, check.names = FALSE)
sites_long   <- read.csv(site_file, check.names = FALSE)


# what if only using 15-year data; if use all data, comment this line out. 
trend_driver <- trend_driver %>% filter(year_end - year_start >= 14)


sites_long <- sites_long %>%
  mutate(source = if_else(source == "AmeriFlux_BASE", "Eddyproc", "Oneflux"))

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
# Step 1: retain the longest forward and backward thinning options for each site
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

# PREVIOUS METHOD: retain every unique forward and backward thinning sequence.
# This is kept for reference but is no longer used.
# make_site_options <- function(site_data, min_gap = 2) {
#   site_id <- unique(site_data$site_ID)
#   years <- sort(unique(site_data$year_mid[is.finite(site_data$year_mid)]))
# 
#   if (length(site_id) != 1) stop("make_site_options received more than one site.")
#   if (length(years) == 0) stop("Site ", site_id, " has no finite year_mid.")
# 
#   groups <- c(
#     directional_groups(years, min_gap, "forward"),
#     directional_groups(years, min_gap, "backward")
#   )
# 
#   keys <- vapply(groups, function(x) paste(x, collapse = "|"), character(1))
#   groups <- groups[!duplicated(keys)]
# 
#   tibble(
#     site_ID = site_id,
#     option_id = seq_along(groups),
#     year_mid_option = groups,
#     row_id_option = map(
#       groups,
#       ~ site_data$.row_id[match(.x, site_data$year_mid)]
#     )
#   )
# }

# NEW METHOD: retain only the sequence with the most observations from each
# direction. If the longest forward and backward sequences are identical, store
# a single option labelled "forward+backward".
make_site_options <- function(site_data, min_gap = 2) {
  site_id <- unique(site_data$site_ID)
  years <- sort(unique(site_data$year_mid[is.finite(site_data$year_mid)]))

  if (length(site_id) != 1) stop("make_site_options received more than one site.")
  if (length(years) == 0) stop("Site ", site_id, " has no finite year_mid.")

  longest_group <- function(direction) {
    groups <- directional_groups(years, min_gap, direction)
    if (length(groups) == 0) return(numeric(0))
    n_observations <- vapply(
      groups,
      function(x) sum(site_data$year_mid %in% x),
      integer(1)
    )
    groups[[which.max(n_observations)]]
  }

  forward <- longest_group("forward")
  backward <- longest_group("backward")

  if (identical(forward, backward)) {
    groups <- list(forward)
    directions <- "forward+backward"
  } else {
    groups <- list(forward, backward)
    directions <- c("forward", "backward")
  }

  tibble(
    site_ID = site_id,
    option_id = seq_along(groups),
    thin_direction = directions,
    n_observations = vapply(
      groups,
      function(x) sum(site_data$year_mid %in% x),
      integer(1)
    ),
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
  select(site_ID, option_id, thin_direction, n_observations, n_years, years)

write.csv(
  option_summary,
  file.path(output_dir, "site_longest_forward_backward_thinning_options.csv"),
  row.names = FALSE
)

site_option_counts <- option_summary %>%
  count(site_ID, name = "n_options")

write.csv(
  site_option_counts,
  file.path(output_dir, "site_longest_forward_backward_option_counts.csv"),
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

# Prepare the original data and valid thinning options separately for each
# climate-vegetation class. This is done once, outside the bootstrap loop.
stratum_data <- imap(analysis_specs, function(spec, climate_veg) {
  x <- spec$subset(trend_driver)
  site_igbp <- x %>% distinct(site_ID, IGBP)

  if (anyNA(site_igbp$site_ID) || anyNA(site_igbp$IGBP)) {
    stop("Missing site_ID or IGBP in climate-vegetation class ", climate_veg, ".")
  }
  if (anyDuplicated(site_igbp$site_ID) > 0) {
    stop("At least one site has multiple IGBP values in ", climate_veg, ".")
  }

  x
})

stratum_options <- map(stratum_data, function(x) {
  x %>%
    filter(!is.na(site_ID), is.finite(year_mid)) %>%
    group_split(site_ID, .keep = TRUE) %>%
    map_dfr(make_site_options, min_gap = thin_interval)
})

# Assign each site to an observation-count class using its mean retained count
# across the longest forward/backward options. In the current data, both
# directions retain the same count for every site, but the mean also handles
# future cases in which the counts differ.
make_observation_registry <- function(data, options, climate_veg) {
  thinning_size <- options %>%
    group_by(site_ID) %>%
    summarise(
      min_thinned_observations = min(n_observations),
      mean_thinned_observations = mean(n_observations),
      max_thinned_observations = max(n_observations),
      .groups = "drop"
    )

  registry <- data %>%
    distinct(site_ID, IGBP) %>%
    left_join(thinning_size, by = "site_ID") %>%
    mutate(
      climate_veg = climate_veg,
      observation_class = case_when(
        mean_thinned_observations <= 1 ~ "1",
        mean_thinned_observations <= 3 ~ "2-3",
        TRUE ~ "4+"
      ),
      observation_class = factor(
        observation_class,
        levels = c("1", "2-3", "4+")
      )
    )

  if (anyNA(registry$mean_thinned_observations) ||
      anyNA(registry$observation_class)) {
    stop("Missing thinning size or observation class in ", climate_veg, ".")
  }

  registry
}

observation_registry <- imap_dfr(stratum_data, function(x, climate_veg) {
  make_observation_registry(
    data = x,
    options = stratum_options[[climate_veg]],
    climate_veg = climate_veg
  )
})

observation_stratum_summary <- observation_registry %>%
  count(climate_veg, IGBP, observation_class, name = "n_sites", .drop = FALSE) %>%
  filter(n_sites > 0)

write.csv(
  observation_registry,
  file.path(output_dir, "bootstrap_site_observation_classes.csv"),
  row.names = FALSE
)

write.csv(
  observation_stratum_summary,
  file.path(output_dir, "bootstrap_IGBP_observation_class_counts.csv"),
  row.names = FALSE
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
    ".row_id", ".source_site_ID", ".bootstrap_draw", "site_ID",
    "IGBP", "Climate.class", "category",
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

sample_sites_within_igbp_and_thin <- function(data, options, climate_veg) {
  site_registry <- make_observation_registry(data, options, climate_veg)

  if (nrow(site_registry) < 2) {
    stop("Fewer than two sites are available in ", climate_veg, ".")
  }
  if (anyNA(site_registry$site_ID) || anyNA(site_registry$IGBP)) {
    stop("Missing site_ID or IGBP in ", climate_veg, ".")
  }
  if (anyDuplicated(site_registry$site_ID) > 0) {
    stop("At least one site has multiple IGBP values in ", climate_veg, ".")
  }

  # Within every IGBP x observation-class stratum, draw exactly the original
  # number of sites with replacement. The proportions of singleton, moderately
  # repeated, and well-repeated sites therefore remain fixed in every iteration.
  draw_table <- site_registry %>%
    group_by(IGBP, observation_class) %>%
    group_modify(~ tibble(
      source_site_ID = sample(.x$site_ID, size = nrow(.x), replace = TRUE),
      draw_within_stratum = seq_len(nrow(.x))
    )) %>%
    ungroup() %>%
    mutate(bootstrap_draw = row_number())

  sampled_copies <- pmap(
    draw_table,
    function(IGBP, observation_class, source_site_ID,
             draw_within_stratum, bootstrap_draw) {
      source_options <- options %>% filter(site_ID == source_site_ID)
      if (nrow(source_options) == 0) {
        stop("No thinning option exists for sampled site ", source_site_ID, ".")
      }

      # Independently choose one valid two-year thinning pattern for this
      # sampled copy: longest forward or longest backward. Identical forward
      # and backward sequences were deduplicated when options were constructed.
      selected_option <- source_options %>% slice_sample(n = 1)
      selected_row_ids <- selected_option$row_id_option[[1]]
      site_copy <- data[match(selected_row_ids, data$.row_id), , drop = FALSE]

      if (anyNA(site_copy$.row_id)) {
        stop("A sampled row ID was not found for site ", source_site_ID, ".")
      }
      if (any(site_copy$IGBP != IGBP)) {
        stop("IGBP changed within sampled site ", source_site_ID, ".")
      }

      site_copy$.source_site_ID <- source_site_ID
      site_copy$.bootstrap_draw <- bootstrap_draw
      site_copy$.bootstrap_IGBP <- IGBP
      site_copy$.bootstrap_observation_class <- as.character(observation_class)

      # Repeated draws of the same source site must be separate clusters in
      # the random effect and temporal correlation structure.
      site_copy$site_ID <- paste0(
        climate_veg, "__", IGBP, "__n", as.character(observation_class),
        "__", source_site_ID,
        "__bootstrap_",
        sprintf("%03d", bootstrap_draw)
      )

      site_copy
    }
  )

  selected_data <- bind_rows(sampled_copies)

  if (n_distinct(selected_data$site_ID) != nrow(site_registry)) {
    stop("The bootstrap sample does not contain the expected number of site copies.")
  }

  list(
    data = selected_data,
    draw_table = draw_table,
    n_sampled_sites = nrow(draw_table),
    n_unique_source_sites = n_distinct(draw_table$source_site_ID),
    n_igbp_strata = n_distinct(draw_table$IGBP),
    n_igbp_observation_strata = nrow(
      distinct(draw_table, IGBP, observation_class)
    )
  )
}

fit_one_bootstrap <- function(iteration, stratum_data, stratum_options,
                              specs, model_specs) {
  map_dfr(names(specs), function(climate_veg) {
    sampled <- sample_sites_within_igbp_and_thin(
      data = stratum_data[[climate_veg]],
      options = stratum_options[[climate_veg]],
      climate_veg = climate_veg
    )
    analysis_data <- scale_for_analysis(sampled$data)

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
            weight_sd_var = model_weight_var
          )
          list(estimate = estimate_i, error = NA_character_)
        },
        error = function(e) {
          list(estimate = NA_real_, error = conditionMessage(e))
        }
      )

      tibble(
        bootstrap_id = iteration,
        n_sampled_sites = sampled$n_sampled_sites,
        n_unique_source_sites = sampled$n_unique_source_sites,
        n_igbp_strata = sampled$n_igbp_strata,
        n_igbp_observation_strata = sampled$n_igbp_observation_strata,
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
    stratum_data = stratum_data,
    stratum_options = stratum_options,
    specs = analysis_specs,
    model_specs = model_specs
  )

  if (bootstrap_id %% checkpoint_every == 0 || bootstrap_id == n_boot) {
    checkpoint <- bind_rows(bootstrap_results[seq_len(bootstrap_id)])
    # saveRDS(
    #   checkpoint,
    #   file.path(output_dir, "traditional_site_bootstrap_total_effect_checkpoint.rds")
    # )
  }
}

bootstrap_results <- bind_rows(bootstrap_results)

# The requested compact result: one estimate for every fitted effect/iteration.
bootstrap_estimates <- bootstrap_results %>%
  select(
    bootstrap_id, n_sampled_sites, n_unique_source_sites, n_igbp_strata,
    n_igbp_observation_strata,
    climate_veg, exposure, outcome, estimate
  )

bootstrap_failures <- bootstrap_results %>%
  filter(!is.na(error)) %>%
  select(
    bootstrap_id, n_sampled_sites, n_unique_source_sites, n_igbp_strata,
    n_igbp_observation_strata,
    climate_veg, exposure, outcome, error
  )

# write.csv(
#   bootstrap_estimates,
#   file.path(output_dir, "total_effect_traditional_site_bootstrap_estimates.csv"),
#   row.names = FALSE
# )
# 
# write.csv(
#   bootstrap_failures,
#   file.path(output_dir, "total_effect_traditional_site_bootstrap_failures.csv"),
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
  file.path(
    output_dir,
    paste0(
      "total_effect_stratified_site_IGBP_observation_class_bootstrap_longest_",
      "forward_backward_interval", thin_interval, "_15years.csv"
    )
  ),
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
#       bootstrap_design = paste(
#         "within each climate-vegetation class, sites sampled with replacement",
#         "within IGBP x observation-class strata; repeated copies assigned",
#         "unique cluster IDs;",
#         "one thinning option sampled per site copy"
#       )
#     ),
#     site_options = site_options,
#     estimates = bootstrap_estimates,
#     failures = bootstrap_failures,
#     summary = bootstrap_summary
#   ),
#   file.path(output_dir, "total_effect_traditional_site_bootstrap_complete.rds")
# )

message("Finished. Results were written to: ", output_dir)
