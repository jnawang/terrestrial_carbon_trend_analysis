# estimate drivers of NEE trend of each climate-vegetation class
# August 11, 2026; Junna Wang

rm(list=ls())
# change this location where the raw data will be saved accordingly
# dir_rawdata <- '/Users/jw2946/Documents/stability_project/SiteData/'
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/SiteData'
dir_rawdata <- '/Volumes/WZZ_disk/Carbon_trend/SiteData'
#
library(librarian)
shelf(tidyverse, ggpubr, lme4, dagitty, lmerTest, brms, randomForest, glmmTMB, performance, MuMIn)


trend_driver <- read.csv(file.path("data2026", "potential_trend_driver.csv"))
sites_long <- read.csv(file.path('data2026', 'long_term_ec_sites.csv'))

sites_long <- sites_long %>% mutate(source = ifelse(source == "AmeriFlux_BASE", "Eddyproc", "Oneflux"))

# adding IGBP, climate, and category into trend driver
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source")], by = "site_ID")
}

#
trend_driver_site <- trend_driver %>% group_by(site_ID) %>% summarise(NEE_at_sd_min = min(NEE_at_sd, na.rm = T))
trend_driver <- left_join(trend_driver, trend_driver_site, by = "site_ID") %>% mutate(NEE_at_sd = NEE_at_sd / (NEE_at_sd_min * sqrt(trend_driver$nyear)))
variable_weight <- "NEE_at_sd"
#
# trend_driver$nyear_sd <- 1 / sqrt(trend_driver$nyear)
# variable_weight <- "nyear_sd"

# variable_weight <- NULL

# settings for the thinning sensitivity analysis
thin_interval <- 2
n_boot <- 2000
random_seed <- 20260918

#--------------------------------------Define DAGs------------------------------
# define and validate DAGs once, before thinning and bootstrapping
DAG_cold_evergreen_forest <- dagitty("dag {
  
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
DAG_cold_deciduous_forest <- dagitty("dag {
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
DAG_cold_non_forest <- dagitty("dag {
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

  DIST_at  -> TA_ngs_at
  
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

dag_list <- list(cold_evergreen_forest = DAG_cold_evergreen_forest,
                 cold_deciduous_forest = DAG_cold_deciduous_forest,
                 cold_non_forest = DAG_cold_non_forest,
                 dry = DAG_dry, `wet-warm` = DAG_wet)
dag_acyclic <- tibble(climate_veg = names(dag_list),
                      is_acyclic = vapply(dag_list, isAcyclic, logical(1)))
print(dag_acyclic)
if (!all(dag_acyclic$is_acyclic)) stop("At least one DAG is cyclic; inspect 'dag_acyclic'.")

#----------------------------------test data-DAG consistency--------------------
# test data-DAG consistency once using the original, unthinned data
test_dag_consistency <- function(dag, data, climate_veg) {
  data_test <- data %>% dplyr::select(-any_of(c("site_ID", "IGBP", "Climate.class",
                                                "category", "NEEpartition", "source", ".row_id")))
  r <- localTests(x = dag, data = data_test, type = "cis", abbreviate.names = FALSE)
  if (nrow(r) > 0) {
    plotLocalTestResults(r)
    message("Data-DAG consistency test: ", climate_veg)
    print(subset(r, p.value < 0.05))
  }
  invisible(r)
}

data_cold_evergreen_forest_original <- trend_driver %>%
  filter(IGBP %in% c("ENF", "EBF"), category == "cold")
data_cold_deciduous_forest_original <- trend_driver %>%
  filter(IGBP %in% c("DBF", "DNF", "MF"), category == "cold", !is.na(DIST_at), !is.na(Lgs_at))
data_cold_non_forest_original <- trend_driver %>%
  filter(!IGBP %in% c("ENF", "EBF", "DBF", "DNF", "MF"), category == "cold")
data_dry_original <- trend_driver %>%
  filter(category == "dry", !is.na(LE_gs_at), !is.na(P_ngs_at))
data_wet_original <- trend_driver %>%
  filter(category == "wet", !is.na(LE_gs_at))

dag_test_cold_evergreen_forest <- test_dag_consistency(DAG_cold_evergreen_forest, data_cold_evergreen_forest_original, "cold_evergreen_forest")
dag_test_cold_deciduous_forest <- test_dag_consistency(DAG_cold_deciduous_forest, data_cold_deciduous_forest_original, "cold_deciduous_forest")
dag_test_cold_non_forest <- test_dag_consistency(DAG_cold_non_forest, data_cold_non_forest_original, "cold_non_forest")
dag_test_dry <- test_dag_consistency(DAG_dry, data_dry_original, "dry")
dag_test_wet <- test_dag_consistency(DAG_wet, data_wet_original, "wet-warm")


##----------------prepare for stratified site bootstrap and effect calculation----------------
# Generate two-year thinning sequences, then retain only the longest sequence
# from each direction for every site.
directional_groups <- function(years, min_gap = 2, direction = c("forward", "backward")) {
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
  if(length(site_id) != 1) stop("make_site_options received more than one site.")
  if(length(years) == 0) stop("Site ",site_id," has no finite year_mid.")

  longest_group <- function(direction) {
    groups <- directional_groups(years,min_gap,direction)
    if(length(groups) == 0) return(numeric(0))
    n_observations <- vapply(groups,function(x) sum(site_data$year_mid %in% x),integer(1))
    groups[[which.max(n_observations)]]
  }

  forward <- longest_group("forward")
  backward <- longest_group("backward")
  if(identical(forward,backward)) {
    groups <- list(forward)
    directions <- "forward+backward"
  } else {
    groups <- list(forward,backward)
    directions <- c("forward","backward")
  }

  tibble(site_ID=site_id,option_id=seq_along(groups),thin_direction=directions,
         n_observations=vapply(groups,function(x) sum(site_data$year_mid %in% x),integer(1)),
         years=groups)
}

trend_driver$.row_id <- seq_len(nrow(trend_driver))
trend_driver_all <- trend_driver

# Keep the five climate-vegetation groups separate throughout bootstrap sampling.
climate_group_data <- list(
  cold_evergreen_forest = trend_driver_all %>%
    filter(IGBP %in% c("ENF", "EBF"), category == "cold"),
  cold_deciduous_forest = trend_driver_all %>%
    filter(IGBP %in% c("DBF", "DNF", "MF"), category == "cold", !is.na(DIST_at), !is.na(Lgs_at)),
  cold_non_forest = trend_driver_all %>%
    filter(!IGBP %in% c("ENF", "EBF", "DBF", "DNF", "MF"), category == "cold"),
  dry = trend_driver_all %>%
    filter(category == "dry", !is.na(LE_gs_at), !is.na(P_ngs_at)),
  `wet-warm` = trend_driver_all %>%
    filter(category == "wet", !is.na(LE_gs_at))
)

climate_group_options <- map(climate_group_data, function(data) {
  data %>% filter(!is.na(site_ID), is.finite(year_mid)) %>%
    group_split(site_ID, .keep = TRUE) %>%
    map_dfr(make_site_options, min_gap = thin_interval)
})

# A site's observation class is fixed before bootstrapping, using its mean
# retained count across its valid forward/backward thinning options.
make_observation_registry <- function(data, options, climate_veg) {
  site_igbp <- data %>% distinct(site_ID, IGBP)
  if(anyNA(site_igbp$site_ID) || anyNA(site_igbp$IGBP))
    stop("Missing site_ID or IGBP in ", climate_veg, ".")
  if(any(duplicated(site_igbp$site_ID)))
    stop("At least one site has multiple IGBP values in ", climate_veg, ".")

  thinning_size <- options %>% group_by(site_ID) %>% summarise(
    min_thinned_observations = min(n_observations),
    mean_thinned_observations = mean(n_observations),
    max_thinned_observations = max(n_observations), .groups = "drop"
  )

  registry <- site_igbp %>% left_join(thinning_size, by = "site_ID") %>% mutate(
    climate_veg = climate_veg,
    observation_class = case_when(
      mean_thinned_observations <= 1 ~ "1",
      mean_thinned_observations <= 3 ~ "2-3",
      TRUE ~ "4+"
    ),
    observation_class = factor(observation_class, levels = c("1", "2-3", "4+"))
  )
  if(anyNA(registry$mean_thinned_observations) || anyNA(registry$observation_class))
    stop("Missing thinning size or observation class in ", climate_veg, ".")
  registry
}

observation_registry <- imap(climate_group_data, function(data, climate_veg) {
  make_observation_registry(data, climate_group_options[[climate_veg]], climate_veg)
})

observation_class_summary <- bind_rows(observation_registry) %>%
  count(climate_veg, IGBP, observation_class, name = "n_sites", .drop = FALSE) %>%
  filter(n_sites > 0)
# write.csv(bind_rows(observation_registry),
#           file.path("data2026", "bootstrap_site_observation_classes_interval2.csv"), row.names = FALSE)
# write.csv(observation_class_summary,
#           file.path("data2026", "bootstrap_IGBP_observation_class_counts_interval2.csv"), row.names = FALSE)

# Traditional clustered bootstrap: within each IGBP x observation-class stratum,
# draw the original number of sites with replacement. A sampled site contributes
# one randomly selected valid thinned time series. Repeated draws receive unique
# site_ID values so mixed models treat them as independent bootstrap clusters.
sample_stratified_sites <- function(data, options, registry, climate_veg) {
  draw_table <- registry %>% group_by(IGBP, observation_class) %>%
    group_modify(~tibble(
      source_site_ID = sample(.x$site_ID, size = nrow(.x), replace = TRUE),
      draw_within_stratum = seq_len(nrow(.x))
    )) %>% ungroup() %>% mutate(bootstrap_draw = row_number())

  sampled_copies <- pmap(draw_table, function(IGBP, observation_class, source_site_ID,
                                               draw_within_stratum, bootstrap_draw, ...) {
    source_options <- options %>% filter(site_ID == source_site_ID)
    if(nrow(source_options) == 0) stop("No thinning option for site ", source_site_ID, ".")
    selected_option <- source_options %>% slice_sample(n = 1)
    selected_years <- selected_option$years[[1]]
    site_copy <- data[data$site_ID == source_site_ID & data$year_mid %in% selected_years,,drop = FALSE]
    if(nrow(site_copy) == 0) stop("No selected rows for site ", source_site_ID, ".")
    site_copy$site_ID <- paste(climate_veg, IGBP, as.character(observation_class),
                              source_site_ID, bootstrap_draw, sep = "__")
    site_copy
  })
  bind_rows(sampled_copies)
}

df_total_effect <- NULL
df_direct_effect <- NULL


############################a function to estimate total and direct effects
run_adjusted_model_lme <- function(
    dag,
    exposure,
    outcome,
    data,
    group_vars = "site_ID",
    confounder = NULL,
    optional_fixed = NULL,
    time_var = "year_mid",
    
    weight_sd_var = NULL,  # NEW: column containing row-level SD
    
    use_cor = TRUE,
    cor_type = c("CAR1", "AR1"),
    method = "REML",
    effect_type = c("direct", "total")
) {
  
  effect_type <- match.arg(effect_type)
  cor_type <- match.arg(cor_type)
  
  if (length(group_vars) != 1) {
    stop("nlme::lme version allows exactly one random effect group.")
  }
  
  adjustment_set <- if (effect_type == "direct") {
    adjustmentSets(
      dag,
      exposure = exposure,
      outcome = outcome,
      type = "minimal",
      effect = "direct"
    )
  } else {
    adjustmentSets(
      dag,
      exposure = exposure,
      outcome = outcome,
      type = "minimal",
      effect = "total"
    )
  }
  
  adjust_vars <- if (length(adjustment_set) == 0) {
    character(0)
  } else {
    unlist(adjustment_set[[1]])
  }
  
  optional_fixed <- optional_fixed[
    vapply(optional_fixed, function(v) {
      v %in% names(data) &&
        dplyr::n_distinct(data[[v]], na.rm=TRUE) > 1
    }, logical(1))
  ]

  rhs_terms <- unique(c(
    exposure,
    adjust_vars,
    confounder,
    optional_fixed
  ))
  
  formula_str <- paste(
    outcome,
    "~",
    paste(rhs_terms, collapse = " + ")
  )
  
  random_str <- paste0("(1 | ", group_vars, ")")
  
  cor_str <- if (use_cor) {
    paste0(
      "cor",
      cor_type,
      "(~",
      time_var,
      " | ",
      group_vars,
      ")"
    )
  } else {
    "no autocorrelation"
  }
  
  # NEW: determine whether weighting is requested
  use_weights <- !is.null(weight_sd_var) &&
    length(weight_sd_var) == 1 &&
    nzchar(weight_sd_var)
  
  # NEW: describe weighting in the printed formula
  weight_str <- if (use_weights) {
    paste0(
      "inverse-variance weights from ",
      weight_sd_var,
      "^2"
    )
  } else {
    NULL
  }
  
  formula_full_str <- paste(
    c(
      formula_str,
      random_str,
      cor_str,
      weight_str
    ),
    collapse = " + "
  )
  
  # message(formula_full_str)
  
  data <- data %>%
    dplyr::arrange(
      .data[[group_vars]],
      .data[[time_var]]
    )
  
  # NEW: default is no variance weighting
  weight_struct <- NULL
  
  # NEW: prepare inverse-variance weighting
  if (use_weights) {
    if (!(weight_sd_var %in% names(data))) {
      stop(
        paste(
          "Weight SD column not found:",
          weight_sd_var
        )
      )
    }
    
    # Use a plain data.frame for nlme
    data <- as.data.frame(data)
    
    # Force the SD column to be a plain numeric vector
    data[[weight_sd_var]] <- suppressWarnings(
      as.numeric(as.character(data[[weight_sd_var]]))
    )
    
    # Variables that must be complete before fitting
    model_vars <- unique(c(
      outcome,
      rhs_terms,
      group_vars,
      if (use_cor) time_var,
      weight_sd_var
    ))
    
    missing_vars <- setdiff(model_vars, names(data))
    
    if (length(missing_vars) > 0) {
      stop(
        paste(
          "Model columns not found:",
          paste(missing_vars, collapse = ", ")
        )
      )
    }
    
    # Remove rows that lme would otherwise omit after weights
    # have already been constructed
    data <- data[
      complete.cases(data[,model_vars,drop=FALSE]),
      ,
      drop=FALSE
    ]
    
    # Remove invalid SD values
    data <- data[
      is.finite(data[[weight_sd_var]]) &
        data[[weight_sd_var]] > 0,
      ,
      drop=FALSE
    ]
    
    if (nrow(data) == 0) {
      stop("No rows remain after filtering invalid SD values.")
    }
    
    # Known sampling variance = SD squared
    data$.sampling_variance <-
      unname(data[[weight_sd_var]])^2
    
    # varFixed expects variance, not inverse variance
    weight_struct <- nlme::varFixed(
      ~.sampling_variance
    )
    
    # Reorder after filtering
    data <- data[
      order(
        data[[group_vars]],
        data[[time_var]]
      ),
      ,
      drop=FALSE
    ]
    
    rownames(data) <- NULL
    
    # Verify that weights and rows align
    initialized_weights <- nlme::Initialize(
      weight_struct,
      data = data
    )
    
    if (
      length(nlme::varWeights(initialized_weights)) !=
      nrow(data)
    ) {
      stop(
        "The variance-weight vector does not match the model rows."
      )
    }
  }
  
  cor_struct <- NULL
  
  if (use_cor) {
    cor_formula <- as.formula(
      paste0(
        "~",
        time_var,
        " | ",
        group_vars
      )
    )
    
    cor_struct <- if (cor_type == "CAR1") {
      nlme::corCAR1(form = cor_formula)
    } else {
      nlme::corAR1(form = cor_formula)
    }
  }
  
  fit <- nlme::lme(
    fixed = as.formula(formula_str),
    random = as.formula(
      paste0("~1 | ", group_vars)
    ),
    correlation = cor_struct,
    
    weights = weight_struct,  # NEW: row-level variance structure
    
    data = data,
    method = method,
    
    # NEW: weighted data are already filtered, preventing
    # row/weight misalignment
    na.action = if (use_weights) {
      stats::na.fail
    } else {
      stats::na.omit
    },
    
    control = nlme::lmeControl(
      opt = "optim",
      maxIter = 200,
      msMaxIter = 200,
      niterEM = 50,
      returnObject = TRUE
    )
  )
  
  tt <- summary(fit)$tTable
  
  r2_overall <- cor(
    model.response(model.frame(fit)),
    fitted(fit, level = 1),
    use = "complete.obs"
  )^2
  
  if (!(exposure %in% rownames(tt))) {
    stop(
      paste(
        "Coefficient",
        exposure,
        "not found in model summary"
      )
    )
  }
  
  est <- tt[exposure, "Value"]
  se <- tt[exposure, "Std.Error"]
  df <- tt[exposure, "DF"]
  
  cor_coef <- if (use_cor) {
    as.numeric(
      coef(
        fit$modelStruct$corStruct,
        unconstrained = FALSE
      )
    )
  } else {
    NA_real_
  }
  
  data.frame(
    exposure = exposure,
    outcome = outcome,
    effect_type = effect_type,
    estimate = est,
    se = se,
    p_value = tt[exposure, "p-value"],
    ci_lower_2.5 = est + qt(0.025, df) * se,
    ci_lower_5 = est + qt(0.05, df) * se,
    ci_upper_95 = est + qt(0.95, df) * se,
    ci_upper_97.5 = est + qt(0.975, df) * se,
    autocor_coef = cor_coef,
    r2_without_autocor = r2_overall,
    
    # NEW: record weighting information
    weighted = use_weights,
    weight_sd_variable = if (use_weights) {
      weight_sd_var
    } else {
      NA_character_
    },
    
    formula = formula_full_str,
    row.names = NULL,
    check.names = FALSE
  )
}


run_adjusted_model_lme_bootstrap <- function(dag, exposure, outcome, effect_type, ...) {
  tryCatch({
    result <- run_adjusted_model_lme(dag = dag, exposure = exposure, outcome = outcome,
                                     effect_type = effect_type, ...)
    result$error <- NA_character_
    result
  }, error = function(e) {
    data.frame(exposure = exposure, outcome = outcome, effect_type = effect_type,
               estimate = NA_real_, se = NA_real_, p_value = NA_real_,
               ci_lower_2.5 = NA_real_, ci_lower_5 = NA_real_,
               ci_upper_95 = NA_real_, ci_upper_97.5 = NA_real_,
               autocor_coef = NA_real_, r2_without_autocor = NA_real_,
               weighted = !is.null(list(...)$weight_sd_var),
               weight_sd_variable = ifelse(is.null(list(...)$weight_sd_var), NA_character_, list(...)$weight_sd_var),
               formula = NA_character_, error = conditionMessage(e), check.names = FALSE)
  })
}

##-----------------------------------------start the bootstrap loop--------------
set.seed(random_seed)
for (bootstrap_id in seq_len(n_boot)) {
  message("Bootstrap iteration ", bootstrap_id, " of ", n_boot)
  trend_driver <- imap_dfr(climate_group_data, function(data, climate_veg) {
    sample_stratified_sites(
      data = data,
      options = climate_group_options[[climate_veg]],
      registry = observation_registry[[climate_veg]],
      climate_veg = climate_veg
    )
  })

################## cold evergreen forests in cold climate
data_cold_evergreen_forest <- trend_driver %>% 
  filter(IGBP %in% c('ENF', 'EBF') & category == "cold")

#--------calculate direct effect-----------
# scale all the variables using Z score
for (i in 1:ncol(data_cold_evergreen_forest)) {
  if (!names(data_cold_evergreen_forest)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end", "nyear_sd", "NEE_at_sd")) {
    data_cold_evergreen_forest[, i] <- scale(data_cold_evergreen_forest[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_cold_evergreen_forest)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  # r1 <- run_adjusted_model(dag = DAG_cold_evergreen_forest, exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct",
  #                          data = data_cold_evergreen_forest, group_vars = c("site_ID", "IGBP", "Climate.class", "source"))
  
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_cold_evergreen_forest, data = data_cold_evergreen_forest,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),  
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)  
  #
  df_direct_effect <- rbind(df_direct_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'cold_evergreen_forest', r1))
}
# df_direct_effect

#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "P_gs_at", "SW_gs_at", "LE_gs_at", "SWC_gs_at", "Lgs_at", "TA_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  #
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_cold_evergreen_forest, data = data_cold_evergreen_forest,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)
  #
  df_total_effect <- rbind(df_total_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'cold_evergreen_forest', r1))
}
# df_total_effect

##################################### deciduous forests in cold climate
# something unexpected: CO2 -> Lgs_at; CO2_at -> Lgs_at; CO2 -> LE_gs_at; CO2 -> SWC_gs_at; 
# something unexpected: TA in non-growing season affects almost everything on growing-season metero conditions.
# the broad effects of CO2 and TA really surprised me!
data_cold_deciduous_forest <- trend_driver %>% 
  filter(IGBP %in% c('DBF', "DNF", 'MF') & category == "cold")  %>% filter(!is.na(DIST_at) & !is.na(Lgs_at))


#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_cold_deciduous_forest)) {
  if (!names(data_cold_deciduous_forest)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end", "nyear_sd", "NEE_at_sd")) {
    data_cold_deciduous_forest[, i] <- scale(data_cold_deciduous_forest[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_cold_deciduous_forest)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
   r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_cold_deciduous_forest, data = data_cold_deciduous_forest,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight) 
  
  df_direct_effect <- rbind(df_direct_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'cold_deciduous_forest', r1))
}
# df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "P_gs_at", "SW_gs_at", "LE_gs_at", "SWC_gs_at", "Lgs_at", "TA_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_cold_deciduous_forest, data = data_cold_deciduous_forest,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight) 
  #
  df_total_effect <- rbind(df_total_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'cold_deciduous_forest', r1))
}
# df_total_effect

############################## non-forests in cold climate

data_cold_non_forest <- trend_driver %>% 
  filter(!IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold")

#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_cold_non_forest)) {
  if (!names(data_cold_non_forest)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end", "nyear_sd", "NEE_at_sd")) {
    data_cold_non_forest[, i] <- scale(data_cold_non_forest[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_cold_non_forest)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_cold_non_forest, data = data_cold_non_forest,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)

  df_direct_effect <- rbind(df_direct_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'cold_non_forest', r1))
}
# df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "P_gs_at", "SW_gs_at", "LE_gs_at", "SWC_gs_at", "Lgs_at", "TA_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_cold_non_forest, data = data_cold_non_forest,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)
  
  df_total_effect <- rbind(df_total_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'cold_non_forest', r1))
}
# df_total_effect

##################################### dry climate

data_dry <- trend_driver %>% filter(category == "dry") %>% dplyr::filter(!is.na(LE_gs_at) & !is.na(P_ngs_at))


#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_dry)) {
  if (!names(data_dry)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end", "nyear_sd", "NEE_at_sd")) {
    data_dry[, i] <- scale(data_dry[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_dry)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_dry, data = data_dry,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)
  
  df_direct_effect <- rbind(df_direct_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'dry', r1))
}
# df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "SWC_gs_at", "SW_gs_at", "P_gs_at", "LE_gs_at", "TA_ngs_at", "P_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_dry, data = data_dry,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP", "source"), 
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)
   
  df_total_effect <- rbind(df_total_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'dry', r1))
}
# df_total_effect

######################################## wet climate - actual trend
data_wet <- trend_driver %>% 
  filter(category == "wet") %>% dplyr::filter(!is.na(LE_gs_at))


#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_wet)) {
  if (!names(data_wet)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end", "nyear_sd", "NEE_at_sd")) {
    data_wet[, i] <- scale(data_wet[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_wet)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_wet, data = data_wet,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)
  #
  df_direct_effect <- rbind(df_direct_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'wet-warm', r1))
}
# df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "SWC_at", "SW_at", "P_at", "LE_at", "TA_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme_bootstrap(dag = DAG_wet, data = data_wet,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid", weight_sd_var = variable_weight)
  #
  df_total_effect <- rbind(df_total_effect, data.frame(bootstrap_id = bootstrap_id, climate_veg = 'wet-warm', r1))
}
# df_total_effect

} # end bootstrap loop

## write out results
safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  unname(quantile(x, p, names = FALSE))
}

summarise_bootstrap <- function(x) {
  x %>% group_by(climate_veg, exposure, outcome, effect_type) %>%
    summarise(n_boot = n(),
              n_success = sum(is.finite(estimate)),
              failure_rate = 1 - n_success / n_boot,
              estimate_median = ifelse(n_success > 0, median(estimate, na.rm = TRUE), NA_real_),
              estimate_mean = ifelse(n_success > 0, mean(estimate, na.rm = TRUE), NA_real_),
              bootstrap_sd = ifelse(n_success > 1, sd(estimate, na.rm = TRUE), NA_real_),
              ci_lower_2.5 = safe_quantile(estimate, 0.025),
              ci_upper_97.5 = safe_quantile(estimate, 0.975),
              ci_lower_5 = safe_quantile(estimate, 0.05),
              ci_upper_95 = safe_quantile(estimate, 0.95),
              .groups = "drop")
}

df_total_effect_bootstrap <- summarise_bootstrap(df_total_effect)
df_direct_effect_bootstrap <- summarise_bootstrap(df_direct_effect)

write.csv(df_total_effect_bootstrap, file.path("data2026", "total_effect_bootstrap_95CI.csv"), row.names = FALSE)
write.csv(df_direct_effect_bootstrap, file.path("data2026", "direct_effect_bootstrap_95CI.csv"), row.names = FALSE)
