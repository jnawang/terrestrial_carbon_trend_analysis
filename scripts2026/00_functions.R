library(tidyverse)

data_length_thresh <- 10

# This function organizes different sources of raw flux data; updated on Jan, 2026
# It takes a few mins to run this function
get_rawflux_file_structure <- function(dir_rawdata){
  
  files.unzip.FLUXNET2015 <- list.files(file.path(dir_rawdata, "FLUXNET2015", "unzip_select"), full.names = T)
  files.unzip.FLUXNET2015_shuttle <- list.files(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_connect2015"), full.names = T)
  files.unzip.FLUXNET2025_shuttle <- list.files(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_connect2025"), full.names = T)
  files.unzip.FLUXNET_shuttle <- list.files(file.path(dir_rawdata, "Fluxnet_shuttle", "unzip_select"), full.names = T)
  files.unzip.Ameri <- list.files(file.path(dir_rawdata, "Ameriflux_FLUXNET", "unzip"), pattern =  "_FLUXNET_FULLSET_(HR|HH|DD|YY)", full.names = T)
  files.Ameri.ReddyProcGapFill <- list.files(file.path(dir_rawdata, "Ameriflux_BASE", "ReddyProcGapFill"), full.names=T)
  files.csv.HRBFlux2022 <- list.files(file.path(dir_rawdata, "HRBFlux2022", "csv"), full.names=T)
  #
  df_file_structure <- data.frame(sources = c("FLUXNET2015", "FLUXNET2015_shuttle", "FLUXNET2025_shuttle", 
                                              "FLUXNET_shuttle", "AmeriFlux_FLUXNET", "AmeriFlux_BASE", "HRBFlux2022"), 
                                  folders = I(list(list(files.unzip.FLUXNET2015), list(files.unzip.FLUXNET2015_shuttle), 
                                                   list(files.unzip.FLUXNET2025_shuttle), list(files.unzip.FLUXNET_shuttle),
                                                   list(files.unzip.Ameri), 
                                                   list(files.Ameri.ReddyProcGapFill), list(files.csv.HRBFlux2022))), 
                                  HHfiles = c("_FLUXNET2015_FULLSET_(HR|HH)", "_FLUXNET_FLUXMET_(HR|HH)", 
                                              "_FLUXNET_FLUXMET_(HR|HH)", "_FLUXNET_FLUXMET_(HR|HH)", 
                                              "_FLUXNET_FULLSET_(HR|HH)", 
                                              "_EDDYPROC_FULLSET_HH", "_ReddyProcRF_FULLSET_HH"), 
                                  DDfiles = c("_FLUXNET2015_FULLSET_(DD)", "_FLUXNET_FLUXMET_(DD)", 
                                              "_FLUXNET_FLUXMET_(DD)", "_FLUXNET_FLUXMET_(DD)", 
                                              "_FLUXNET_FULLSET_(DD)",
                                              "_EDDYPROC_FULLSET_DD", "_ReddyProcRF_FULLSET_DD"),
                                  YYfiles = c("_FLUXNET2015_FULLSET_(YY)", "_FLUXNET_FLUXMET_(YY)", 
                                              "_FLUXNET_FLUXMET_(YY)", "_FLUXNET_FLUXMET_(YY)", 
                                              "_FLUXNET_FULLSET_(YY)", 
                                              "_EDDYPROC_FULLSET_YY", "_ReddyProcRF_FULLSET_YY"))
  return(df_file_structure)
}

#---------------functions to calculate absolute or relative slopes--------------
# use new methods covering irregular time series
library(mblm)
cal_sens_slope <- function(time, y) {
  id.nonna <- !is.na(y) & !is.na(time)
  y <- y[id.nonna]
  time <- time[id.nonna]
  #
  model <- mblm(y ~ time, repeated = FALSE)
  # summary(model)
  return(suppressWarnings(summary(model)$coefficients[2, 1]))
}


cal_sens_slope_NEP <- function(time, y, nboot=1000) {
  ok <- complete.cases(time, y)
  time <- time[ok]
  y <- y[ok]
  
  if (length(y) < 3) return(c(slope=NA, slope_sd=NA))
  
  slope <- unname(coef(mblm::mblm(y~time, repeated=FALSE))[2])
  
  set.seed(123)
  boot_slopes <- replicate(nboot, {
    i <- sample.int(length(y), replace=TRUE)
    y_sub <- y[i]
    time_sub <- time[i]
    tryCatch(
      unname(coef(mblm::mblm(y_sub~time_sub, repeated=FALSE))[2]),
      error=function(e) NA_real_
    )
  })
  
  c(slope=slope, slope_sd=sd(boot_slopes, na.rm=TRUE))
}


cal_sens_slope_scale <- function(time, y, yall) {
  id.nonna <- !is.na(y) & !is.na(time)
  y <- y[id.nonna]
  time <- time[id.nonna]
  #
  if (is.null(yall)) {
    ymean <- mean(y)
    ystd <- sd(y)
  } else {
    ymean <- mean(yall, na.rm=T)
    ystd <- sd(yall, na.rm=T)
  }
  # scale y
  y <- (y - ymean) / ystd
  model <- mblm(y ~ time, repeated = FALSE)
  # summary(model)

  result <- suppressWarnings(summary(model)$coefficients[2, c(1, 4)])

  return(as.numeric(result))
}

sites_grow_yearlong <- c('BR-Ma2', 'BR-Sa1', 'US-Elm', 'US-Esm', 'US-LL1', 'US-SP1', 'US-SP3', 
                         'GF-Guy', 'ID-PaD', 'AU-Rob', "AU-Tum", "AU-Cum", "AU-Cow")

# match climate class and its coding ID
df_koppen_climate <- tribble(
  ~Climate.class, ~climate_id,
  "Af", 1,  "Am", 2,  "Aw", 3,
  "BWh", 4, "BWk", 5, "BSh", 6, "BSk", 7,
  "Csa", 8, "Csb", 9, "Csc", 10,
  "Cwa", 11, "Cwb", 12, "Cwc", 13,
  "Cfa", 14, "Cfb", 15, "Cfc", 16,
  "Dsa", 17, "Dsb", 18, "Dsc", 19, "Dsd", 20,
  "Dwa", 21, "Dwb", 22, "Dwc", 23, "Dwd", 24,
  "Dfa", 25, "Dfb", 26, "Dfc", 27, "Dfd", 28,
  "ET", 29, "EF", 30
)


df_IGBP <- tribble(
  ~IGBP, ~IGBP_id,
  "ENF", 1, "EBF", 2, "DNF", 3, 
  "DBF", 4, "MF", 5,  "CSH", 6,
  "OSH", 7, "WSA", 8, "SAV", 9, 
  "GRA", 10, "WET", 11, 
  "Croplands", 12, "Urban", 13, 
  "CropMosaic", 14, "Ice", 15, 
  "Barren", 16, "Water", 17, 
  "Water", 0)


##########check if a column exists in a dataframe
is_col_exist <- function(df, colname) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(colname), length(colname) == 1)
  
  if (!(colname %in% names(df))) {
    return(FALSE)
  }
  
  x <- df[[colname]]
  
  if (all(is.na(x) | x == -9999)) {
    return(FALSE)
  }
  
  TRUE
}


################a function to merge ERA and FLUX DATA
# merge ERA and FLUX data
merge_ERA_FLUX <- function(years2update, col_year, col_flux, col_era) {
  df <- data.frame(year=col_year, flux=col_flux, era=col_era)
  id2update <- df$year %in% years2update
  #
  df.sub <- df[!id2update, ]
  mod <- lm(data=df.sub, flux ~ era)
  #
  df$flux[id2update] <- predict(mod, newdata=df[id2update, ])
  return(df$flux)
}


############ this function calculates day of year for the southern hemisphere using July 1st as doy = 1
southern_hemisphere_doy <- function(dates) {
  year_start <-
    ymd(paste(year(dates), "07-01", sep = "-"))  # July 1 of the same year
  doy <- as.numeric(dates - year_start + 1)  # Compute DOY
  
  # Adjust for dates before July 1
  doy[doy < 1] <-
    doy[doy < 1] + 365 + leap_year(year(dates[doy < 1]))
  
  return(doy)
}



##################### the whole site bootstrap
library(dplyr)
library(purrr)
library(ggplot2)

whole_site_bootstrap <- function(
    data,
    site_col = "site",
    year_col = "year",
    value_col = "value",
    n_boot = 500,
    min_sites_per_year = 10,
    seed = 123
) {
  
  #---------------------------------------------
  # 1. Standardize column names and clean data
  #---------------------------------------------
  dat <- data %>%
    transmute(
      site = as.character(.data[[site_col]]),
      year = .data[[year_col]],
      value = as.numeric(.data[[value_col]])
    ) %>%
    filter(
      !is.na(site),
      !is.na(year),
      is.finite(value)
    ) %>%
    group_by(site, year) %>%
    summarise(
      # Ensures one observation per site per year
      value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  site_ids <- unique(dat$site)
  n_original_sites <- length(site_ids)
  
  if (n_original_sites < 2) {
    stop("At least two sites are required.")
  }
  
  #---------------------------------------------
  # 2. Calculate the observed annual mean
  #---------------------------------------------
  observed <- dat %>%
    group_by(year) %>%
    summarise(
      n_sites = n_distinct(site),
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  #---------------------------------------------
  # 3. Bootstrap entire sites
  #---------------------------------------------
  set.seed(seed)
  
  bootstrap_results <- map_dfr(
    seq_len(n_boot),
    function(iteration) {
      
      # Select the same number of sites as in the original
      # dataset, with replacement
      sampled_sites <- sample(
        site_ids,
        size = n_original_sites,
        replace = TRUE
      )
      
      # Count how many times each site was selected
      site_weights <- tibble(site = sampled_sites) %>%
        count(site, name = "bootstrap_weight")
      
      # Attach the bootstrap weight to every available
      # year from each selected site
      dat %>%
        inner_join(site_weights, by = "site") %>%
        group_by(year) %>%
        summarise(
          bootstrap_mean = weighted.mean(
            value,
            w = bootstrap_weight,
            na.rm = TRUE
          ),
          
          # Number of site copies contributing that year
          bootstrap_n = sum(bootstrap_weight),
          
          .groups = "drop"
        ) %>%
        mutate(iteration = iteration)
    }
  )
  
  #---------------------------------------------
  # 4. Calculate percentile 95% confidence limits
  #---------------------------------------------
  bootstrap_ci <- bootstrap_results %>%
    group_by(year) %>%
    summarise(
      n_boot_effective = sum(is.finite(bootstrap_mean)),
      
      bootstrap_median = median(
        bootstrap_mean,
        na.rm = TRUE
      ),
      
      ci_lower = quantile(
        bootstrap_mean,
        probs = 0.025,
        na.rm = TRUE,
        names = FALSE
      ),
      
      ci_upper = quantile(
        bootstrap_mean,
        probs = 0.975,
        na.rm = TRUE,
        names = FALSE
      ),
      
      .groups = "drop"
    )
  
  #---------------------------------------------
  # 5. Join observed means and confidence limits
  #---------------------------------------------
  annual_summary <- observed %>%
    left_join(bootstrap_ci, by = "year") %>%
    mutate(
      sufficient_sites = n_sites >= min_sites_per_year,
      
      # Do not report a CI for extremely sparse years
      ci_lower = if_else(
        sufficient_sites,
        ci_lower,
        NA_real_
      ),
      
      ci_upper = if_else(
        sufficient_sites,
        ci_upper,
        NA_real_
      )
    )
  
  return(
    list(
      summary = annual_summary,
      bootstrap_samples = bootstrap_results,
      cleaned_data = dat
    )
  )
}


##################################use random forest to gap-fill daily GPP, ER, and NEE
library(ranger)

fit_rf_gapfill <- function(
    data,
    response,
    predictors,
    train_condition,
    gap_condition,
    transform = c("none", "log1p"),
    lower_bound = -Inf,
    upper_bound = Inf,
    num_trees = 1000,
    mtry = NULL,
    min_node_size = 5,
    sample_fraction = 1,
    seed = 123,
    num_threads = max(1, parallel::detectCores() - 1),
    minimum_training_rows = 100
) {
  
  #------------------------------------------------------------
  # 1. Check arguments
  #------------------------------------------------------------
  
  transform <- match.arg(transform)
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  
  if (!response %in% names(data)) {
    stop("Response column not found: ", response)
  }
  
  missing_predictors <- setdiff(predictors, names(data))
  
  if (length(missing_predictors) > 0) {
    stop(
      "Predictor columns not found: ",
      paste(missing_predictors, collapse = ", ")
    )
  }
  
  if (length(train_condition) != nrow(data)) {
    stop(
      "'train_condition' must have one TRUE/FALSE value ",
      "for every row in data."
    )
  }
  
  if (length(gap_condition) != nrow(data)) {
    stop(
      "'gap_condition' must have one TRUE/FALSE value ",
      "for every row in data."
    )
  }
  
  # Treat NA conditions as FALSE
  train_condition[is.na(train_condition)] <- FALSE
  gap_condition[is.na(gap_condition)] <- FALSE
  
  response_original <- data[[response]]
  
  if (!is.numeric(response_original)) {
    stop("The response variable must be numeric.")
  }
  
  #------------------------------------------------------------
  # 2. Identify usable training rows
  #------------------------------------------------------------
  
  predictors_complete <- complete.cases(
    data[, predictors, drop = FALSE]
  )
  
  training_index <-
    train_condition &
    !is.na(response_original) &
    predictors_complete
  
  n_training <- sum(training_index)
  
  if (n_training < minimum_training_rows) {
    stop(
      "Too few usable training observations for ",
      response,
      ": ",
      n_training,
      ". Required minimum: ",
      minimum_training_rows
    )
  }
  
  #------------------------------------------------------------
  # 3. Transform response
  #------------------------------------------------------------
  
  response_for_model <- response_original
  
  if (transform == "log1p") {
    
    if (any(
      response_original[training_index] < 0,
      na.rm = TRUE
    )) {
      stop(
        response,
        " contains negative training values. ",
        "Use transform = 'none', or inspect the data."
      )
    }
    
    response_for_model <- log1p(response_original)
    
    inverse_transform <- function(x) {
      expm1(x)
    }
    
  } else {
    
    inverse_transform <- function(x) {
      x
    }
  }
  
  #------------------------------------------------------------
  # 4. Prepare training table
  #------------------------------------------------------------
  
  model_data <- data[, predictors, drop = FALSE]
  
  model_data$.response_rf <- response_for_model
  
  training_data <- model_data[
    training_index,
    c(".response_rf", predictors),
    drop = FALSE
  ]
  
  if (is.null(mtry)) {
    mtry <- max(
      1,
      floor(sqrt(length(predictors)))
    )
  }
  
  mtry <- min(mtry, length(predictors))
  
  #------------------------------------------------------------
  # 5. Fit random forest
  #------------------------------------------------------------
  
  set.seed(seed)
  
  rf_model <- ranger(
    formula = reformulate(
      predictors,
      response = ".response_rf"
    ),
    data = training_data,
    num.trees = num_trees,
    mtry = mtry,
    min.node.size = min_node_size,
    sample.fraction = sample_fraction,
    importance = "permutation",
    keep.inbag = TRUE,
    write.forest = TRUE,
    seed = seed,
    num.threads = num_threads
  )
  
  #------------------------------------------------------------
  # 6. Identify rows that can be gap-filled
  #------------------------------------------------------------
  
  prediction_index <-
    gap_condition &
    is.na(response_original) &
    predictors_complete
  
  predictions <- rep(
    NA_real_,
    nrow(data)
  )
  
  #------------------------------------------------------------
  # 7. Predict missing values
  #------------------------------------------------------------
  
  if (any(prediction_index)) {
    
    prediction_data <- data[
      prediction_index,
      predictors,
      drop = FALSE
    ]
    
    prediction_transformed <- predict(
      rf_model,
      data = prediction_data,
      num.threads = num_threads
    )$predictions
    
    prediction_original_scale <-
      inverse_transform(prediction_transformed)
    
    prediction_original_scale <- pmax(
      prediction_original_scale,
      lower_bound
    )
    
    prediction_original_scale <- pmin(
      prediction_original_scale,
      upper_bound
    )
    
    predictions[prediction_index] <-
      prediction_original_scale
  }
  
  #------------------------------------------------------------
  # 8. Preserve observed values and fill only selected gaps
  #------------------------------------------------------------
  
  filled_values <- response_original
  
  filled_values[prediction_index] <-
    predictions[prediction_index]
  
  #------------------------------------------------------------
  # 9. Calculate out-of-bag diagnostics
  #------------------------------------------------------------
  
  oob_predictions <- inverse_transform(
    rf_model$predictions
  )
  
  oob_predictions <- pmax(
    oob_predictions,
    lower_bound
  )
  
  oob_predictions <- pmin(
    oob_predictions,
    upper_bound
  )
  
  observed_training <- response_original[
    training_index
  ]
  
  oob_residuals <-
    oob_predictions - observed_training
  
  oob_rmse <- sqrt(
    mean(oob_residuals^2, na.rm = TRUE)
  )
  
  oob_mae <- mean(
    abs(oob_residuals),
    na.rm = TRUE
  )
  
  oob_bias <- mean(
    oob_residuals,
    na.rm = TRUE
  )
  
  oob_r2 <- 1 -
    sum(oob_residuals^2, na.rm = TRUE) /
    sum(
      (
        observed_training -
          mean(observed_training, na.rm = TRUE)
      )^2,
      na.rm = TRUE
    )
  
  oob_scale <- mean(abs(observed_training), na.rm = TRUE)
  
  oob_rel_rmse <- oob_rmse / oob_scale
  oob_rel_mae  <- oob_mae / oob_scale
  oob_rel_bias <- oob_bias / oob_scale
  #------------------------------------------------------------
  # 10. Return results
  #------------------------------------------------------------
  
  list(
    response = response,
    
    model = rf_model,
    
    filled = filled_values,
    
    predicted = predictions,
    
    training_index = training_index,
    
    prediction_index = prediction_index,
    
    n_training = n_training,
    
    n_requested_gaps = sum(
      gap_condition & is.na(response_original)
    ),
    
    n_filled = sum(prediction_index),
    
    n_not_filled = sum(
      gap_condition &
        is.na(response_original) &
        !predictors_complete
    ),
    
    transform = transform,
    
    predictors = predictors,
    
    parameters = list(
      num_trees = num_trees,
      mtry = mtry,
      min_node_size = min_node_size,
      sample_fraction = sample_fraction
    ),
    
    diagnostics = data.frame(
      response = response,
      n_training = n_training,
      n_filled = sum(prediction_index),
      OOB_RMSE = oob_rmse,
      OOB_MAE = oob_mae,
      OOB_bias = oob_bias,
      OOB_rel_RMSE = oob_rel_rmse,
      OOB_rel_MAE = oob_rel_mae,
      OOB_rel_bias = oob_rel_bias,      
      OOB_R2 = oob_r2
    ),
    
    variable_importance = sort(
      rf_model$variable.importance,
      decreasing = TRUE
    )
  )
}


######Calculate the maximum gap length
max_gap_length <- function(x) {
  r <- rle(is.na(x))
  if (!any(r$values)) return(0)
  max(r$lengths[r$values])
}

######Extend gap DOY
extend_gap_doy <- function(gap_doy, n_extend = 7, last_doy = 366) {
  gap_doy <- sort(unique(gap_doy))
  
  grp <- cumsum(c(TRUE, diff(gap_doy) > 1))
  
  intervals <- split(gap_doy, grp)
  
  extended <- unlist(lapply(intervals, function(x) {
    start <- min(x)
    end <- max(x)
    
    start_ext <- if (start == 1) start else max(1, start - n_extend)
    end_ext <- if (end == last_doy) end else min(last_doy, end + n_extend)
    
    start_ext:end_ext
  }))
  
  sort(unique(extended))
}


############ adding fake time series 
fill_gap_years_clim <- function(data, date_col = "TIMESTAMP", gpp_col = "GPP",
                                gap_frac_threshold = 0.8) {
  date_sym <- rlang::sym(date_col); gpp_sym <- rlang::sym(gpp_col)
  
  df <- data %>%
    mutate(.date = as.Date(!!date_sym), .year = lubridate::year(.date),
           .doy = lubridate::yday(.date), .gpp = !!gpp_sym)
  
  df_full <- tibble::tibble(.date = seq(min(df$.date, na.rm = TRUE),
                                        max(df$.date, na.rm = TRUE), by = "day")) %>%
    mutate(.year = lubridate::year(.date), .doy = lubridate::yday(.date)) %>%
    left_join(df %>% select(.date, .gpp), by = ".date")
  
  gap_summary <- df_full %>%
    group_by(.year) %>%
    summarise(n_days = n(), n_missing = sum(is.na(.gpp)),
              gap_frac = n_missing / n_days, .groups = "drop")
  
  gap_years <- gap_summary %>%
    filter(gap_frac >= gap_frac_threshold) %>%
    pull(.year)
  
  clim <- df_full %>%
    filter(!.year %in% gap_years, !is.na(.gpp)) %>%
    group_by(.doy) %>%
    summarise(.gpp_clim = median(.gpp, na.rm = TRUE), .groups = "drop")
  
  out <- df_full %>%
    left_join(clim, by = ".doy") %>%
    mutate(.was_filled = .year %in% gap_years & is.na(.gpp),
           GPP_filled = if_else(.was_filled, .gpp_clim, .gpp)) %>%
    transmute(!!date_col := .date, GPP_filled, was_filled = .was_filled)
  
  list(data = out, gap_years = gap_years,
       gap_summary = gap_summary, climatology = clim)
}


###########
## a function to estimate the start and end thaw season
library(dplyr)
library(lubridate)
library(zoo)
soil_temp_season <- function(data, date_col = "TIMESTAMP", temp_col = "TS",
                             min_obs_frac = 0.95, threshold = 0,
                             window = 3) {
  
  date_sym <- rlang::sym(date_col)
  temp_sym <- rlang::sym(temp_col)
  
  df <- data %>%
    transmute(
      date = as.Date(!!date_sym),
      year = year(date),
      doy = yday(date),
      temp = as.numeric(!!temp_sym)
    )
  
  full_df <- tibble(
    date = seq(min(df$date, na.rm = TRUE), max(df$date, na.rm = TRUE), by = "day")
  ) %>%
    mutate(year = year(date), doy = yday(date)) %>%
    left_join(df %>% select(date, temp), by = "date")
  
  full_df %>%
    group_by(year) %>%
    summarise(
      n_days = n(),
      n_obs = sum(!is.na(temp)),
      obs_frac = n_obs / n_days,
      start_doy = {
        if (obs_frac < min_obs_frac) {
          NA_integer_
        } else {
          roll_mean <- zoo::rollapply(
            temp,
            width = window,
            FUN = mean,
            align = "right",
            fill = NA,
            na.rm = FALSE
          )
          idx <- which(roll_mean > threshold)
          if (length(idx) == 0) NA_integer_ else doy[min(idx) - window + 1]
        }
      },
      end_doy = {
        if (obs_frac < min_obs_frac) {
          NA_integer_
        } else {
          roll_mean <- zoo::rollapply(
            temp,
            width = window,
            FUN = mean,
            align = "right",
            fill = NA,
            na.rm = FALSE
          )
          idx <- which(roll_mean > threshold)
          if (length(idx) == 0) NA_integer_ else doy[max(idx)]
        }
      },
      .groups = "drop"
    )
}


######################################a function to combine SWC measurements#########
combine_swc_layers <- function(
    data,
    n_layers = 3,
    min_qc = 0.8,
    min_overlap = 365,
    min_overlap_years = 2,
    min_r2 = 0.70,
    keep_without_qc = TRUE,
    allow_extrapolation = FALSE
) {
  
  # Validate n_layers
  if (
    length(n_layers) != 1 ||
    !is.finite(n_layers) ||
    n_layers < 1
  ) {
    stop("n_layers must be one positive integer.")
  }
  
  n_layers <- as.integer(n_layers)
  
  # Candidate layers: 1 through n_layers
  candidate_layers <- paste0(
    "SWC_F_MDS_",
    seq_len(n_layers)
  )
  
  available_layers <- intersect(
    candidate_layers,
    names(data)
  )
  
  reference_layer <- "SWC_F_MDS_1"
  
  # Initialize final output columns
  data$SWC <- NA_real_
  data$SWC_source <- NA_character_
  data$SWC_calibrated <- FALSE
  
  # Layer 1 must be available because it is the reference
  if (!(reference_layer %in% available_layers)) {
    return(list(
      data = data,
      diagnostics = data.frame(
        layer = reference_layer,
        role = "reference",
        n_valid = 0L,
        n_overlap = NA_integer_,
        n_overlap_years = NA_integer_,
        intercept = NA_real_,
        slope = NA_real_,
        r_squared = NA_real_,
        rmse = NA_real_,
        accepted = FALSE,
        reason = "Reference layer is absent",
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # Determine calibration year
  if ("year" %in% names(data)) {
    
    calibration_year <- suppressWarnings(
      as.integer(data$year)
    )
    
  } else if ("TIMESTAMP" %in% names(data)) {
    
    calibration_year <- lubridate::year(
      data$TIMESTAMP
    )
    
  } else {
    
    calibration_year <- rep(
      NA_integer_,
      nrow(data)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Apply QC filtering separately to every available layer
  # ---------------------------------------------------------------------------
  
  for (layer in available_layers) {
    
    qc_column <- paste0(layer, "_QC")
    
    values <- suppressWarnings(
      as.numeric(as.character(data[[layer]]))
    )
    
    values[
      !is.finite(values) |
        values == -9999
    ] <- NA_real_
    
    if (qc_column %in% names(data)) {
      
      qc <- suppressWarnings(
        as.numeric(as.character(data[[qc_column]]))
      )
      
      qc[
        !is.finite(qc) |
          qc == -9999
      ] <- NA_real_
      
      # Larger QC values represent better daily coverage
      values[
        is.na(qc) |
          qc < min_qc
      ] <- NA_real_
      
    } else if (!keep_without_qc) {
      
      # Strict option: remove a layer when its QC column is absent
      values[] <- NA_real_
    }
    
    # Replace the original layer with its QC-filtered values
    data[[layer]] <- values
  }
  
  reference <- data[[reference_layer]]
  
  if (all(is.na(reference))) {
    return(list(
      data = data,
      diagnostics = data.frame(
        layer = reference_layer,
        role = "reference",
        n_valid = 0L,
        n_overlap = NA_integer_,
        n_overlap_years = NA_integer_,
        intercept = NA_real_,
        slope = NA_real_,
        r_squared = NA_real_,
        rmse = NA_real_,
        accepted = FALSE,
        reason = "No reliable reference values after QC",
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # ---------------------------------------------------------------------------
  # Rank secondary layers by valid-data abundance
  # ---------------------------------------------------------------------------
  
  secondary_layers <- setdiff(
    available_layers,
    reference_layer
  )
  
  if (length(secondary_layers) > 0) {
    
    secondary_counts <- vapply(
      secondary_layers,
      function(layer) {
        sum(is.finite(data[[layer]]))
      },
      numeric(1)
    )
    
    secondary_layers <- secondary_layers[
      order(secondary_counts, decreasing = TRUE)
    ]
  }
  
  diagnostics <- list(
    data.frame(
      layer = reference_layer,
      role = "reference",
      n_valid = sum(is.finite(reference)),
      n_overlap = NA_integer_,
      n_overlap_years = NA_integer_,
      intercept = NA_real_,
      slope = NA_real_,
      r_squared = NA_real_,
      rmse = NA_real_,
      accepted = TRUE,
      reason = "Reference layer",
      stringsAsFactors = FALSE
    )
  )
  
  accepted_calibrated_columns <- character(0)
  
  # ---------------------------------------------------------------------------
  # Calibrate secondary layers to layer 1
  # ---------------------------------------------------------------------------
  
  for (secondary_layer in secondary_layers) {
    
    secondary <- data[[secondary_layer]]
    
    overlap <- (
      is.finite(reference) &
        is.finite(secondary)
    )
    
    n_overlap <- sum(overlap)
    
    overlap_years <- unique(
      calibration_year[
        overlap &
          !is.na(calibration_year)
      ]
    )
    
    n_overlap_years <- length(overlap_years)
    
    calibrated_column <- paste0(
      secondary_layer,
      "_CAL_TO_1"
    )
    
    data[[calibrated_column]] <- NA_real_
    
    accepted <- FALSE
    reason <- "Calibration not evaluated"
    intercept <- NA_real_
    slope <- NA_real_
    r_squared <- NA_real_
    rmse <- NA_real_
    
    if (n_overlap < min_overlap) {
      
      reason <- paste0(
        "Insufficient overlap: ",
        n_overlap,
        " < ",
        min_overlap
      )
      
    } else if (n_overlap_years < min_overlap_years) {
      
      reason <- paste0(
        "Insufficient overlap years: ",
        n_overlap_years,
        " < ",
        min_overlap_years
      )
      
    } else {
      
      calibration_data <- data.frame(
        reference = reference[overlap],
        secondary = secondary[overlap]
      )
      
      calibration_model <- tryCatch(
        lm(
          reference ~ secondary,
          data = calibration_data
        ),
        error = function(e) NULL
      )
      
      if (is.null(calibration_model)) {
        
        reason <- "Regression failed"
        
      } else {
        
        model_coefficients <- coef(
          calibration_model
        )
        
        intercept <- unname(
          model_coefficients["(Intercept)"]
        )
        
        slope <- unname(
          model_coefficients["secondary"]
        )
        
        r_squared <- summary(
          calibration_model
        )$r.squared
        
        rmse <- sqrt(
          mean(
            residuals(calibration_model)^2,
            na.rm = TRUE
          )
        )
        
        if (!is.finite(slope) || slope <= 0) {
          
          reason <- "Regression slope is not positive"
          
        } else if (
          !is.finite(r_squared) ||
          r_squared < min_r2
        ) {
          
          reason <- paste0(
            "R-squared below threshold: ",
            round(r_squared, 3),
            " < ",
            min_r2
          )
          
        } else {
          
          calibrated <- (
            intercept +
              slope * secondary
          )
          
          # Avoid predictions outside the observed calibration range
          if (!allow_extrapolation) {
            
            secondary_range <- range(
              secondary[overlap],
              na.rm = TRUE
            )
            
            outside_range <- (
              secondary < secondary_range[1] |
                secondary > secondary_range[2]
            )
            
            calibrated[outside_range] <- NA_real_
          }
          
          data[[calibrated_column]] <- calibrated
          
          accepted <- TRUE
          reason <- "Accepted"
          
          accepted_calibrated_columns <- c(
            accepted_calibrated_columns,
            calibrated_column
          )
        }
      }
    }
    
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      layer = secondary_layer,
      role = "secondary",
      n_valid = sum(is.finite(secondary)),
      n_overlap = n_overlap,
      n_overlap_years = n_overlap_years,
      intercept = intercept,
      slope = slope,
      r_squared = r_squared,
      rmse = rmse,
      accepted = accepted,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Coalesce in the required order:
  # layer 1 -> calibrated most-abundant secondary -> remaining layers
  # ---------------------------------------------------------------------------
  
  combine_columns <- c(
    reference_layer,
    accepted_calibrated_columns
  )
  
  data$SWC <- do.call(
    dplyr::coalesce,
    data[combine_columns]
  )
  
  # Record which layer supplied each final daily value
  for (column in combine_columns) {
    
    rows_to_fill <- (
      is.na(data$SWC_source) &
        is.finite(data[[column]])
    )
    
    original_layer <- sub(
      "_CAL_TO_1$",
      "",
      column
    )
    
    data$SWC_source[rows_to_fill] <- original_layer
    
    if (column != reference_layer) {
      data$SWC_calibrated[rows_to_fill] <- TRUE
    }
  }
  
  list(
    data = data,
    diagnostics = dplyr::bind_rows(diagnostics)
  )
}
  

