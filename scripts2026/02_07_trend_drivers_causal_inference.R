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
# trend_driver <- trend_driver %>% filter(nyear >= 10)

# what if only using sites with >= 15-year data
# sites_15yr <- trend_driver %>% group_by(site_ID) %>% summarise(n = n()) %>% filter(n > 1) %>% pull(site_ID)
# trend_driver <- trend_driver %>% filter(site_ID %in% sites_15yr | nyear >= 13)


# # what if I use relative values
# trend_driver$NEE_at <- trend_driver$NEE_t
# trend_driver$TA_at <- trend_driver$TA_t
# trend_driver$TA_gs_at <- trend_driver$TA_gs_t
# trend_driver$TA_ngs_at <- trend_driver$TA_ngs_t
# trend_driver$P_at <- trend_driver$P_t
# trend_driver$P_gs_at <- trend_driver$P_gs_t
# trend_driver$P_ngs_at <- trend_driver$P_ngs_t
# trend_driver$SW_at <- trend_driver$SW_t
# trend_driver$SW_gs_at <- trend_driver$SW_gs_t
# trend_driver$SW_ngs_at <- trend_driver$SW_ngs_t
# trend_driver$SWC_at <- trend_driver$SWC_t
# trend_driver$SWC_gs_at <- trend_driver$SWC_gs_t
# trend_driver$SWC_ngs_at <- trend_driver$SWC_ngs_t
# trend_driver$VPD_at <- trend_driver$VPD_t
# trend_driver$VPD_gs_at <- trend_driver$VPD_gs_t
# trend_driver$VPD_ngs_at <- trend_driver$VPD_ngs_t
# trend_driver$LE_at <- trend_driver$LE_t
# trend_driver$LE_gs_at <- trend_driver$LE_gs_t
# trend_driver$LE_ngs_at <- trend_driver$LE_ngs_t
# trend_driver$DIST_at <- trend_driver$DIST_t
# trend_driver$Lgs_at <- trend_driver$Lgs_t



# adding IGBP, climate, and category into trend driver
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source")], by = "site_ID")
}

df_total_effect <- data.frame(climate_veg = character(), exposure = character(), outcome = character(), effect_type = character(), 
                              estimate = double(), ci_lower_2.5 = double(), ci_lower_5 = double(), ci_upper_95 = double(), ci_upper_97.5 = double(), 
                              formula = character())

df_direct_effect <- df_total_effect


############################a function to estimate total and direct effects
run_adjusted_model_lme <- function(
    dag, exposure, outcome, data,
    group_vars = "site_ID",
    confounder = NULL,
    optional_fixed = NULL,
    time_var = "year_mid",
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
    adjustmentSets(dag, exposure = exposure, outcome = outcome,
                   type = "minimal", effect = "direct")
  } else {
    adjustmentSets(dag, exposure = exposure, outcome = outcome,
                   type = "minimal", effect = "total")
  }
  
  adjust_vars <- if (length(adjustment_set) == 0) character(0) else unlist(adjustment_set[[1]])
  
  rhs_terms <- unique(c(exposure, adjust_vars, confounder, optional_fixed))
  formula_str <- paste(outcome, "~", paste(rhs_terms, collapse = " + "))
  
  random_str <- paste0("(1 | ", group_vars, ")")
  cor_str <- if (use_cor) paste0("cor", cor_type, "(~ ", time_var, " | ", group_vars, ")") else "no autocorrelation"
  formula_full_str <- paste(formula_str, "+", random_str, "+", cor_str)
  message(formula_full_str)
  
  data <- data %>%
    dplyr::arrange(.data[[group_vars]], .data[[time_var]])
  
  cor_struct <- NULL
  if (use_cor) {
    cor_formula <- as.formula(paste0("~ ", time_var, " | ", group_vars))
    cor_struct <- if (cor_type == "CAR1") {
      nlme::corCAR1(form = cor_formula)
    } else {
      nlme::corAR1(form = cor_formula)
    }
  }
  
  fit <- nlme::lme(
    fixed = as.formula(formula_str),
    random = as.formula(paste0("~ 1 | ", group_vars)),
    correlation = cor_struct,
    data = data,
    method = method,
    na.action = na.omit,
    control = nlme::lmeControl(
      opt = "optim",
      maxIter = 200,
      msMaxIter = 200,
      niterEM = 50,
      returnObject = TRUE
    )
  )
  
  tt <- summary(fit)$tTable
  
  # r2 <- MuMIn::r.squaredGLMM(fit)
  r2_overall <- cor(
    model.response(model.frame(fit)),
    fitted(fit, level = 1),
    use = "complete.obs"
  )^2
  
  
  if (!(exposure %in% rownames(tt))) {
    stop(paste("Coefficient", exposure, "not found in model summary"))
  }
  
  est <- tt[exposure, "Value"]
  se <- tt[exposure, "Std.Error"]
  df <- tt[exposure, "DF"]
  
  cor_coef <- if (use_cor) {
    as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
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
    # r2_marginal <- as.numeric(r2[1, 1]),
    # r2_conditional <- as.numeric(r2[1, 2]),
    r2_without_autocor = r2_overall, 
    formula = formula_full_str,
    row.names = NULL,
    check.names = FALSE
  )
}


################## 
## r DAGs of cold evergreen forests in cold climate
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

plot(graphLayout(DAG_cold_evergreen_forest))

isAcyclic(DAG_cold_evergreen_forest)

# test data-DAG consistency; 
#----------------------------------------------------------------------------------
data_cold_evergreen_forest <- trend_driver %>% 
  filter(IGBP %in% c('ENF', 'EBF') & category == "cold")

data_cold_evergreen_forest_test <- data_cold_evergreen_forest %>% dplyr::select(-c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source"))

r <- localTests(x = DAG_cold_evergreen_forest, data=data_cold_evergreen_forest_test, type='cis', abbreviate.names = F)
if (nrow(r) > 0) {
  plotLocalTestResults(r)
  print(subset(r, p.value < 0.05))
}

#----------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_cold_evergreen_forest)) {
  if (!names(data_cold_evergreen_forest)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end")) {
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
  
  r1 <- run_adjusted_model_lme(dag = DAG_cold_evergreen_forest, data = data_cold_evergreen_forest,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")  
  #
  df_direct_effect <- rbind(df_direct_effect, data.frame(climate_veg = 'cold_evergreen_forest', r1))
}
df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "P_gs_at", "SW_gs_at", "LE_gs_at", "SWC_gs_at", "Lgs_at", "TA_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  #
  r1 <- run_adjusted_model_lme(dag = DAG_cold_evergreen_forest, data = data_cold_evergreen_forest,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")
  #
  df_total_effect <- rbind(df_total_effect, data.frame(climate_veg = 'cold_evergreen_forest', r1))
}
df_total_effect

##################################### 
#### DAGs of deciduous forests in cold climate
# something unexpected: CO2 -> Lgs_at; CO2_at -> Lgs_at; CO2 -> LE_gs_at; CO2 -> SWC_gs_at; 
# something unexpected: TA in non-growing season affects almost everything on growing-season metero conditions.
# the broad effects of CO2 and TA really surprised me!

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

plot(graphLayout(DAG_cold_deciduous_forest))

isAcyclic(DAG_cold_deciduous_forest)

# test data-DAG consistency; 
#----------------------------------------------
data_cold_deciduous_forest <- trend_driver %>% 
  filter(IGBP %in% c('DBF', "DNF", 'MF') & category == "cold")  %>% filter(!is.na(DIST_at))

data_cold_deciduous_forest_test <- data_cold_deciduous_forest %>% dplyr::select(-c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source"))

r <- localTests(x = DAG_cold_deciduous_forest, data=data_cold_deciduous_forest_test, type='cis', abbreviate.names = F)
if (nrow(r) > 0) {
  plotLocalTestResults(r)
  print(subset(r, p.value < 0.05))
}

#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_cold_deciduous_forest)) {
  if (!names(data_cold_deciduous_forest)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end")) {
    data_cold_deciduous_forest[, i] <- scale(data_cold_deciduous_forest[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_cold_deciduous_forest)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  # r1 <- run_adjusted_model(dag = DAG_cold_deciduous_forest, exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct",
  #                          data = data_cold_deciduous_forest, group_vars = c("site_ID", "IGBP", "Climate.class", "source"))
  
  r1 <- run_adjusted_model_lme(dag = DAG_cold_deciduous_forest, data = data_cold_deciduous_forest,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")  
  
  df_direct_effect <- rbind(df_direct_effect, data.frame(climate_veg = 'cold_deciduous_forest', r1))
}
df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "P_gs_at", "SW_gs_at", "LE_gs_at", "SWC_gs_at", "Lgs_at", "TA_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme(dag = DAG_cold_deciduous_forest, data = data_cold_deciduous_forest,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")  
  #
  df_total_effect <- rbind(df_total_effect, data.frame(climate_veg = 'cold_deciduous_forest', r1))
}
df_total_effect

#-------------------------------------------------------------------------------
############################## 
## DAGs of non-forests in cold climate
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

plot(graphLayout(DAG_cold_non_forest))

isAcyclic(DAG_cold_non_forest)

####non-growing season's effects. 
# P_ngs_at -> NEE_at
# SWC_ngs_at -> NEE_at  

# test data-DAG consistency; 
#----------------------------------------------
data_cold_non_forest <- trend_driver %>% 
  filter(!IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold")

data_cold_non_forest_test <- data_cold_non_forest %>% dplyr::select(-c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source"))

r <- localTests(x = DAG_cold_non_forest, data=data_cold_non_forest_test, type='cis', abbreviate.names = F)
if (nrow(r) > 0) {
  plotLocalTestResults(r)
  print(subset(r, p.value < 0.05))
}


#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_cold_non_forest)) {
  if (!names(data_cold_non_forest)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end")) {
    data_cold_non_forest[, i] <- scale(data_cold_non_forest[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_cold_non_forest)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  r1 <- run_adjusted_model_lme(dag = DAG_cold_non_forest, data = data_cold_non_forest,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")

  df_direct_effect <- rbind(df_direct_effect, data.frame(climate_veg = 'cold_non_forest', r1))
}
df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "P_gs_at", "SW_gs_at", "LE_gs_at", "SWC_gs_at", "Lgs_at", "TA_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme(dag = DAG_cold_non_forest, data = data_cold_non_forest,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")
  
  df_total_effect <- rbind(df_total_effect, data.frame(climate_veg = 'cold_non_forest', r1))
}
df_total_effect

#####################################
## r dry climate

#  I have to figure out why LE_ngs_at is so important. 
#  Growing season length is controled by different factors. 
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

plot(graphLayout(DAG_dry))
isAcyclic(DAG_dry)

# test data-DAG consistency; 

#----------------------------------------------
data_dry <- trend_driver %>% filter(category == "dry") %>% dplyr::filter(!is.na(LE_gs_at) & !is.na(Lgs_at))

data_dry_test <- data_dry %>% dplyr::select(-c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source"))

r <- localTests(x = DAG_dry, data=data_dry_test, type='cis', abbreviate.names = F)
if (nrow(r) > 0) {
  plotLocalTestResults(r)
  print(subset(r, p.value < 0.05))
}

#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_dry)) {
  if (!names(data_dry)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end")) {
    data_dry[, i] <- scale(data_dry[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_dry)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  r1 <- run_adjusted_model_lme(dag = DAG_dry, data = data_dry,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")
  
  df_direct_effect <- rbind(df_direct_effect, data.frame(climate_veg = 'dry', r1))
}
df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "TA_gs_at", "SWC_gs_at", "SW_gs_at", "P_gs_at", "LE_gs_at", "TA_ngs_at", "P_ngs_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme(dag = DAG_dry, data = data_dry,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")
   
  df_total_effect <- rbind(df_total_effect, data.frame(climate_veg = 'dry', r1))
}
df_total_effect

################################################################################
########################################
## r wet climate - actual trend}

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
plot(graphLayout(DAG_wet))
isAcyclic(DAG_wet)

# test data-DAG consistency; 
#----------------------------------------------
data_wet <- trend_driver %>% 
  filter(category == "wet") %>% dplyr::filter(!is.na(LE_gs_at))

data_wet_test <- data_wet %>% dplyr::select(-c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source"))

r <- localTests(x = DAG_wet, data=data_wet_test, type='cis', abbreviate.names = F)
if (nrow(r) > 0) {
  plotLocalTestResults(r)
  print(subset(r, p.value < 0.05))
}

#-------------------------------------------------------------------------------------
#----calculate direct effect
# scale all the variables using Z score
for (i in 1:ncol(data_wet)) {
  if (!names(data_wet)[i] %in% c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "source", "year_mid", "year_start", "year_end")) {
    data_wet[, i] <- scale(data_wet[, i])
  }
}

# get all the direct effects
edges <- edges(DAG_wet)
direct_effects <- data.frame(from = edges$v, to = edges$w)

# call functions
for (i in 1:nrow(direct_effects)) {
  r1 <- run_adjusted_model_lme(dag = DAG_wet, data = data_wet,
                               exposure = direct_effects$from[i], outcome = direct_effects$to[i], effect_type = "direct", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")
  #
  df_direct_effect <- rbind(df_direct_effect, data.frame(climate_veg = 'wet-warm', r1))
}
df_direct_effect

#------------------------------------------------------------------------------------------------
#--------------------calculate total effect size of each driver----------------------------------
exposure_vars <- c("CO2", "SWC_at", "SW_at", "P_at", "LE_at", "TA_at", "DIST_at")
for (exposure_var in exposure_vars) {
  r1 <- run_adjusted_model_lme(dag = DAG_wet, data = data_wet,
                               exposure = exposure_var, outcome = "NEE_at", effect_type = "total", optional_fixed = c("IGBP", "source"),
                               cor_type = "CAR1", group_vars = "site_ID", time_var = "year_mid")
  #
  df_total_effect <- rbind(df_total_effect, data.frame(climate_veg = 'wet-warm', r1))
}
df_total_effect

## write out results
write.csv(row.names = F, df_total_effect, file.path('data2026', 'effect_total_NEE_trend_driver.csv'))
write.csv(row.names = F, df_direct_effect, file.path('data2026', 'effect_direct_NEE_trend_driver.csv'))

