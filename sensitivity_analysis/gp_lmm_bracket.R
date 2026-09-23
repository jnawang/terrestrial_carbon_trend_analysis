# ============================================================================
# gp_lmm_bracket.R   (scripts2026/review_response/)
#
# Bracket the random-forest NEP-trend projection with two extrapolation-capable
# models + an honesty mask, and report how much of Fig. 4 survives.
#
#   RF  : random forest (numeric drivers) -- interpolates, CLAMPS out of domain.
#   LIN : linear tangent model on the same drivers -- extrapolates as a straight
#         line (unbounded); the optimistic-responsiveness arm.
#   GP  : Gaussian process (DiceKriging) -- extrapolates toward the prior mean
#         with a predictive SD that GROWS away from training. That SD is the
#         applicability surface (the mask).
#
# Computed on the CHANGE (future - past) per grid cell, matching the manuscript's
# decision to interpret differences, not absolute trends. Mirrors
# export_projection_frames_for_ad.R so the frames, helpers, unseen-level
# imputation and wet-warm renaming match your pipeline exactly; it also builds
# the PAST frame (which the AD export omitted) and predicts past SWC_gs_at the
# same way future SWC_gs_at is predicted.
#
# Run from the repo root (like the other review_response scripts).
# Needs: ranger, DiceKriging (GP arm; falls back to RF-vs-LIN if absent), terra.
#
# Outputs -> output/review_response/bracket/
#   bracket_cell_<group>.csv   per-cell rf/lin/gp change, gp_sd, in_domain, survives
#   bracket_summary.csv        per-group agreement, bracket width, % surviving
#   {gp_sd,bracket_width,survives,rf_change,lin_change,gp_change}_<group>.tif  maps
# ============================================================================

suppressMessages({ library(tidyverse); library(terra); library(ranger) })
have_dk <- requireNamespace("DiceKriging", quietly = TRUE)

source(file.path("scripts2026", "review_response", "00_review_helpers.R"))  # data_dir, df_koppen_climate, df_IGBP
source(file.path("scripts2026", "00_functions.R"))
set.seed(128)

USE_CO2_RATE <- tolower(Sys.getenv("USE_CO2_RATE", unset = "false")) %in% c("1", "true", "yes")
                        # TRUE expects data_dir/co2_rate_2005-2023.tif,
                        # co2_rate_2026-2045.tif, and the CO2_at training column.
outdir <- file.path("output", "review_response", if (USE_CO2_RATE) "bracket_rate" else "bracket")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# helpers (identical to export_projection_frames_for_ad.R)
# ---------------------------------------------------------------------------
sens_slope_est <- function(x) {
  ok <- is.finite(x); y <- x[ok]; time <- seq_along(x)[ok]
  if (length(y) < 2) return(NA_real_)
  ij <- utils::combn(seq_along(y), 2); dt <- time[ij[2,]] - time[ij[1,]]; dy <- y[ij[2,]] - y[ij[1,]]
  stats::median(dy[dt != 0] / dt[dt != 0], na.rm = TRUE)
}
cell_values <- function(r, cells) { v <- terra::values(r); if (is.matrix(v)||is.data.frame(v)) as.vector(v[cells,1]) else as.vector(v[cells]) }
trend_raster <- function(path, transform = identity) { r <- transform(terra::rast(path)); terra::app(r, function(x) ifelse(any(is.na(x)), NA, sens_slope_est(x))) }
flatten_frame <- function(df) { df[] <- lapply(df, function(x) if (!is.null(dim(x))) as.vector(x) else x); as.data.frame(df) }

# ---------------------------------------------------------------------------
# rasters (same block as the export script)
# ---------------------------------------------------------------------------
message("Reading rasters / computing trend rasters ...")
nee_baseline <- terra::rast(file.path(data_dir,"baseline_nee_2005-2021.tif"))
gpp_baseline <- terra::rast(file.path(data_dir,"baseline_gpp_2005-2021.tif"))
et_baseline  <- terra::rast(file.path(data_dir,"baseline_et_2005-2021.tif")); le_baseline <- et_baseline/3600*2.45*1e6
mat_baseline <- terra::rast(file.path(data_dir,"baseline_mat.tif"))
map_baseline <- terra::rast(file.path(data_dir,"baseline_map.tif"))
gpp_at_baseline <- terra::rast(file.path(data_dir,"CEDAR_GPP_Sen_slope_ML.tif"))
land_cover <- terra::rast(file.path(data_dir,"968_Land_Cover_Class_0.25degree_layer2.tif"))
koppen_climate <- terra::rast(file.path(data_dir,"koppen_geiger_0p00833333.tif"))
aridity <- terra::rast(file.path(data_dir,"Global-AI.tif"))
zone_group <- terra::rast(file.path(data_dir,"zones_climate_vegetation_group.tif"))
koppen_climate <- terra::resample(koppen_climate, land_cover, method="near")
aridity <- terra::resample(aridity, land_cover)/10000

co2_past <- terra::rast(file.path(data_dir,"co2_mean_2005-2023.tif"))
co2_fut  <- terra::rast(file.path(data_dir,"co2_mean_2025-2045.tif"))
if (USE_CO2_RATE) {
  co2_past <- terra::rast(file.path(data_dir,"co2_rate_2005-2023.tif"))
  co2_fut  <- terra::rast(file.path(data_dir,"co2_rate_2026-2045.tif"))
}

ta_gs_at_past  <- trend_raster(file.path(data_dir,"t2m_gs_2005_2023.tif"),  function(r) r-273.15)
ta_ngs_at_past <- trend_raster(file.path(data_dir,"t2m_ngs_2005_2023.tif"), function(r) r-273.15)
ta_gs_at_fut   <- trend_raster(file.path(data_dir,"tas_gs_2026_2050.tif"),  function(r) r-273.15)
ta_ngs_at_fut  <- trend_raster(file.path(data_dir,"tas_ngs_2026_2050.tif"), function(r) r-273.15)
tp_gs_at_past  <- trend_raster(file.path(data_dir,"tp_gs_2005_2023.tif"),  function(r) r*1000*24)
tp_ngs_at_past <- trend_raster(file.path(data_dir,"tp_ngs_2005_2023.tif"), function(r) r*1000*24)
tp_gs_at_fut   <- trend_raster(file.path(data_dir,"pr_gs_2026_2050.tif"),  function(r) r*86400)
tp_ngs_at_fut  <- trend_raster(file.path(data_dir,"pr_ngs_2026_2050.tif"), function(r) r*86400)
sw_gs_at_past  <- trend_raster(file.path(data_dir,"ssrd_gs_2005_2023.tif"), function(r) r/3600)
sw_gs_at_fut   <- trend_raster(file.path(data_dir,"rsds_gs_2026_2050.tif"))
et_gs_past <- terra::rast(file.path(data_dir,"et_gs_2005_2023.tif"))
et_at_change <- terra::rast(file.path(data_dir,"et_at_past_fut_change.tif"))
le_gs_at_past <- trend_raster(file.path(data_dir,"et_gs_2005_2023.tif"), function(r) r/3600*2.45*1e6)
le_gs_at_fut  <- le_gs_at_past + et_at_change/3600*2.45*1e6
DIST_at_wet_past <- terra::rast(file.path(data_dir,"DIST_at_wet_past.tif")); DIST_at_wet_fut <- DIST_at_wet_past

trend_layers <- list(ta_gs_at_past=ta_gs_at_past, ta_gs_at_fut=ta_gs_at_fut, ta_ngs_at_past=ta_ngs_at_past, ta_ngs_at_fut=ta_ngs_at_fut,
  tp_gs_at_past=tp_gs_at_past, tp_gs_at_fut=tp_gs_at_fut, tp_ngs_at_past=tp_ngs_at_past, tp_ngs_at_fut=tp_ngs_at_fut,
  sw_gs_at_past=sw_gs_at_past, sw_gs_at_fut=sw_gs_at_fut, le_gs_at_past=le_gs_at_past, le_gs_at_fut=le_gs_at_fut,
  DIST_at_wet_past=DIST_at_wet_past, DIST_at_wet_fut=DIST_at_wet_fut)
trend_layers <- lapply(trend_layers, terra::resample, y=land_cover); list2env(trend_layers, envir=environment())

trend_driver <- read.csv(file.path(data_dir,"potential_trend_driver.csv"))
sites_long <- read.csv(file.path(data_dir,"long_term_ec_sites.csv")) %>% mutate(source=ifelse(source=="AmeriFlux_BASE","Eddyproc","Oneflux"))
if (!"category" %in% names(trend_driver))
  trend_driver <- trend_driver %>% left_join(sites_long[,c("site_ID","IGBP","Climate.class","category","NEEpartition","source","AI")], by="site_ID")
if (USE_CO2_RATE) {  # swap the CO2 level column for the rate, if the refactor table exists
  cr <- read.csv(file.path(data_dir,"potential_trend_driver_co2rate.csv"))
  if ("CO2_at" %in% names(cr)) {
    trend_driver <- trend_driver %>%
      select(-any_of("CO2_at")) %>%
      left_join(cr[,c("site_ID","year_mid","CO2_at")], by=c("site_ID","year_mid")) %>%
      mutate(CO2 = CO2_at)
  }
}

climate_vegs <- c("cold_ENF","cold_DBF","cold_nonforest","dry","wet")
frame_names <- c(cold_ENF="cold_evergreen_forest", cold_DBF="cold_deciduous_forest", cold_nonforest="cold_non_forest", dry="dry", wet="wet-warm")
zone_value  <- c(cold_ENF=2, cold_DBF=1, cold_nonforest=3, dry=4, wet=5)
# numeric predictor set per group (matches ad_per_variable); CO2 is level or rate per USE_CO2_RATE
num_preds <- list(
  cold_ENF       = c("CO2","TA_gs_at","MAP","MAT","AI","GPP_at_mean"),
  cold_DBF       = c("CO2","TA_ngs_at","MAP","MAT","AI","GPP_at_mean"),
  cold_nonforest = c("SWC_gs_at","MAP","MAT","AI","GPP_at_mean"),
  dry            = c("CO2","TA_gs_at","P_gs_at","SWC_gs_at","MAP","MAT","AI","GPP_at_mean"),
  wet            = c("CO2","TA_at","LE_at","DIST_at","MAP","MAT","AI","GPP_at_mean"))

z_apply <- function(m, mu, sdv) sweep(sweep(as.matrix(m),2,mu,"-"),2,sdv,"/")
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}

build_frames <- function(cv, data_group) {
  cells <- which(terra::values(zone_group) == zone_value[[cv]])
  base <- data.frame(cell_id=cells,
    CO2=cell_values(co2_past,cells), TA_gs_at=cell_values(ta_gs_at_past,cells), TA_ngs_at=cell_values(ta_ngs_at_past,cells),
    P_gs_at=cell_values(tp_gs_at_past,cells), P_ngs_at=cell_values(tp_ngs_at_past,cells), LE_gs_at=cell_values(le_gs_at_past,cells),
    SW_gs_at=cell_values(sw_gs_at_past,cells), DIST_at=cell_values(DIST_at_wet_past,cells), MAP=cell_values(map_baseline,cells),
    MAT=cell_values(mat_baseline,cells), NEE=cell_values(nee_baseline,cells), GPP=cell_values(gpp_baseline,cells),
    LE=cell_values(le_baseline,cells), AI=cell_values(aridity,cells), GPP_at_mean=cell_values(gpp_at_baseline,cells),
    climate_id=cell_values(koppen_climate,cells), IGBP_id=cell_values(land_cover,cells)) %>%
    left_join(df_koppen_climate, by="climate_id") %>% left_join(df_IGBP, by="IGBP_id") %>% flatten_frame()
  # impute unseen categorical levels to the modal class (as the export script does)
  nc <- setdiff(unique(base$Climate.class), unique(data_group$Climate.class))
  ni <- setdiff(unique(base$IGBP), unique(data_group$IGBP))
  if (length(nc)) base$Climate.class[base$Climate.class %in% nc] <- names(which.max(table(base$Climate.class)))
  if (length(ni)) base$IGBP[base$IGBP %in% ni] <- names(which.max(table(base$IGBP)))
  past <- base
  futu <- base %>% mutate(CO2=cell_values(co2_fut,cells), TA_gs_at=cell_values(ta_gs_at_fut,cells), TA_ngs_at=cell_values(ta_ngs_at_fut,cells),
    P_gs_at=cell_values(tp_gs_at_fut,cells), P_ngs_at=cell_values(tp_ngs_at_fut,cells), LE_gs_at=cell_values(le_gs_at_fut,cells),
    SW_gs_at=cell_values(sw_gs_at_fut,cells), DIST_at=cell_values(DIST_at_wet_fut,cells)) %>% flatten_frame()
  # predict SWC_gs_at for BOTH periods (cold_nonforest & dry use it), same sub-model as export
  if (cv %in% c("cold_nonforest","dry")) {
    swc_rf <- ranger(SWC_gs_at ~ TA_gs_at + SW_gs_at + P_gs_at + CO2 + MAP + MAT + IGBP + Climate.class + AI,
                     data=data_group, mtry=2, min.node.size=5)
    past$SWC_gs_at <- predict(swc_rf, data=past)$predictions
    futu$SWC_gs_at <- predict(swc_rf, data=futu)$predictions
  }
  keep_na <- if (cv!="wet") function(d) d %>% filter(if_all(-c(DIST_at), ~!is.na(.))) else function(d) d %>% filter(if_all(-c(TA_ngs_at,P_ngs_at), ~!is.na(.)))
  # align rows across periods on cell_id after NA filtering
  past <- keep_na(past); futu <- keep_na(futu)
  if (cv=="wet") { ren <- c(TA_at="TA_gs_at", SW_at="SW_gs_at", P_at="P_gs_at", LE_at="LE_gs_at")
    for (nn in names(ren)) { names(past)[names(past)==ren[[nn]]] <- nn; names(futu)[names(futu)==ren[[nn]]] <- nn } }
  common <- intersect(past$cell_id, futu$cell_id)
  list(past=past %>% filter(cell_id %in% common) %>% arrange(cell_id),
       futu=futu %>% filter(cell_id %in% common) %>% arrange(cell_id))
}

# ---------------------------------------------------------------------------
# per group: fit RF / LIN / GP, predict past+future, bracket the CHANGE
# ---------------------------------------------------------------------------
summ <- list()
for (cv in climate_vegs) {
  message("=== ", cv, " ===")
  data_group <- switch(cv,
    cold_ENF = trend_driver %>% filter(IGBP %in% c("ENF","EBF") & category=="cold"),
    cold_DBF = trend_driver %>% filter(IGBP %in% c("DBF","DNF","MF") & category=="cold"),
    cold_nonforest = trend_driver %>% filter(!IGBP %in% c("ENF","EBF","DBF","DNF","MF") & category=="cold"),
    dry = trend_driver %>% filter(category=="dry"),
    wet = trend_driver %>% filter(category=="wet"))

  fr <- build_frames(cv, data_group)
  preds <- num_preds[[cv]]
  # training matrix
  tk <- complete.cases(data_group[,c(preds,"NEE_at")])
  Xtr <- data_group[tk,preds,drop=FALSE]; ytr <- data_group$NEE_at[tk]
  mu <- sapply(Xtr,mean); sdv <- sapply(Xtr,sd); sdv[sdv==0] <- 1
  Ztr <- z_apply(Xtr,mu,sdv); dtr <- data.frame(NEE_at=ytr, Ztr)
  Zp <- z_apply(fr$past[,preds,drop=FALSE], mu,sdv); Zf <- z_apply(fr$futu[,preds,drop=FALSE], mu,sdv)

  rf  <- ranger(NEE_at ~ ., data=dtr, num.trees=500, mtry=2, min.node.size=5, seed=128)
  lin <- lm(NEE_at ~ ., data=dtr)
  rf_ch  <- predict(rf, data.frame(Zf))$predictions - predict(rf, data.frame(Zp))$predictions
  lin_ch <- as.numeric(predict(lin, data.frame(Zf)) - predict(lin, data.frame(Zp)))

  gp_ch <- gp_sd <- rep(NA_real_, nrow(Zf)); gp_ok <- FALSE; gp_covtype <- NA_character_
  if (have_dk) {
    for (ct in c("matern5_2","matern3_2","gauss","exp")) {   # fall back through covtypes if km fails
      gp <- try(DiceKriging::km(design=as.data.frame(Ztr), response=ytr, covtype=ct,
                                nugget.estim=TRUE, control=list(trace=FALSE)), silent=TRUE)
      if (!inherits(gp,"try-error")) {
        pf <- DiceKriging::predict(gp, as.data.frame(Zf), type="UK", checkNames=FALSE)
        pp <- DiceKriging::predict(gp, as.data.frame(Zp), type="UK", checkNames=FALSE)
        gp_ch <- pf$mean - pp$mean; gp_sd <- pf$sd; gp_ok <- TRUE; gp_covtype <- ct; break
      }
    }
    if (!gp_ok) message("  GP failed for ", cv, " (all covtypes); bracket uses RF+LIN, mask uses Mahalanobis only.")
  }

  # Applicability masks computed on the CURRENT predictor set (self-consistent with
  # USE_CO2_RATE) -- NOT read from the stale CO2-level AD flags. Three components:
  #   (a) Mahalanobis distance of each future cell to the training centroid vs chi-sq(0.975)
  #   (b) kNN mean distance <= the 95th percentile of training self-distances
  #   (c) GP predictive SD <= the 97.5th percentile of leave-one-out training SD.
  #       LOO uncertainty is the appropriate reference; fitted-point SD is artificially tiny.
  Smat <- tryCatch(if (requireNamespace("MASS", quietly=TRUE)) MASS::cov.rob(Ztr)$cov else stats::cov(Ztr),
                   error=function(e) stats::cov(Ztr))
  Sinv <- tryCatch(solve(Smat), error=function(e) MASS::ginv(Smat))
  d2   <- stats::mahalanobis(Zf, center=colMeans(Ztr), cov=Sinv, inverted=TRUE)
  in_dom_maha <- d2 <= stats::qchisq(0.975, df=ncol(Ztr))

  in_dom_knn <- rep(NA, nrow(Zf)); knn_threshold <- NA_real_
  if (requireNamespace("FNN", quietly=TRUE)) {
    k <- min(5L, nrow(Ztr) - 1L)
    self_dist <- FNN::get.knn(Ztr, k=k+1L)$nn.dist[, -1, drop=FALSE]
    knn_threshold <- unname(stats::quantile(rowMeans(self_dist), 0.95, na.rm=TRUE))
    future_dist <- FNN::get.knnx(Ztr, Zf, k=k)$nn.dist
    in_dom_knn <- rowMeans(future_dist) <= knn_threshold
  }

  gp_sd_threshold <- NA_real_; loo_sd_median <- NA_real_
  if (gp_ok) {
    loo <- try(DiceKriging::leaveOneOut.km(gp, type="UK"), silent=TRUE)
    if (!inherits(loo, "try-error") && any(is.finite(loo$sd))) {
      loo_sd_median <- stats::median(loo$sd, na.rm=TRUE)
      gp_sd_threshold <- unname(stats::quantile(loo$sd, 0.975, na.rm=TRUE))
    } else {
      message("  GP LOO calibration failed for ", cv, "; GP mask omitted.")
    }
  }
  in_dom_gp <- if (!is.finite(gp_sd_threshold)) rep(NA, nrow(Zf)) else gp_sd <= gp_sd_threshold
  in_dom <- in_dom_maha & (is.na(in_dom_knn) | in_dom_knn) & (is.na(in_dom_gp) | in_dom_gp)

  gp_finite_frac <- mean(is.finite(gp_ch))
  gp_change_sd <- stats::sd(gp_ch, na.rm=TRUE)
  if (gp_ok && (!is.finite(gp_change_sd) || gp_change_sd < sqrt(.Machine$double.eps)))
    message("  GP change is missing or effectively constant for ", cv, " (covtype=", gp_covtype, ").")

  agree3 <- (sign(rf_ch)==sign(lin_ch)) & (is.na(gp_ch) | sign(rf_ch)==sign(gp_ch))
  width  <- apply(cbind(rf_ch,lin_ch,gp_ch), 1, function(v){v<-v[is.finite(v)]; if(length(v)<2) NA else max(v)-min(v)})
  survives <- (is.na(in_dom) | in_dom) & agree3

  cells <- tibble(cell_id=fr$futu$cell_id, rf_change=rf_ch, lin_change=lin_ch, gp_change=gp_ch,
                  gp_sd=gp_sd, in_domain_maha=in_dom_maha, in_domain_knn=in_dom_knn, in_domain_gp=in_dom_gp,
                  in_domain=in_dom, sign_agree=agree3, bracket_width=width, survives=survives)
  write.csv(cells, file.path(outdir, paste0("bracket_cell_", frame_names[[cv]], ".csv")), row.names=FALSE)

  # maps
  for (v in c("gp_sd","bracket_width","survives","rf_change","lin_change","gp_change")) {
    r <- zone_group; terra::values(r) <- NA; terra::values(r)[cells$cell_id] <- as.numeric(cells[[v]])
    terra::writeRaster(r, file.path(outdir, sprintf("%s_%s.tif", v, frame_names[[cv]])), overwrite=TRUE)
  }

  summ[[cv]] <- tibble(group=frame_names[[cv]], n_cells=nrow(cells),
    cor_rf_lin=safe_cor(rf_ch,lin_ch),
    cor_rf_gp =safe_cor(rf_ch,gp_ch),
    sign_agree_rf_lin=mean(sign(rf_ch)==sign(lin_ch),na.rm=TRUE), sign_agree_all3=mean(agree3,na.rm=TRUE),
    median_bracket_width=median(width,na.rm=TRUE), median_abs_rf_change=median(abs(rf_ch),na.rm=TRUE),
    gp_ok=gp_ok, gp_covtype=gp_covtype, gp_finite_frac=gp_finite_frac, gp_change_sd=gp_change_sd,
    loo_sd_median=loo_sd_median, gp_sd_threshold=gp_sd_threshold, knn_threshold=knn_threshold,
    frac_in_domain_maha=mean(in_dom_maha,na.rm=TRUE),
    frac_in_domain_knn=if (all(is.na(in_dom_knn))) NA_real_ else mean(in_dom_knn,na.rm=TRUE),
    frac_in_domain_gp=if (all(is.na(in_dom_gp))) NA_real_ else mean(in_dom_gp,na.rm=TRUE),
    frac_in_domain=mean(in_dom,na.rm=TRUE), frac_survives=mean(survives,na.rm=TRUE))
}

out <- bind_rows(summ)
write.csv(out, file.path(outdir,"bracket_summary.csv"), row.names=FALSE)
cat("\n================ BRACKET SUMMARY ================\n"); print(as.data.frame(out), digits=3)
cat("\nfrac_survives = share of projected cells in-domain AND sign-agreed across RF/LIN/GP",
    "\n  -> the honest 'how much of Fig. 4 holds up'. median_bracket_width = the projection's structural uncertainty.\n")
if (!have_dk) cat("\n[note] install.packages('DiceKriging') to add the GP arm + SD mask (ran RF vs LIN only).\n")
if (!USE_CO2_RATE) cat("[note] USE_CO2_RATE=FALSE: this is the CURRENT (CO2-level) projection. Set TRUE after building the CO2-rate rasters to see the post-fix picture.\n")
