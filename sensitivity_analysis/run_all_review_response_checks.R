# Run all review-response diagnostics in sequence.

scripts <- c(
  "01_cv_diagnostics.R",
  "02_co2_identification_diagnostics.R",
  "03_multiplicity_fdr_summary.R",
  "04_window_overlap_sensitivity.R",
  "05_rf_extrapolation_range_check.R",
  "06_swc_sensitivity_check.R",
  "07_reproducibility_audit.R",
  "08_classification_stability_check.R"
)

for (script in scripts) {
  path <- file.path("scripts2026", "review_response", script)
  cat("\n\n================", script, "================\n")
  source(path, local = new.env(parent = globalenv()))
}
