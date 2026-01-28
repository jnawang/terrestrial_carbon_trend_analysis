# A function to compare two SWC products (measured vs modelled)
# 1/28/2026, Junna Wang

# data I needed: site, yearStart, yearEnd, trend_obsv, trend_anls
SWC_trend <- data.frame(site = character(), yearStart = integer(), yearEnd = integer(), 
                        trend_obsv = double(), trend_anls = double())

# function to calculate sens' slope
cal_sens_slope <- function(x) {
  x <- x[!is.na(x)]
  # x <- scale(x)
  return(as.numeric(sens.slope(x)[c(1,3)]))
}

cal_sens_slope_scale <- function(x) {
  x <- x[!is.na(x)]
  x <- scale(x)
  return(as.numeric(sens.slope(x)[c(1,3)]))
}


# Find all the sites with FLUXNET products with good measured SWC values



# compare these sites with climate reanalysis products


# the last thing is to ask models to obtain more results. 




