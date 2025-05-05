##' .. this function is designed to calculate UTC offset given a local time . ..
##'
##' @title utc_offset, this function will produce a hour offset number converting local to UTC.
##' @param input lat, lon
##' @return UTC offset
##' @note Offset is the difference between UTC and local time without considering Daylight Saving Time
##' 
##' @title date_conversion, this function will produce a converted time, either utc -> local, or local -> utc
##' @param input lat, lon, time (string format), flag (flag == 0 local to utc, flag == 1 utc to local)
##' @note IT TAKES INTO ACCOUNT Daylight Saving Time!!
##' 
##' @author Boya ("Paul") Zhang
##' 

library(lutz)
library(lubridate)
#library(tmaptools)
library(sf)
#library(tzlookup)
#library(httr)
library(data.table)
#library(timechange)
#library(purrr)

#rm(list = ls())
dt <- data.table(  
  local_time = c("2018-01-16 22:02:37"),  
  utc_time =c("2018-01-16 22:02:37"),
  lat = c(25.0433),  
  long = c(-80.1918)  
)

lat <- 25.0433
lon <- -80.1918
time <- "2018-01-16 22:02:37"
flag <- 0 ### local to utc

check_DST <- function(lat, lon, timepoint) {
  # Get the timezone based on lat, lon
  timezone <- tz_lookup_coords(lat = lat, lon = lon, method = "accurate")
  
  # Convert the provided timepoint to the correct time zone
  timepoint <- as.POSIXct(timepoint, tz = timezone)
  
  # Check if the given time point is in Daylight Saving Time
  is_DST <- dst(timepoint)  # TRUE if in DST, FALSE if not
  
  # If it is in DST, adjust the time to the standard (non-DST) time
  if (is_DST) {
    # Subtract the DST offset to get standard time (non-DST)
    standard_time <- timepoint - hours(1)
  } else {
    # If not in DST, no adjustment is needed
    standard_time <- timepoint
  }
  
  # Return a list with DST status and the standard time if DST is in effect
  return(list(is_DST = is_DST, standard_time = standard_time))
}

timepoint <- now()  # This will get the current time

# Check if the current time is DST
is_DST <- check_DST(lat, lon, timepoint)

utc_offset <- function(lat, lon) {
  # Get timezone for the given coordinates
  timezone <- tz_lookup_coords(lat = lat, lon = lon, method = "accurate")
  
  if (is.na(timezone)) {
    stop("Timezone not found for given coordinates.")
  }
  
  # Get the current time in UTC
  time_utc <- now(tz = "UTC")  # No DST
  #time_utc <- as.POSIXct("2025-06-21 12:00:00", tz = "UTC") ## WITH DST
  #print(time_utc)
  
  # Convert UTC time to local time with proper timezone
  local_time_result <- with_tz(time_utc, tzone = timezone)  # Keeps time zone
  #print(local_time_result)
  
  # Check if the current time is in Daylight Saving Time (DST)
  result <- check_DST(lat, lon, local_time_result)
  #print(paste("Is it DST? ", result$is_DST))
  #print(paste("Standard time (if not in DST):", result$standard_time))
  local_revised <- result$standard_time
  
  conversion <- force_tz(local_revised, "UTC")
  #print(conversion)
  # Calculate the difference in hours, respecting time zones
  offset_hours <- as.numeric(difftime(conversion, time_utc, units = "hours"))
  
  print(offset_hours)
  return(offset_hours)
}



offset <- utc_offset(lat, lon)


date_conversion <- function(lat, lon, time, flag) {
 
  # Get timezone for the given coordinates
  timezone <- tz_lookup_coords(lat = lat, lon = lon, method = "accurate")
  print(timezone)
  if (flag == 1) {
    # Convert UTC to Local Time
    time <- as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    local_time_result <- with_tz(time, timezone)
    return(local_time_result)
  } else if (flag == 0) {
    # Convert Local Time to UTC
    time <- as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S")
    utc_time_result <- with_tz(time, "UTC")
    return(utc_time_result)
  } else {
    stop("Invalid flag. Please use 'utc_to_local' or 'local_to_utc'.")
  }
}

utc_result_winter <- date_conversion(lat, lon, "2018-01-16 22:02:37", 0)
utc_result_summer <- date_conversion(lat, lon, "2018-07-16 22:02:37", 0)

local_result <- date_conversion(lat, lon, "2018-01-16 22:02:37", 1)

####################### Block below used for dataframe with multiple rows ###########
############### converting utc to local ########
dt[, utc_time := as.POSIXct(utc_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]

# Get timezone for each coordinate
dt[, timezone := tz_lookup_coords(lat = lat, lon = long, method = "accurate")]

# Convert UTC time to local time using the determined timezone
dt[, local_time_result := with_tz(utc_time, timezone)]

####### converting local to utc #####
dt[, timezone := tz_lookup_coords(lat = lat, lon = long, method = "accurate")]  

dt[, local_time := as.POSIXct(local_time, format = "%Y-%m-%d %H:%M:%S", tz = timezone)]

dt[, utc_time_result := with_tz(local_time, "UTC")]

