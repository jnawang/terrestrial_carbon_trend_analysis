# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:30:27 2025

@author: ammar
"""

import cdsapi
import datetime
import os

# Set the API credentials as environment variables
os.environ['CDSAPI_URL'] = 'https://cds.climate.copernicus.eu/api'
os.environ['CDSAPI_KEY'] = '******'  # get that key from your CDS account; Junna's account


dataset = "derived-era5-single-levels-daily-statistics"
request = {
    "product_type": "reanalysis",
    "variable": ["2m_temperature"],
    "year": "2012",
    "month": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12"
    ],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15",
        "16", "17", "18",
        "19", "20", "21",
        "22", "23", "24",
        "25", "26", "27",
        "28", "29", "30",
        "31"
    ],
    "daily_statistic": "daily_mean",
    "time_zone": "utc+00:00",
    "frequency": "1_hourly"
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()


##########################################################################################################
# surface_solar_radiation_downwards; total_precipitation; 2m_temperature
###R code to generate year series and coordinates for ERA5 API. 
# cat(paste0('"', 1998:2024, '"', collapse = ", "))
# 
# coor <- round(c(45.8059,	-90.0799), 2)
# 
# paste(as.character(round(c(coor[1] + 0.006, coor[2]-0.006, coor[1] - 0.006, coor[2]+0.006), 2)), collapse = ", ")

