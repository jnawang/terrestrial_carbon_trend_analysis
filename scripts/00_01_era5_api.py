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
os.environ['CDSAPI_KEY'] = '47eed9ac-8876-4a98-93b1-d0a00c32e849' # get that key from your CDS account; Junna's account
# os.environ['CDSAPI_KEY'] = '8938d21f-2336-42c4-8fe5-5783996193b0'   # this is David's account

# download_dir = "/Users/jw2946/Downloads/"  # Change this to your desired directory

# Ensure the directory exists
# os.makedirs(download_dir, exist_ok=True)

# Define the output file path
# output_file = os.path.join(download_dir, "US-CMW_2011-2022.nc")  # Change the filename if needed

# Initialize CDS API client
client = cdsapi.Client()

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "Total precipitation"
    ],
    "year": [
        "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011"
    ],
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
    "time": [
        "00:00", "01:00", "02:00",
        "03:00", "04:00", "05:00",
        "06:00", "07:00", "08:00",
        "09:00", "10:00", "11:00",
        "12:00", "13:00", "14:00",
        "15:00", "16:00", "17:00",
        "18:00", "19:00", "20:00",
        "21:00", "22:00", "23:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [-2.85, -54.96, -2.86, -54.95]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download(target=output_file)

##########################################################################################################

###R code to generate year series and coordinates for ERA5 API. 
# cat(paste0('"', 1998:2024, '"', collapse = ", "))
# 
# coor <- round(c(45.8059,	-90.0799), 2)
# 
# paste(as.character(round(c(coor[1] + 0.006, coor[2]-0.006, coor[1] - 0.006, coor[2]+0.006), 2)), collapse = ", ")

