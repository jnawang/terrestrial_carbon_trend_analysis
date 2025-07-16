# Long-term trends and drivers of net carbon uptake of global terrestrial ecosystems

This project uses global eddy covariance CO2 flux data from AmeriFlux, FLUXNET, and ICOS to examine the long-term CO2 flux trends and their drivers of the terrestrial biosphere.

## Objectives

-   Calculate decadal trends of CO2 fluxes over the past two decades across a wide range of terrestrial ecosystems.
-   Identify the drivers of CO2 trends for different ecosystems.
-   Estimate the changes in CO2 flux trends at two-decade-long sites.

## Structure of R and Python scripts

-   00_download_data
    -   00_download_flux_data
    -   00_01_era5_api_daily.py (API for downloading daily global climate data from ERA5)
    -   00_01_era5_api.py (API for downloading hourly climate data at certain locations from ERA5)
    -   00_02_ERA5_csv.rmd (convert data from ERA5 to csv files)
-   01_data_harmonization
    -   01_01_Ustar-filter_gap-fill_partition_AmeriFlux-BASE_sites.Rmd
    -   01_02_connect_European_FLUXNET_recent_ICOS.Rmd
    -   01_03_unify_four_flux_data.Rmd
-   02_quantify_NEE_trends_significant_trend_drivers
    -   02_01_growing_season_detection.Rmd
    -   02_02_trend_abrupt_change_sen.Rmd (calculate NEE trends by sen slope and their potential changes)
    -   02_03_trend_drivers_sen.Rmd (calculate potential driver variables)
    -   02_04_trend_drivers_causal_inference_brm.Rmd
-   03_scale_up
    -   03_00_download_future_climate.py
    -   03_01_upscale_global_NEE_trends.Rmd
