# Long-term trends and drivers of carbon dynamics of the terrestrial biosphere
This project uses global eddy covariance CO2 flux data from AmeriFlux, FLUXNET, and ICOS to examine the long-term CO2 flux dynamics of the terrestrial biosphere.  

## Objectives
-   Calculate decadal trends of CO2 fluxes over the past three decades across a wide range of terrestrial ecosystems.
-   Identify the drivers of CO2 trends for different ecosystems.
-   Estimate the frequency of abrupt change in CO2 fluxes across different periods

## Structure of Rscript
-   00_download_flux_data
-   01_data_harmanization
    -   01_01_Ustar-filtering_gap-filling_AmeriFlux-BASE_sites
    -   01_02_process_qualified_FLUXNET_ICOS_sites 
-   02_Quantify_abrupt_change_in_CO2_fluxes_time-series
-   03_Calculate_trend_of_CO2_fluxes_and_potential_drivers
-   04_Identify_drivers_of_CO2_fluxes_across_ecosystem_types
