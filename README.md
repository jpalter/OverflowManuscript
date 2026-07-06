# OverflowManuscript
Codes used to create figures and other analyses for the manuscript exploring the rate of Overflow Transports out of the Nordic Seas.

1) Price_and_Barringer_FBC_loopy.m runs the streamtube model, configured for the Faroe Bank Channel, and looping over all the temperature and salinities.  Used to make the colors in Figure 4.  It can be modified for Denmark Strait location.

2) Initial_model_loop.ipynb - first step to get a dictionary of the CMIP6 models to be analyzed.  It interfaces with the XMIP tools and pangeo toolbox.

3)  Extract_profile_clean.ipynb - Loop through the dictionary of CMIP6 model output to extract profiles at locations that are on the Nordic Seas side of the Greenland Scotland Ridge, likely to provide source waters for overflow parameterizations of the DSOW (at Denmark Strait) and ISOW (near Faroe Bank Channel)

4)  Plot_temperature_map_Nordic.ipynb - Creates maps of the temperature for CMIP6 models in the vertical model level closest to 800 m for the Nordic Seas and subpolar North Atlantic, used to make Figure S1. 
   
5)  Plot_TS_DS.ipynb - Computes the mean temperature and salinity at the Denmark Strait location for CMIP6 model potential densities greater than 1027.8 (or between 350 and 700 m, if there are no waters dense enough in the profile).  Adds the model results to the T-S plots (i.e. the dots and arrows in Figure 4a).

6)  Plot_TS_FSC.ipynb - Computes the mean temperature and salinity at the Faroe Bank location for CMIP6 model potential densities greater than 1027.8 (or between 350 and 800 m, if there are no waters dense enough in the profile).  Adds the model results to the T-S plots (i.e. the dots and arrows in Figure 4b).

