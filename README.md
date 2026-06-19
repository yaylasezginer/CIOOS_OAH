CIOOS OAH README


Code collection for model evaluations provided within CIOOS Pacific Ocean Chemistry monitoring applicaiton: https://oa.cioospacificlabs.ca/ 
Model evaluations compare regional ocean biogeochemistry models, LiveOcean (https://faculty.washington.edu/pmacc/LO/LiveOcean.html) and SalishSeaCast (https://salishsea.eos.ubc.ca/nemo/) against available in-situ sensor data. 

__________________________________________________________________________________________________________________

MODEL EVALUATION WORKFLOW
1) Potential in-situ sensors are identified from CIOOS Data Explorer (https://explore.cioos.ca/?lang=en) and are downloaded from either CIOOS Pacific ERDDAP or directly from provider data source 
2) SensorNetork.csv is manually updated with new sensor information:
     SiteID,
     Latitude (degrees),
     Longitude (degrees),
     Depth (m or 'profile'),
     params (desired parameters using controlled vocabulary),
     params_sensor (parameter names given by sensor data provider),
     units (units given by sensor data provider).
   All proceeding functions read from this SensorNetwork.csv. Other fields in SensorNetwork.csv are updated autonomously
3) 'virtual moorings' are extracted from LiveOcean and SalishSeaCast models
        A) LiveOcean: (files referenced avialable in CIOOS_OAH>code>LiveOcean)
             Fill in extraction request form with virual mooring geographic coords and time range to be extracted. Pass extraction form to LiveOcean contacts
             Once available, manually update SensorNetwork.csv with local directory for LiveOcean virtual mooring
        B) SalishSeaCast: (files referenced available in: CIOOS_OAH>code>SalishSeaCast)
             SSC_extraction_shell.m reads SensorNetork.csv. For partially filled entries, the shell function calls: 
               i) use latlon2salishgrid.m to convert virtual mooring geographic coords to SalishSeaCast x and y grid points. Updates SesnsorNetwork.csv with x and y coordinates
               ii) Calls parallel_erddap.sh to automate virtual mooring extraction from SalishSeaCast ERDDAP. Hourly data is downloaded in monthly chunks to avoind overloading server
               iii) Call compile_SSC_sensor_network.m to aggregate monthly chunks into a single .nc virtual mooring file. Update SensorNetork with local directory for SalishSeaCast virtual mooring
4) ExportFilesCreate.m (CIOOS_OAH>ModelEvaluation) is used to:
             a) Read in-situ and virtual mooring data files
             b) Rewrite all parameters with controlled vocabulary
             c) Convert all data to standardized units
             d) Interpolate sensor data to have matching time stamps with LiveOcean and SalishSeaCast data - avoids interpolating across large data gaps associated with instrument down time
             e) Generate and save comparison figures
             f) Save model performance statistics (RMSE, residuals, R2)
             g) Save a .csv file that compiles unmodified LiveOcean, SalishSeaCast, and Sensor data in a single data matrix for user download
             h) Save a .nc file with the interpolated sensor data for internal use.  


__________________________________________________________________________________________________________________

DIRECTORY STRUCTURE

README.md - readme
config - contains SensorNetork.csv. This file tracks metadata associated with each site and tracks data processing 
code - 
	carbonate_system- Matlab code associated with carbonate system calculations and unit conversions
	SalishSeaCast- Matlab + linux code required to download and format SalishSeaCast virtual mooring data
	ModelEvaluation - comparison of LiveOcean, SalishSeaCast, and in-situ data. Produces final outputs seen on https://salishsea.eos.ubc.ca/nemo/


__________________________________________________________________________________________________________________

Disclaimer: code written in MatlabR2025b. Some functions may not work for older versions. Some functions may need to be downloaded if not already (e.g. setup_nctoolbox)

Last updated: June 2026
contact: yayla.sezginer@cioospacific.ca 

For further details on individual functions, call 'help [function_name_here]' to see internal details. 