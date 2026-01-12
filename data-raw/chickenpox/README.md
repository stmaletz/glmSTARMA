## Description
This folder contains raw data and code to process the raw data into the final datasets included in the `glmSTARMA` package.

The files `hungary_chickenpox.csv` and `hungary_county_edges.csv` are the files downloaded from the UCI Machine Learning Repository  
(https://archive.ics.uci.edu/dataset/580/hungarian+chickenpox+cases).
The file `stadat-nep0034-22.1.2.1-hu.csv` was downloaded from the Hungarian Central Statistical Office (https://www.ksh.hu/stadat_files/nep/en/nep0034.html) and contains the annual population data for the regions of Hungary.

The R script `process_chickenpox_data.R` processes these raw data files and creates the final dataset `chickenpox`, which is included in the `glmSTARMA` package.


