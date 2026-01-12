## Description
This folder contains raw data and code to process the raw data into the final datasets included in the `glmSTARMA` package.

The subdirectory `raw` contains the original files downloaded from the RKI SurvStat@RKI application.

The file `germany_county_shapes.json` contains the shapefiles for the German counties, which were copied from https://github.com/ostojanovic/BSTIM.

The Census 2022 data for Berlin's districts were downloaded from https://ergebnisse.zensus2022.de. They are provided in the file `census_berlin.csv`.

The annual population data for the German counties were downloaded from the Federal Statistical Office of Germany (Destatis) at https://www.destatis.de. They are provided in the file `12411-0015_de.csv`.

Data processing scripts are provided in the files `rota.R`. This script requires the R packages `data.table`, `stringdist`, `spdep`, `dplyr`, and `sf` for data manipulation and processing.