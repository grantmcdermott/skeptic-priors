# Data

This folder contains the main dataset used in the paper — i.e. `climate.csv` —
as well as the *R* code that was used to construct the dataset in the first place. In addition, the file `rcp_flowchart.odp` provides a schematic of the RCP forcings data, which details how individual forcings (solar, GHGs, etc.) combine to form the total radiative forcing (RF) series that is used in the paper. Lastly, the PAGE09 folder contains estimates of the social cost of carbon (`scc.csv`), which result from running the PAGE09 integrated assement model on various posterior estimates of climate sensitivity. Since PAGE09 is proprietary, I am only able to provide the final output of the model.

## Data sources

The sources for constructing the `climate.csv` dataset are described in the paper. Briefly, they are as follows:

| Variable | Product | Description | Period | URL |
|----------|---------|-------------|--------|------|
| GMST     | HadCRUT4| Global mean surface temperature. Primary series. Compiled by the UK Met Office and the Climatic Research Unit at the University of East Anglia. | 1850-2015 | [Link](http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html) |
|          | CW2014  | Secondary series. Compiled by Cowtan and Way (2014a,b). Corrects for coverage bias in HadCRUT4. | 1850-2015 | [Link](http://www-users.york.ac.uk/~kdc3/papers/coverage2013/series.html) |
|          | GISTEMP | Secondary series. Compiled by the NASA Goddard Institute for Space Studies. | 1880–2015 | [Link](http://data.giss.nasa.gov/gistemp/) |
| RF       | RCP     | Total radiative forcing due to anthropogenic and natural factors (excluding volcanic aerosols). Compiled by Meinshausen et al. (2011). Historical data until 2005, simulated scenarios thereafter. | 1765–2300 | [Link](http://www.pik-potsdam.de/~mmalte/rcps/) |
| AER      | RCP     | Radiative forcing due to volcanic stratospheric aerosols. Compiled by Meinshausen et al. (2011). | 1750-2005 | [Link](http://www.pik-potsdam.de/~mmalte/rcps/) |
| AMO      | NOAA    | Atlantic Multidecadal Oscillation. | 1856-2015 | [Link](http://www.esrl.noaa.gov/psd/data/timeseries/AMO/) |
| SOI      | NCAR    | Southern Oscillation Index. | 1866-2015 | [Link](http://www.cgd.ucar.edu/cas/catalog/climind/soi.html) |
