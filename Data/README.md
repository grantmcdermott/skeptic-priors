# Data

This folder contains the main `climate.csv` dataset used in the paper, as well as the *R* code that was used to construct the dataset in the first place (`R/climate.R`). The file `rcp_flowchart.pdf` summarizes the RCP forcings data (see further below). 

## Data sources

The sources for constructing the `climate.csv` dataset are described in the paper. Briefly, they are as follows:

| Variable | Product | Description | Period | Source |
|----------|---------|-------------|--------|--------|
| GMST     | HadCRUT4| Global mean surface temperature. Primary series. Compiled by the UK Met Office and the Climatic Research Unit at the University of East Anglia. | 1850-2015 | [Link](http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html) |
|          | CW2014  | Secondary series. Compiled by Cowtan and Way (2014a,b). Corrects for coverage bias in HadCRUT4. | 1850-2015 | [Link](http://www-users.york.ac.uk/~kdc3/papers/coverage2013/series.html) |
|          | GISTEMP | Secondary series. Compiled by the NASA Goddard Institute for Space Studies. | 1880–2015 | [Link](http://data.giss.nasa.gov/gistemp/) |
| RF       | RCP     | Total radiative forcing due to anthropogenic and natural factors (excluding volcanic aerosols). Compiled by Meinshausen et al. (2011). Historical data until 2005, simulated scenarios thereafter. | 1765–2300 | [Link](http://www.pik-potsdam.de/~mmalte/rcps/) |
| AER      | RCP     | Radiative forcing due to volcanic stratospheric aerosols. Compiled by Meinshausen et al. (2011). | 1750-2005 | [Link](http://www.pik-potsdam.de/~mmalte/rcps/) |
| AMO      | NOAA    | Atlantic Multidecadal Oscillation. | 1856-2015 | [Link](http://www.esrl.noaa.gov/psd/data/timeseries/AMO/) |
| SOI      | NCAR    | Southern Oscillation Index. | 1866-2015 | [Link](http://www.cgd.ucar.edu/cas/catalog/climind/soi.html) |

### RCP forcings metadata

The file `rcp_flowchart.pdf` details how individual RCP forcings (solar irradience, GHGs, etc.) combine to form the total radiative forcing series that is used in the paper. This summarises the metadata that can be found within each RCP data file (e.g. [here](http://www.pik-potsdam.de/~mmalte/rcps/data/RCP3PD_MIDYEAR_RADFORCING.DAT)), partly reproduced below for convenience:

| Column | Variable | Description |
|--------|----------|-------------|
| 1 | TOTAL_INCLVOLCANIC_RF | Total anthropogenic and natural radiative forcing |
| 2 | VOLCANIC_ANNUAL_RF | Annual mean volcanic stratospheric aerosol forcing |
| 3 | SOLAR_RF | Solar irradience forcing |
| 4 | TOTAL_ANTHRO_RF | Total anthropogenic forcing |
| 5 | GHG_RF | Total greenhouse gas forcing (CO2, CH4, N2O, HFCs, PFCs, SF6, and Montreal Protocol gases) |
| 6 | KYOTOGHG_RF | Total forcing from greenhouse gases controlled under the Kyoto Protocol (CO2, CH4, N2O, HFCs, PFCs, SF6) |
| 7 | CO2CH4N2O_RF | Total forcing from CO2, methan and nitrous oxide. |
| 8 | CO2_RF | CO2 Forcing |
| 9 | CH4_RF | Methane Forcing |
| 10 | N2O_RF | Nitrous Oxide Forcing |
| 11 | FGASSUM_RF | Total forcing from all flourinated gases controlled under the Kyoto Protocol (HFCs, PFCs, SF6; i.e. columns 13-24) |
| 12 | MHALOSUM_RF | Total forcing from all gases controlled under the Montreal Protocol (columns 25-40) |
| 13-24 |  | Flourinated gases controlled under the Kyoto Protocol |
| 25-40 |  | Ozone Depleting Substances controlled under the Montreal Protocol |
| 41 | TOTAER_DIR_RF | Total direct aerosol forcing (aggregating columns 42 to 47) |
| 42 | OCI_RF | Direct fossil fuel aerosol (organic carbon) |
| 43 | BCI_RF | Direct fossil fuel aerosol (black carbon) |
| 45 | NOXI_RF | Direct nitrate aerosol |
| 47 | MINERALDUST_RF | Direct Forcing from mineral dust aerosol |
| 48 | CLOUD_TOT_RF | Cloud albedo effect |
| 49 | STRATOZ_RF | Stratospheric ozone forcing |
| 50 | TROPOZ_RF | Tropospheric ozone forcing |
| 51 | CH4OXSTRATH2O_RF | Stratospheric water-vapour from methane oxidisation |
| 52 | LANDUSE_RF | Landuse albedo |
| 53 | BCSNOW_RF | Black carbon on snow |
