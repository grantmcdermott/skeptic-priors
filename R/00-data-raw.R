library(data.table)
library(here)

# GMST --------------------------------------------------------------------

# * HadCRUT4 --------------------------------------------------------------

## Website: http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html
## Base period: 1961-1990
had =	
	fread(
		"https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.6.0.0.annual_ns_avg.txt",
		col.names = c("year","med","bias_025","bias_975","me_025","me_975",
									"cvrg_025","cvrg_975", "me_bias_025","me_bias_975",
									"uncert_025", "uncert_975")
				)
fwrite(had, here('data/raw/had.csv'))

# * Cowtan & Way (2014) ---------------------------------------------------

## Website: http://www-users.york.ac.uk/~kdc3/papers/coverage2013/series.html
## Base period: 1961-1990
## Note: Use updated version 2 of their estimates
cw14 =
	fread(
		"http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_krig_annual_v2_0_0.txt",
		col.names = c("year", "cw", "cw_1sigma", "cw_cvrg_1sigma", "cw_ensmbl_1sigma")
		)
fwrite(cw14, here('data/raw/cw14.csv'))

# * GISTEMP ---------------------------------------------------------------

## Website: http://data.giss.nasa.gov/gistemp/
## Base period: 1951-1980 (Note difference compared to HadCRUT4 and CW2014)
giss =
	fread(
		"https://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.csv",
		na.strings = c("***", "****"), skip = 1
		)
fwrite(giss, here('data/raw/giss.csv'))

# Forcings ----------------------------------------------------------------

# * RCPs ------------------------------------------------------------------

## Download and combine the four different RCP forcing scenarios
url1 = "http://www.pik-potsdam.de/~mmalte/rcps/data/RCP"
url2 = "_MIDYEAR_RADFORCING.DAT"
rcps =
	rbindlist(
		mapply(function(x, y) {
			dat = fread(paste0(url1, x, url2), skip = 59, header = TRUE)
			dat[, rcp := y]
			return(dat)
		},
		c("3PD", "45", "6", "85"), ## "x" variable vector
		paste0("rcp", c(26, 45, 60, 85)), ## "y" variable vector
		SIMPLIFY = F
		)
		)
fwrite(rcps, here('data/raw/rcps.csv'))


# * Dessler and Forster (2018) --------------------------------------------

## Forcing data via Hausfather et al. (2019)
df18_url = 'https://github.com/hausfath/OldModels/raw/master/data/raw/forcing_data/1750-Oct2017_forcings.idlsave'
download.file(df18_url, here('data/raw/df18.idlsave'))


# Ocean-atmospheric phenomena ---------------------------------------------

# * SOI -------------------------------------------------------------------

## Website: http://www.cgd.ucar.edu/cas/catalog/climind/soi.html
## Note: using the SOI signal (Standardized Tahiti -- Standardized Darwin) series.
soi =
	fread(
		"http://www.cgd.ucar.edu/cas/catalog/climind/SOI.signal.ascii",
		na.strings = c("-99.9"), col.names = c("year", tolower(month.abb))
		)
fwrite(soi, here('data/raw/soi.csv'))

# * AMO -------------------------------------------------------------------

## Website: http://www.esrl.noaa.gov/psd/data/timeseries/AMO/
## Note: using unsmoothed monthly series. Will take annual averages anyway.
amo =
	fread(
		"http://www.esrl.noaa.gov/psd/data/correlation/amon.us.long.data",
		skip = 1, na = "-99.990", nrows = max(had$year) - 1856 + 1,
		col.names = c("year", tolower(month.abb))
		)
fwrite(amo, here('data/raw/amo.csv'))
