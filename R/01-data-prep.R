library(data.table)
library(magrittr)
library(reticulate)
library(fst)
library(here)

# GMST --------------------------------------------------------------------

# * HadCRUT4 --------------------------------------------------------------

## Base period: 1961-1990
had =	fread(here('data/raw/had.csv'))
had = had[year <= year(Sys.time())-1, 
					.(year, had = med, had_025 = uncert_025, had_975 = uncert_975)]

# * Cowtan & Way (2014) ---------------------------------------------------

## Base period: 1961-1990
cw14 = fread(here('data/raw/cw14.csv'))
cw14 = cw14[year <= max(had$year), .(year, cw, cw_1sigma)]

# * GISTEMP ---------------------------------------------------------------

## Base period: 1951-1980 (Note difference compared to HadCRUT4 and CW2014)
giss = fread(here('data/raw/giss.csv'))

## Subset data
giss = giss[Year <= max(had$year), .(year = Year, giss = `J-D`)] 

## Adjust to match HadCRUT4 base period (i.e 1961-1990)
giss_rebase = giss[year >= 1961 & year <= 1990, mean(giss, na.rm=TRUE)]
giss[, giss := giss - round(giss_rebase, 2)]
rm(giss_rebase)



# * Merged GMST dataset ---------------------------------------------------

## Comparison plot
plot(had$year, had$had, type = "l")
lines(cw14$year, cw14$cw, col = "red")
lines(giss$year, giss$giss, col = "blue")
abline(v = 1961, lty = 4)
abline(v = 1990, lty = 4)
abline(h = 0, lty = 2)
abline(v = 1871, lty = 4, col = "green")
abline(v = 1900, lty = 4, col = "green")
dev.off()

## Merge
gmst = merge(merge(had, cw14, all.x = TRUE), giss, all.x = TRUE)

## Rebase all GMST products down by HadCRUT4 pre-industrial temperature
## Using 1871-1900, but makes no material difference if chose, say, 1851-1880.
had_rebase = gmst[year >= 1871 & year <= 1900, mean(had, na.rm=TRUE)]
gmst_cols = c('had', 'had_025', 'had_975', 'cw', 'giss')
gmst[, (gmst_cols) := lapply(.SD, function(x) x - had_rebase), .SDcols = gmst_cols]
rm(had_rebase, gmst_cols)

plot(gmst$year, gmst$had, type = "l")
lines(gmst$year, gmst$cw, col = "red")
lines(gmst$year, gmst$giss, col = "blue")
abline(v = 1961, lty = 4)
abline(v = 1990, lty = 4)
abline(v = 1871, lty = 4, col = "green")
abline(v = 1900, lty = 4, col = "green")
abline(h = 0, lty = 2, col = "green")
dev.off()

# Forcings ----------------------------------------------------------------

# * RCPs ------------------------------------------------------------------

rcps = fread(here('data/raw/rcps.csv'))

colnames(rcps) = gsub("_rf", "", tolower(colnames(rcps)))

## Will only be using some variables, so select the revelant columns plus a few
## extra. (Also, see the RCP flowchart in the ./Data folder to get a sense of how 
## the individual forcings add up to the total radiative forcings column.)
rcps = rcps[, .(year = years, 
								trf = total_inclvolcanic - volcanic_annual, 
								volc = volcanic_annual, 
								solar,
								anthro = total_anthro, 
								ghg, co2ch4n2o, co2, fgassum, mhalosum, 
								aerosols = totaer_dir,
								other = cloud_tot + stratoz + tropoz + ch4oxstrath2o + landuse + bcsnow,
								rcp)]


# Ocean-atmospheric phenomena ---------------------------------------------

# * SOI -------------------------------------------------------------------

soi =	fread(here('data/raw/soi.csv'))
soi[, soi := rowMeans(.SD), .SDcols = 2:13]
soi = soi[year <= max(gmst$year), .(year, soi)]

# * AMO -------------------------------------------------------------------

amo = fread(here('data/raw/amo.csv'))
amo[, amo := rowMeans(.SD), .SDcols = 2:13]
amo = amo[, .(year, amo)]


# Full merged dataset -----------------------------------------------------

climate =
	merge(gmst, rcps, all = TRUE) %>%
	merge(soi, all.x = TRUE) %>%
	merge(amo, all.x = TRUE)

front_cols = c('year', 'rcp', 'had', 'had_025', 'had_975', 'cw', 'cw_1sigma', 'giss')
other_cols = setdiff(colnames(climate), front_cols)
setcolorder(climate, c(front_cols, other_cols))
setorder(climate, rcp, year)

## Subset the data to a common period for analysis
comm_min_year = max(min(gmst$year), min(rcps$year), min(soi$year), min(amo$year)) 
climate = climate[year >= comm_min_year & year <= 2100]

## Lastly, adjust temperature series to match end of historic rcp data (i.e. 2005)
## and create "mean" versions of volc, soi, and amo variables for future projections.
climate =	
	climate[, ':=' (had_full = had,	cw_full = cw,	giss_full = giss)] %>%
	.[year > 2005, ':=' (had = NA, cw = NA, giss = NA)] %>%
	.[volc == 0, volc := NA] %>%
	.[, ':=' (volc_mean = fifelse(year <= 2006, volc, mean(volc, na.rm = TRUE)), ## Note: Not 2005
						soi_mean = fifelse(!is.na(soi), soi, mean(soi, na.rm = TRUE)),
						amo_mean = fifelse(!is.na(amo), amo, mean(amo, na.rm = TRUE)))]

## Write the merged dataset to disk
fwrite(climate, here("data/climate.csv"))


# Priors data frame -------------------------------------------------------

## Define priors for our group of sceptics. For the noninformative prior, I'll
## actually follow the Stan recommendations -- e.g. see rstanarm documentation --
## for "weakly informative" priors with mean 0 (for centered data) and a std
## deviation of 2.5 * sd(y) / sd(x)
ni_sigma = climate[rcp=='rcp26' & year <= 2005, 2.5 * (sd(had)/sd(trf))]
priors_df = 
	data.table(
		mu = c(0, 1, 1, 0, 0),
		sigma = c(ni_sigma, 0.25, 0.065, 0.25, 0.065),
		prior_type = c("ni", "luke", "luke", "den", "den"),
		convic_type = c("", "mod", "strong", "mod", "strong")
		)

fwrite(priors_df, here("data/priors.csv"))


# Dessler and Forster (2018) ----------------------------------------------

## We'll prep and save the 1,000 simulation ensemble from DF18 separately. 

## Need spio Python library (calling via reticulate)
spio = import('scipy.io')

df18_file = here('data/raw/df18.idlsave')
df18_idl = spio$readsav(df18_file, python_dict = TRUE, verbose = FALSE)

## The main columns that we're interested in, i.e. TRF and VOLC, are split
## across different parts of the df18_idl object. So we'll have to extract
## each one separately (over all of the 1,000 simulation runs) and then merge
## together at the end.

ar5 = df18_idl$ar5all; df18_idl$ar5all = NULL
ar5 = as.data.table(ar5)
setnames(ar5, old = c('V1', 'V2', 'V3', 'value'), new = c('forcing', 'year', 'sim', 'rf'))

ar5[, year := year+1750-1]

## From Zeke (originally from Piers). Note: RFari is aerosol-radiation interaction
# dict = c('CO2'=1, 'Other_WMGHG'=2, 'O3_trop'=3, 'O3_strat'=4, 'RFari'=5, 
#          'totalaerosolERF'=6, 'LUCERF'=7, 'StratWVapour'=8, 'BCsnow'=9, 
#          'Contrail_cirrus'=10, 'Solar'=11, 'Volcanic'=12)
dict = c('co2'=1, 'ghg_other'=2, 'o3_trop'=3, '03_strat'=4, 'rfari'=5, 'aer'=6, 
				 'luc'=7, 'h20'=8, 'bc'=9, 'contrails'=10, 'solar'=11, 'volc'=12)

ar5[, forcing := names(dict[forcing])]
ar5[, forcing := factor(forcing, levels = names(dict))]

rf = rbindlist(lapply(c('rf_total', 'rf_anthro', 'rf_nat'),
											function(i) {
												d = as.data.table(data.frame(df18_idl[[i]]))
												d[, `:=` (year = df18_idl$year, id = i)]
												return(d)
											}))
rf = melt(rf, id = c('id', 'year'), variable.name = 'sim', value.name = 'rf')
rf[, sim := as.integer(gsub("X", "", sim))]

## Merge into single Dessler and Forster (2018) DT
df18 = merge(rf[id=='rf_total', .(year, sim, rf_tot = rf)], 
						 ar5[forcing=='volc', .(year, sim, volc = rf)])
## Use "trf" to be consistent with main data set
df18[, trf := rf_tot - volc][, rf_tot := NULL] 
## Optional subsetting, rearranging and ordering
df18 = df18[year >= 1866 & year <= 2005, .(sim, year, trf, volc)]
setorder(df18, sim)

## Use .fst file format since this is a large object.
write_fst(df18, here('data/df18.fst'))

rm(df18, df18_file, df18_idl, ar5, rf)