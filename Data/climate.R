rm(list = ls()) # Clear data

library(readr)
library(dplyr)
require(tidyr)
library(ggplot2)


#####################################
## GLOBAL MEAN SURFACE TEMPERATURE ##
#####################################

### HadCRUT4 ###
## Website: http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html
## Base period: 1961-1990
had <- 
  tbl_df(
    read.table("http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.4.0.0.annual_ns_avg.txt",
               col.names = 
                 c("year", "med", 
                   "bias_025", "bias_975",
                   "me_025", "me_975",
                   "cvrg_025", "cvrg_975",
                   "me_bias_025", "me_bias_975",
                   "uncert_025", "uncert_975")
               )
    )

## Rename some columns and subset the data as needed
had <- 
  had %>%
  rename(had = med, had_025 = uncert_025, had_975 = uncert_975) %>%
  filter(year <= 2015) %>% ## Last year with full annual data as of writing
  select(year, had, had_025, had_975)


### Cowtan & Way (2014) ###
## Website: http://www-users.york.ac.uk/~kdc3/papers/coverage2013/series.html
## Base period: 1961-1990
## Note: Use updated version 2 of their estimates

cw <- 
  tbl_df(
    read.table("http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_krig_annual_v2_0_0.txt",
               col.names = 
                 c("year", "cw", "cw_1sigma", "cw_cvrg_1sigma", "cw_ensmbl_1sigma")
               )
    )
  
cw <-
  cw %>%
  filter(year <= max(had$year)) %>%
  select(year, cw, cw_1sigma)


## GISTEMP ##
## Website: http://data.giss.nasa.gov/gistemp/
## Base period: 1951-1980 (Note difference compared to HadCRUT4 and CW2014)

giss <- 
  tbl_df(
    read.csv("http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.csv",
             na.strings = c("***", "****"), stringsAsFactors = F, skip = 1
             )
    )

## Subset data
giss <- 
  giss %>%
  rename(year = Year, giss = J.D) %>%
  # mutate(giss = 0.01 * giss) %>%
  filter(year <= max(had$year)) %>%
  select(year, giss)

## Adjust to match HadCRUT4 base period (i.e 1961-1990)
giss_rebase <- mean((giss %>% filter(year >= 1961 & year <= 1990))$giss)
giss <- 
  giss %>%
  mutate(giss = giss - round(giss_rebase, 2))
rm(giss_rebase)

plot(had$year, had$had, type = "l")
lines(cw$year, cw$cw, col = "red")
lines(giss$year, giss$giss, col = "blue")
abline(v = 1961, lty = 4)
abline(v = 1990, lty = 4)
abline(h = 0, lty = 2)
abline(v = 1871, lty = 4, col = "green")
abline(v = 1900, lty = 4, col = "green")


## Merged GMST dataset
gmst <- left_join(left_join(had, cw), giss)

## Rebase all GMST products down by HadCRUT4 pre-industrial temperature
## Using 1871-1900, but makes no material difference if chose, say, 1851-1880.
had_rebase  <- mean((gmst %>% filter(year >= 1871 & year <= 1900))$had)
gmst <-
  gmst %>%
  gather(key, value, -year) %>%
  mutate(value = ifelse(key != "cw_1sigma", value - had_rebase, value)) %>%
  spread(key, value)
rm(had_rebase)

plot(gmst$year, gmst$had, type = "l")
lines(gmst$year, gmst$cw, col = "red")
lines(gmst$year, gmst$giss, col = "blue")
abline(v = 1961, lty = 4)
abline(v = 1990, lty = 4)
abline(v = 1871, lty = 4, col = "green")
abline(v = 1900, lty = 4, col = "green")
abline(h = 0, lty = 2, col = "green")
# abline(v = 1851, lty = 4, col = "grey")
# abline(v = 1880, lty = 4, col = "grey")
# abline(h = mean((gmst %>% filter(year >= 1851 & year <= 1880))$had), lty = 2, col = "grey")


################################
### RADIATIVE FORCINGS (RCP) ###
################################

## Download and combine the four different RCP forcing scenarios

url1 <- "http://www.pik-potsdam.de/~mmalte/rcps/data/RCP"
url2 <- "_MIDYEAR_RADFORCING.DAT"

rcps <-
  bind_rows(
    mapply(function(x, y) {
      tbl_df(read.table(paste0(url1, x, url2), skip = 59, header = T)) %>%
        mutate(rcp = y)
      },
      c("3PD", "45", "6", "85"), ## "x" variable vector
      paste0("rcp", c(26, 45, 60, 85)), ## "y" variable vector
      SIMPLIFY = F
      )
    )

colnames(rcps) <- gsub("_rf", "", tolower(colnames(rcps)))

## Will only be using some variables, so select the revelant columns plus a few
## extra. (Also, see the RCP flowchart in the ./Data folder to get a sense of how 
## the individual forcings add up to the total radiative forcings column.)
rcps <-
  rcps %>%
  rename(year = years, 
         trf_inclvolc = total_inclvolcanic, 
         volc = volcanic_annual, 
         anthro = total_anthro,
         aerosols = totaer_dir) %>%
  mutate(trf = trf_inclvolc - volc,
         other = cloud_tot + stratoz + tropoz + ch4oxstrath2o + landuse + bcsnow) %>%
  select(year, trf, volc, solar, anthro, ghg, co2ch4n2o, 
         fgassum, mhalosum, aerosols, other, rcp) %>% 
  # select(year, trf, volc, solar, anthro, co2, rcp) %>%
  group_by(rcp)


#####################################
### OCEAN-ATMOSPHERIC PHENOMENA ###
#####################################

### SOI ###
## Website: http://www.cgd.ucar.edu/cas/catalog/climind/soi.html
## Note: using the SOI signal (Standardized Tahiti -- Standardized Darwin) series.
soi <- 
  tbl_df(
    read.table("http://www.cgd.ucar.edu/cas/catalog/climind/SOI.signal.ascii",
               na.strings = c("-99.9"),
               col.names = c("year", tolower(month.abb))
               )
    )

soi <- 
  soi %>%
  mutate(soi = rowMeans(.[2:13])) %>%
  filter(year <= max(gmst$year)) %>%
  select(year, soi)

### AMO ###
## Website: http://www.esrl.noaa.gov/psd/data/timeseries/AMO/
## Note: using unsmoothed monthly series. Will take annual averages anyway.

amo <- 
  read_table("http://www.esrl.noaa.gov/psd/data/correlation/amon.us.long.data",
             skip = 1, na = "-99.990", n_max = max(gmst$year) - 1856 + 1,
             col_names = c("year", tolower(month.abb))
             )

amo <- 
  amo %>%
  mutate(amo = rowMeans(.[2:13])) %>%
  select(year, amo)



######################
### MERGED DATASET ###
######################

climate <-
  gmst %>%
  full_join(rcps) %>%
  left_join(soi) %>%
  left_join(amo) %>% 
  select(year, rcp, had, had_025, had_975, cw, cw_1sigma, giss, everything()) %>%
  arrange(rcp, year)

## Subset the data to a common period for analysis
climate <-
  climate %>%
  filter(year >= 1866) %>%
  filter(year <= 2100)

## Lastly, adjust temperature series to match end of historic rcp data (i.e. 2005)
## and create "mean" versions of volc, soi, and amo variables for future projections.

climate <- 
  climate %>%
  mutate(had_full = had,
         cw_full = cw,
         giss_full = giss) %>%
  mutate(had = ifelse(year <= 2005, had, NA),
         cw = ifelse(year <= 2005, cw, NA),
         giss = ifelse(year <= 2005, giss, NA)) %>%
  mutate(volc_mean = ifelse(volc != 0, volc, NA)) %>%
  mutate(volc_mean = ifelse(year <= 2006, volc_mean, mean(volc_mean, na.rm = T)), ## Note: not 2005
         soi_mean = ifelse(!is.na(soi), soi, mean(soi, na.rm = T)),
         amo_mean = ifelse(!is.na(amo), amo, mean(amo, na.rm = T)))

## Finally, write to disk
write_csv(climate, "./Data/climate.csv")

