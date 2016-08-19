rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
## Total chain length will thus be chain_length * n_chains
chain_length <- 10000
n_chains <- detectCores() - 1 

rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855)

##### 
## Marvel et al. (2016), hereafter MEA2016, provide the following scaled 
## efficacies (see table S1 in SI):
## http://data.giss.nasa.gov/modelforce/Marvel_etal2015.html
# Aerosols: 1.55* (1.05, 2.05)
# GHG:      1.17* (1.06, 1.27)
# Land use: 4.27 (-2.42, 10.95)
# Ozone:    0.66* (0.34, 0.98)
# Solar:    1.68 (-1.27, 4.63)
# Volc:     0.61* (0.33, 0.89)
## historical: 1.00 (0.83, 1.16)

## Note: efficacies are relative to CO2 (= 1)
# "The uncertainty in the efficacies is estimated from individual members of the
# single-forcing ensembles (Figure S2). Confidence intervals on the sample mean are
# constructed using a student-t distribution with 4 degrees of freedom (5 in the case
# of the 6-member historical ensemble)."

################################
# The main climate.csv data file only has aggregate forcings (e.g. landuse is 
# contained within the "other" series). So we need to read in the RCP data again
# and adjust the individual forcings to their effect values as per MEA16, before 
# joining these effective forcings to the main climate dataset. 

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

## Subset the data and change some column names for convenience.
rcps <-
  rcps %>%
  rename(year = years, 
         trf_inclvolc = total_inclvolcanic, 
         volc = volcanic_annual, 
         anthro = total_anthro,
         aerosols = totaer_dir) %>%
  mutate(trf = trf_inclvolc - volc) 

## Adjust the relevant columns to match up with the MEA2016 forcings. Note that 
# the forcings in Marvel at el. don't match up exactly to those in the RCP
# dataset. For example, they have a forcing called "ozone", which I assume 
# corresponds to the RCP forcing "mhalosum" (i.e. gases controlled under the 
# Montreal Protocol). See the RCP metadata and flowchart in the ./Data folder 
# for more information. Summary version: We construct a new effective TRF variable 
# using the following variables:
# TRF_eff (excl. volc) = 
#   1.68*solar +
#   1.17*co2ch4n2o + ## "well-mixed GHGs" in MEA2016
#   fgassum + ## flourinated gases, not adjusted
#   0.66*mhaslosum + ## "ozone" in MEA2016
#   1.55*aerosols + ## "direct aerosols" in MEA16
#   other (incl. 4.27*landuse)


####################################################
## 1) Adjust by effectiveness coefficient means only

rcps_eff <-
  rcps %>%
  filter(year >= 1866 & year <= 2100) %>%
  select(year, rcp, trf, volc, solar, co2ch4n2o, fgassum, mhalosum, aerosols,
         cloud_tot, stratoz, tropoz, ch4oxstrath2o, landuse, bcsnow ## "other" category
         ) %>% 
  ## Create effective forcing versions of the necessary variables
  mutate(volc_eff = 0.61*volc,
         solar_eff = 1.68*solar, 
         co2ch4n2o_eff = 1.17*co2ch4n2o, 
         # fgassum_eff = fgassum, 
         mhalosum_eff = 0.66*mhalosum, 
         aerosols_eff = 1.55*aerosols,
         landuse_eff = 4.27*landuse
         ) %>%
  mutate(other = landuse + cloud_tot + stratoz + tropoz + ch4oxstrath2o + bcsnow,
         other_eff = landuse_eff + cloud_tot + stratoz + tropoz + ch4oxstrath2o + bcsnow) %>%
  ## Sum these and (other variables) to get effective TRF
  mutate(trf_eff = solar_eff + co2ch4n2o_eff + fgassum + mhalosum_eff + aerosols_eff + other_eff) %>%
  ## Stricter version of effective TRF that doesn't adjust nonsig. coefs in MEA2016 (solar & landuse)
  mutate(trf_eff_sig = solar + co2ch4n2o_eff + fgassum + mhalosum_eff + aerosols_eff + other)

##
# Join new effective TRF (and volc) variables with main dataset
climate_eff <- 
  bind_cols(
    climate, 
    rcps_eff %>% select(trf_eff, trf_eff_sig, volc_eff)
    )
## Compare
ggplot(climate_eff, aes(x = year, y = trf, col = rcp)) +
  geom_line(aes(lty = "tcr")) +
  geom_line(aes(y = trf_eff, lty = "tcr_eff")) +
  geom_line(aes(y = trf_eff_sig, lty = "tcr_eff_sig")) +
  scale_colour_manual(values = c(rcp_cols)) 

## Run the Bayesian regression on the historical data using a noninformative prior
tcr_eff1 <- 
  lapply(
    c("trf_eff", "trf_eff_sig"),
    function(x){
      
      clim_df <- 
        climate_eff %>%
        filter(rcp == "rcp26") %>%
        filter(year <= 2005) %>%
        select_("had", x, "volc_eff", "soi", "amo") %>%
        rename_("trf_eff" = x)
      
      theta_sample_eff <- blinreg(clim_df$had, 
                                  cbind(alpha = 1, 
                                        beta = clim_df$trf_eff,
                                        gamma = clim_df$volc_eff, 
                                        delta = clim_df$soi, 
                                        eta = clim_df$amo), 
                                  chain_length*n_chains)
      
      tcr <- as.data.frame(theta_sample_eff[[1]])$Xbeta * rf2x
      
      tcr_df <- data_frame(tcr, series = x)
      
      rm(theta_sample_eff)
      
      return(tcr_df)
    }
  ) %>%
  bind_rows()

tcr_eff1 %>%
  group_by(series) %>%
  summarise(mean = round(mean(tcr), 1),
            q025 = round(quantile(tcr, .025), 1),
            q975 = round(quantile(tcr, .975), 1)) 
#      series  mean  q025  q975
#       <chr> <dbl> <dbl> <dbl>
#     trf_eff   2.1   1.9   2.4
# trf_eff_sig   1.6   1.4   1.7

tcr_eff1 %>%
  ggplot(aes(x = tcr, col = series)) +
  geom_vline(xintercept = 1.55, lwd = .25, lty = 4) +
  geom_line(stat = "density") +
  labs(x = expression(~degree*C), y = "Density") +
  xlim(-1, 3) + 
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf,
           alpha = .2) +
  scale_colour_brewer(palette = "Set1", name = "") +
  scale_linetype(name = "Landuse separate?") +
  theme(
    text = element_text(family = font_type),
    legend.position = "bottom"
  )


#######################################################
## 2) Adjust by effectiveness coefficient distributions

## Rather than adjusting by mean coefficient values (as above), an alternative
## approach is to sample from the relevant t distribution, where the 95% CI is
## given by:
# X +/- t_{df=4, (1-a)=.95}(s/sqrt(n)) = X +/- 2.78(s/sqrt(5))
## Using aerosols as an example, we solve for the std deviation, s, by plugging
## in the mean (1.55) and lower CI (1.05):
# 1.05 = 1.55 - 2.78(s/sqrt(5))
# 0.50 = 2.78(s/sqrt(5))
# 0.50/2.78*sqrt(5) = s
# s = 0.4021705

## Simplifying, we can write this as a function that takes the mean (m) and lower
# 95% CI (lc) as inputs, before returning a single, random draw from the 
# resulting t distribution.
eff_func <- function(m, lc){
  rt(1, df=4)*(m - lc)/2.78 + m
}
# E.g. For aerosols:
eff_func(1.55, 1.05)

## Load climate data
# climate <- read_csv("./Data/climate.csv")

## Rest chain (and rf2x) length to reduce time needed for large loop over the
## t-distribution sampling below
chain_length <- 3000
rf2x <- rf2x[1:chain_length * n_chains]

tcr_eff2 <-
  pblapply(1:10000, function(x){
    
    rcps_eff <-
      rcps %>%
      filter(year >= 1866 & year <= 2100) %>%
      select(year, rcp, trf, volc, solar, co2ch4n2o, fgassum, mhalosum, aerosols,
             cloud_tot, stratoz, tropoz, ch4oxstrath2o, landuse, bcsnow ## "other" category
             ) %>% 
      ## Create effective forcing versions of the necessary variables
      mutate(volc_eff = volc * eff_func(0.61, 0.33),
             solar_eff = solar * eff_func(1.68, -1.27), 
             co2ch4n2o_eff = co2ch4n2o * eff_func(1.17, 1.07),
             # fgassum_eff = fgassum, 
             mhalosum_eff = mhalosum * eff_func(0.66, 0.34),
             aerosols_eff = aerosols * eff_func(1.55, 1.05),
             landuse_eff = landuse * eff_func(3.82, -2.16)
             ) %>%
      mutate(other_eff = landuse_eff + cloud_tot + stratoz + tropoz + ch4oxstrath2o + bcsnow) %>%
      ## Sum these and (other variables) to get effective TRF
      mutate(trf_eff = solar_eff + co2ch4n2o_eff + fgassum + mhalosum_eff + aerosols_eff + other_eff)
    
    ##
    # Join new effective TRF (and volc) variables with main dataset
    clim_df <- 
      bind_cols(
        climate, 
        rcps_eff %>% select(trf_eff, volc_eff)
        ) %>%
      filter(rcp == "rcp26") %>%
      filter(year <= 2005)
    
    theta_sample_eff <- blinreg(clim_df$had, 
                                cbind(alpha = 1, 
                                      beta = clim_df$trf_eff,
                                      gamma = clim_df$volc_eff, 
                                      delta = clim_df$soi, 
                                      eta = clim_df$amo), 
                                chain_length*n_chains)
    
    tcr <- as.data.frame(theta_sample_eff[[1]])$Xbeta * rf2x
    
    tcr_df <- 
      data_frame(tcr, series = x) %>%
      sample_frac(0.05) ## Sample 5% (to reduce size of combined data frame)
    
    rm(theta_sample_eff)
    
    return(tcr_df)
    
  }) %>%
  bind_rows()

tcr_eff2 %>%
  summarise(mean = round(mean(tcr), 1),
            q025 = round(quantile(tcr, .025), 1),
            q975 = round(quantile(tcr, .975), 1))
#  mean  q025  q975
# <dbl> <dbl> <dbl>
#   1.9   0.4   3.4

tcr_eff2 %>%
  # sample_frac(0.1) %>%
  ggplot(aes(x = tcr)) +
  geom_line(stat = "density") +
  labs(x = expression(~degree*C), y = "Density") +
  # xlim(-1, 3) + 
  coord_cartesian(xlim=c(-1, 3)) + 
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf,
           alpha = .2) +
  # stat_function(fun = dnorm, args = list(mean = 0, sd = .065), 
  #               lty=2, col=prior_cols[1]) +
  # stat_function(fun = dnorm, args = list(mean = 0, sd = .25), 
  #               lty=2, col=prior_cols[2]) + 
  # stat_function(fun = dnorm, args = list(mean = 1, sd = .065), 
  #               lty=2, col=prior_cols[3]) +
  # stat_function(fun = dnorm, args = list(mean = 1, sd = .25), 
  #               lty=2, col=prior_cols[4]) + 
  # scale_colour_manual(values = prior_cols) +
  # guides(col = guide_legend(nrow = 2)) +
  theme_tcr 


bind_rows(
  tcr_eff1 %>%
    filter(series == "trf_eff") %>%
    mutate(series = "marvel_a") %>%
    group_by(series) %>%
    summarise(mean = mean(tcr),
              q025 = quantile(tcr, .025),
              q975 = quantile(tcr, .975)),
  tcr_eff2 %>%
    mutate(series = "marvel_b") %>%
    group_by(series) %>%
    summarise(mean = mean(tcr),
              q025 = quantile(tcr, .025),
              q975 = quantile(tcr, .975))
  ) %>%
  write_csv("Robustness/Data/tcr-marvel.csv")

rm(tcr_eff1, tcr_eff2)






# ######################################
# ######################################
# ## Extra:
# 
# ### Compare RCP forcings with that in Marvel et al., who use the historical (1850-2012) 
# # forcings from the GISS CMIP5 simulations as described in Miller et al. (2014).
# ## See: http://data.giss.nasa.gov/modelforce/
# miller <- 
#   tbl_df(read.table("http://data.giss.nasa.gov/modelforce/Fi_Miller_et_al14_upd.txt",
#                     skip = 3, header = T, stringsAsFactors = F))
# ## Keep only ensemble means
# miller <-
#   miller %>%
#   magrittr::set_colnames(gsub("_", "", tolower(colnames(.)))) %>% 
#   rename(ghg = wmghg,
#          aerosols_dir = tropaerdir, 
#          aerosols_ind = tropaerind,
#          mhalosum = ozone,
#          volc = strataer) %>% 
#   filter(year <= 2005)
# 
# miller %>% gather(key, value, -year) %>% ggplot(aes(year, value, col = key)) + geom_line()
# 
# miller$trf <- miller %>% select(-c(year, volc)) %>% rowSums(na.rm=T)
# 
# rcps_base <-
#   rcps %>% 
#   filter(rcp == "rcp26") %>%
#   filter(year >= min(miller$year) & year <= 2005) %>%
#   mutate_if(is.double, funs(. - first(.)))
# 
# # TRF (excl. volc)
# ggplot(rcps_base, aes(year, trf)) + 
#   geom_line(aes(col = "rcps")) + 
#   geom_line(data=miller, aes(col="miller"), lty = 2)
# 
# # Solar
# ggplot(rcps_base, aes(year, solar)) +
#   geom_line(aes(col = "rcps")) + 
#   geom_line(data=miller, aes(col="miller"), lty = 2)
# 
# # Volc
# ggplot(rcps_base, aes(year, volc)) +
#   geom_line(aes(col = "rcps")) + 
#   geom_line(data=miller, aes(col="miller"), lty = 2)
# 
# # GHGs
# ggplot(rcps_base, aes(year, ghg)) +
#   geom_line(aes(col = "rcps_ghg")) + 
#   geom_line(aes(y=kyotoghg, col = "rcps_kyotoghg")) + 
#   geom_line(aes(y=co2ch4n2o, col = "rcps_co2ch4n2o")) +
#   geom_line(data=miller, aes(col="miller"), lty = 2)
# 
# # Land use
# ggplot(rcps_base, aes(year, landuse)) +
#   geom_line(aes(col = "rcps")) + 
#   geom_line(data=miller, aes(col="miller"), lty = 2)
# 
# # Ozone
# ggplot(rcps_base, aes(year, mhalosum)) +
#   geom_line(aes(col = "rcps")) + 
#   geom_line(data=miller, aes(col="miller"), lty = 2)
