# modrun = 'eff'
ptime = proc.time()

# Libraries ---------------------------------------------------------------

library(cmdstanr)
library(posterior)
library(data.table)
library(magrittr)
library(fst)
library(future.apply)
library(progressr)
library(here)

## Progress bar options and output
options(progressr.enable = TRUE) ## For batch mode
handlers(
	handler_progress(
		format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
		width    = 60,
		complete = "+"
		)
	)

# Global variables and functions ------------------------------------------

## Stan model (will compile first time it is run)
mod = cmdstan_model(here('stan/mod.stan'))

# set.seed(123)
# n_chains = 4
# chain_length = 1000
# rf2x = rnorm(n_chains * chain_length, mean = 3.71, sd = 0.1855)

center_pre2005 = function(d) {
	d = as.data.table(d)
	d05 = d[year<=2005, .SD, .SDcols = is.numeric][, !'year']
	cols05 = colnames(d05)
	colsother = setdiff(colnames(d), cols05)
	dnums = d[, ..cols05]
	dnums = sweep(dnums, 2, apply(d05, 2, mean, na.rm = TRUE))
	ret = cbind(d[, ..colsother], dnums)
}


# Load data ---------------------------------------------------------------

## Load climate data.
climate = fread(here("data/climate.csv"))[
	rcp=='rcp60' & year <= 2005, 
	.(year, had, trf, volc_mean, soi_mean, amo_mean)]
setkey(climate, year)

## Load RCP data
rcps = fread(here('data/raw/rcps.csv'))[rcp=='rcp60']
colnames(rcps) = gsub("_rf", "", tolower(colnames(rcps)))
rcps = rcps[, .(year = years, 
								# trf = total_inclvolcanic - volcanic_annual, 
								volc = volcanic_annual,
								solar,
								co2ch4n2o, fgassum, 
								mhalosum, 
								aerosols = totaer_dir,
								cloud_tot, stratoz, tropoz, ch4oxstrath2o, landuse, bcsnow) ## other
						] 


## Load priors data frame
priors_df = fread(here("data/priors.csv"))



# Background: Marvel et al. (2016) ----------------------------------------

## Marvel et al. (2016), hereafter M2016, provide the following scaled 
## efficacies (see table S1 in SI):
## http://data.giss.nasa.gov/modelforce/Marvel_etal2015.html
## Aerosols: 1.55* (1.05, 2.05)
## GHG:      1.17* (1.06, 1.27)
## Land use: 4.27 (-2.42, 10.95)
## Ozone:    0.66* (0.34, 0.98)
## Solar:    1.68 (-1.27, 4.63)
## Volc:     0.61* (0.33, 0.89)
## historical: 1.00 (0.83, 1.16)

## Note: efficacies are relative to CO2 (= 1)
## "The uncertainty in the efficacies is estimated from individual members of the
## single-forcing ensembles (Figure S2). Confidence intervals on the sample mean are
## constructed using a student-t distribution with 4 degrees of freedom (5 in the case
## of the 6-member historical ensemble)."


# 1) Means-only efficacies ------------------------------------------------

## Adjust the relevant columns to match up with the M2016 forcings. Note that 
## the forcings in M2016 don't match up exactly to those in the RCP dataset. For
## example, they have a forcing called "ozone", which I assume corresponds to
## the RCP forcing "mhalosum" (i.e. gases controlled under the Montreal 
## Protocol). See the RCP metadata and flowchart in the data folder for more 
## information. Summary version: We construct a new effective TRF variable 
## using the following variables:
## TRF_eff (excl. volc) = 
##   1.68*solar +
##   1.17*co2ch4n2o + ## "well-mixed GHGs" in MEA2016
##   fgassum + ## flourinated gases, not adjusted
##   0.66*mhaslosum + ## "ozone" in MEA2016
##   1.55*aerosols + ## "direct aerosols" in MEA16
##   other (incl. 4.27*landuse)

# * Globals ---------------------------------------------------------------

modrun = 'eff1'

seed(123)
n_chains = 4
chain_length = 1000
rf2x = rnorm(n_chains * chain_length, mean = 3.71, sd = 0.1855)


# * Adjust forcing efficacies ---------------------------------------------

rcps_eff1 = copy(rcps)  

## Create effective forcing versions of the necessary variables
rcps_eff1[, ':=' (volc_eff      = 0.61 * volc,
							 	  solar_eff     = 1.68 * solar, 
								  co2ch4n2o_eff = 1.17 * co2ch4n2o, 
								  mhalosum_eff  = 0.66 * mhalosum, 
								  aerosols_eff  = 1.55 * aerosols,
								  landuse_eff   = 4.27 * landuse)
				 ] %>%
	.[, other_eff := landuse_eff + cloud_tot + stratoz + tropoz + ch4oxstrath2o + bcsnow] %>%
	## Sum these and (other variables) to get effective TRF
	.[, trf_eff := solar_eff + co2ch4n2o_eff + fgassum + mhalosum_eff + aerosols_eff + other_eff]
setkey(rcps_eff1, year)

## Join new effective TRF (and volc) variables with main dataset
clim_df = merge(climate, 
								rcps_eff1[, .(year, trf_eff, volc_eff)],
								by = 'year')
## Compare
# library(ggplot2)
# melt(clim_df[, .(year, trf, trf_eff)], 
# 		 id.vars = 'year', variable.name = 'key', value.name = 'trf') %>%
# 	ggplot(aes(x = year, y = trf, lty = key)) +
# 	geom_line() 


# * Set up loop function --------------------------------------------------

plan(multisession, workers = floor(availableCores()/n_chains))

priors_loop = function() {
	pb = progressor(along = 1:nrow(priors_df))
	future_lapply(
		1:nrow(priors_df),
		function(j) {
			
			# ** Priors ---------------------------------------------------------      
			
			beta_mu = priors_df$mu[j]
			beta_sigma = round(priors_df$sigma[j], 3)
			prior_type = priors_df$prior_type[j]
			convic_type = priors_df$convic_type[j]
			prior_convic = paste0(prior_type, convic_type)
			
			## Progress bar
			pb(sprintf("Prior = %s", prior_convic))
			
			## adjust subjective TCR priors to climate resistance scale 
			if (j!=1) {
				beta_mu = beta_mu/3.71
				beta_sigma = beta_sigma/3.71
			}
			
			
			# clim_df = climate[, .(year, had, trf, volc_mean, soi_mean, amo_mean)]
			
			## Center the vars according to pre 2005 values
			clim_df = center_pre2005(clim_df)
			
			
			# ** Specify data list and params -----------------------------
			
			N = nrow(clim_df)
			
			## Standard deviation for weakly informative priors. Here I follow
			## the advice of the Stan core team, using 2.5 * sd(y) / sd(x). The
			## prior means of the centered data are ofc zero.
			wi_sigma = function(gmst_var = 'had') {
				sd_cols = setdiff(names(clim_df), c('year', gmst_var))
				sd_cols = setdiff(sd_cols, names(which(!sapply(clim_df, is.numeric))))
				clim_df[year<=2005 , 
								lapply(.SD, function(x) round(2.5*sd(get(gmst_var))/sd(x), 2)), 
								.SDcols = sd_cols]
			}
			prior_sigmas = wi_sigma()
			
			data_list = 
				list(
					'N' = N, 
					'gmst' = clim_df$had,
					'trf' = clim_df$trf_eff,   ## CHANGED
					'volc' = clim_df$volc_eff, ## CHANGED
					'soi' = clim_df$soi_mean, 'amo' = clim_df$amo_mean,
					'beta_mu' = beta_mu, 'beta_sigma' = beta_sigma,
					'gamma_sigma' = prior_sigmas$volc_mean,
					'delta_sigma' = prior_sigmas$soi_mean,
					'eta_sigma' = prior_sigmas$amo_mean
					)
			
			
			# ** Fit the CmdStanR model -----------------------------------
			
			## Model fit
			fit = 
				mod$sample(
					data = data_list,
					seed = 123,
					chains = n_chains,
					parallel_chains = n_chains,
					iter_sampling = chain_length,
					refresh = 0,
					show_messages = FALSE
					)
			
			
			## Params of interest
			params = c("alpha", "beta", "gamma", "delta", "eta", "phi")
			
			
			# ** Posterior draws ------------------------------------------
			
			## TCR
			beta_draws = as.matrix(fit$draws("beta"))[, 1]
			tcr = data.table(tcr = rf2x * beta_draws, prior = prior_convic)
			
			## Summary table
			params_tab = 
				rbind(
					data.table(fit$summary(params, "mean", ~quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE))),
					tcr[, .('tcr', mean(tcr), quantile(tcr, .025, na.rm = TRUE), quantile(tcr, .975, na.rm = TRUE))],
					use.names = FALSE
					)
			setnames(params_tab, 
							 old = c('variable','2.5%', '97.5%'), 
							 new = c('param', 'q025', 'q975'))
			params_tab$prior = prior_convic
			
			
			return(list(params_tab = params_tab, tcr = tcr))
		},
		future.seed = 123L
	)
}


# * Run Stan models (with progress bar) -----------------------------------

system.time(with_progress({res1 = priors_loop()}))

## Recombine the sub-elements of the list based on their common indexes
res1 =
	do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), args = res1)

res1$tcr$run = modrun
res1$params_tab$run = modrun


# * Export results --------------------------------------------------------

res_dir = 'results/sensitivity'

write_fst(res1$tcr, here(res_dir, paste0('tcr-', modrun,'.fst')))
fwrite(res1$params_tab, here(res_dir, paste0('params-', modrun,'.csv')))

rm(clim_df)



# 2) Full distribution efficacies -----------------------------------------


# * Globals ---------------------------------------------------------------

modrun = 'eff2'

set.seed(123)
n_chains = 1 #4
chain_length = 1000
rf2x = rnorm(n_chains * chain_length, mean = 3.71, sd = 0.1855)


# * Adjust forcing efficacies ---------------------------------------------

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
eff_func = 
	function(m, lc){
		rt(1, df=4)*(m - lc)/2.78 + m
	}
## E.g. For aerosols:
# eff_func(1.55, 1.05)

## Bind together in a 1,000 run simulation, taking a new draw from the 
## underlying distribution each time.
rcps_eff2 =
	rbindlist(lapply(
		1:1000,
		function(i) {
			d = copy(rcps)[year >= min(climate$year) & year <= max(climate$year)]
			
			## Create effective forcing versions of the necessary variables
			d[, ':=' (volc_eff      = volc * eff_func(0.61, 0.33),
								solar_eff     = solar * eff_func(1.68, -1.27), 
								co2ch4n2o_eff = co2ch4n2o * eff_func(1.17, 1.07), 
								mhalosum_eff  = mhalosum * eff_func(0.66, 0.34), 
								aerosols_eff  = aerosols * eff_func(1.55, 1.05),
								landuse_eff   = landuse * eff_func(3.82, -2.16))
			] %>%
				.[, other_eff := landuse_eff + cloud_tot + stratoz + tropoz + ch4oxstrath2o + bcsnow] %>%
				## Sum these and (other variables) to get effective TRF
				.[, trf_eff := solar_eff + co2ch4n2o_eff + fgassum + mhalosum_eff + aerosols_eff + other_eff] %>%
				.[, sim := i]
			d
		}
	))

# rcps_eff2
setkey(rcps_eff2, sim, year)


# * Set up simulation function --------------------------------------------

plan(multisession, workers = floor(availableCores()/n_chains))
sims = unique(rcps_eff2[, sim])

sim_loop = function() {
	
	pb = progressor(steps = length(sims)/10)
	
	future_lapply(
		sims,
		function(i) {
			
			## Progress bar (only show every 10th step)
			if (i %% 10 == 0) pb()
			# pb()
			
			## Which simulation?
			trf_eff2 = rcps_eff2[sim == i]
			
			## Join new effective TRF (and volc) variables with main dataset
			clim_df = merge(climate, 
											trf_eff2[, .(year, trf_eff, volc_eff)],
											by = 'year')
			
			## Center the vars according to pre 2005 values
			clim_df = center_pre2005(clim_df)

						
			# Inner prior loop --------------------------------------------------
			
			priors_loop = 
				lapply(
					1:nrow(priors_df),
					function(j) {
						
						# ** Priors ---------------------------------------------------
						
						beta_mu = priors_df$mu[j]
						beta_sigma = round(priors_df$sigma[j], 3)
						prior_type = priors_df$prior_type[j]
						convic_type = priors_df$convic_type[j]
						prior_convic = paste0(prior_type, convic_type)
						
						## adjust subjective TCR priors to climate resistance scale 
						if (j!=1) {
							beta_mu = beta_mu/3.71
							beta_sigma = beta_sigma/3.71
						}
						
						
						# ** Specify data list and params -----------------------------
						
						N = nrow(clim_df)
						
						## Standard deviation for weakly informative priors. Here I follow
						## the advice of the Stan core team, using 2.5 * sd(y) / sd(x). The
						## prior means of the centered data are ofc zero.
						wi_sigma = function(gmst_var = 'had') {
							sd_cols = setdiff(names(clim_df), c('year', gmst_var))
							sd_cols = setdiff(sd_cols, names(which(!sapply(clim_df, is.numeric))))
							clim_df[year<=2005 , 
											lapply(.SD, function(x) round(2.5*sd(get(gmst_var))/sd(x), 2)), 
											.SDcols = sd_cols]
							}
						prior_sigmas = wi_sigma('had')
						
						data_list = 
							list(
								'N' = N, 
								'gmst' = clim_df$had, 
								'trf' = clim_df$trf_eff,   ## CHANGED
								'volc' = clim_df$volc_eff, ## CHANGED 
								'soi' = clim_df$soi_mean, 
								'amo' = clim_df$amo_mean,
								'beta_mu' = beta_mu, 'beta_sigma' = beta_sigma,
								'gamma_sigma' = prior_sigmas$volc_mean,
								'delta_sigma' = prior_sigmas$soi_mean,
								'eta_sigma' = prior_sigmas$amo_mean
								)
						
						
						# ** Fit the CmdStanR model -----------------------------------
						
						## Model fit
						fit = 
							mod$sample(
								data = data_list,
								seed = 123,
								chains = n_chains,
								parallel_chains = n_chains,
								iter_sampling = chain_length,
								refresh = 0,
								show_messages = FALSE
							)
						
						# ** Posterior draws ------------------------------------------
						
						## TCR
						beta_draws = as.matrix(fit$draws("beta"))[, 1]
						tcr = data.table(tcr = rf2x * beta_draws, prior = prior_convic, 
														 sim = i)
						
						return(tcr = tcr)
					}
				)
			
			# Recombine the sub-elements of the list based on their common indexes
			priors_loop = rbindlist(priors_loop)
			
			
			return(priors_loop)
		},
		future.seed = 123L
	)
	}


# * Run Stan models (with progress bar) -----------------------------------

system.time(with_progress({res2 = sim_loop()}))

## Bind internal data tables
res2 = rbindlist(res2)

res2$run = modrun


# * Export results --------------------------------------------------------

res_dir = 'results/sensitivity'

write_fst(res2, here(res_dir, paste0('tcr-', modrun,'.fst')))


# Performance -------------------------------------------------------------

ptime = proc.time() - ptime
mem_gb = as.numeric(system(
	"awk '/MemTotal/ { print $2/1024/1024 }' /proc/meminfo", 	
	intern=TRUE))
perf = data.table(run = modrun, sec = sprintf('%.3f', ptime[[3]]), 
									cores = future::availableCores(),
									mem_gb = round(mem_gb),	os = sessionInfo()$running, 
									arch = sessionInfo()$platform)
fwrite(perf, here('performance', paste0('perf-', modrun, '.csv')))
