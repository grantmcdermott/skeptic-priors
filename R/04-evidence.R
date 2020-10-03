modrun = 'evid'
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

set.seed(123)
n_chains = 4
chain_length = 1000
rf2x = rnorm(n_chains * chain_length, mean = 3.71, sd = 0.1855)

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
climate = fread(here('data/climate.csv'))[
	rcp=='rcp60', 
	.(year, had, trf, volc_mean, soi_mean, amo_mean)]

## Load priors data frame
priors_df = fread(here('data/priors.csv'))[prior_type!='ni']

## Simulated GMST from the predicted posterior draws of the main model run.
gmst_sim = fread(here('results/main/gmst-sim.csv'))

## Regression results from main run. As above, will use the noninformative 
## posterior parameters to simulate "true" future values
pparams = fread(here('results/main/params.csv'))[prior=='ni']


# * Priors range ----------------------------------------------------------

mu = seq(max(priors_df$mu), min(priors_df$mu), length = 11)
sigma = round(seq(max(priors_df$sigma), min(priors_df$sigma), length = 11), 3)

priors_range = 
	expand.grid(mu = mu, sigma = sigma) %>%
	as.data.table() %>%
	.[, max_mu := mu == max(mu), by = sigma] ## For 'resetting the clock' below

rm(mu, sigma)


# * Simulate future climate data ------------------------------------------

## Center the vars according to pre 2005 values
climate = center_pre2005(climate)

## Add model noise
climate = merge(climate, gmst_sim, by = 'year')

## Replace simulated values with historical ones where appropriate
climate[year<=2005, gmst_sim := had]

# plot(climate$year, climate$had_sim, type = 'l')

# Evid function -----------------------------------------------------------

plan(multisession, workers = 2)

## Loop counter variables and thresholds
threshold = c(1.3, 1.5)
yr_start_1.3C = min(climate$year) + 60
yr_start_1.5C = min(climate$year) + 90
yrs_start = list(yr_start_1.3C, yr_start_1.5C)
yrs = yrs_start
yr_max = 2100

## DESCRIPTION:
## The function iterates over a range of priors (characterised by mu and sigma). 
## Specifically, it moves along a sequence of (decreasing) sigma values for a 
## given mu. Once it reaches the end of the sigma sequence for that mu, it moves 
## on to the next, lower mu value and repeats the procedure. It does this for 
## two TCR threshold levels: 1.3 °C and 1.5°C. (Note that the function will need 
## to be adjusted  if different thresholds are chosen.) The mean posterior TCR 
## values must be at least as big as the relevant threshold to qualify as fully 
## converged with mainstream. 
evid_func = 
	function() {
		
		pb = progressor(along = 1:nrow(priors_range))
		
    # * Outer loop over thresholds ----------------------------------------
		
		future_lapply(
			seq_along(yrs_start),
			function(i) {
				
				thresh = threshold[i]
				

				# * Inner loop over priors ----------------------------------------
				
				priors_loop = 
					lapply(
						1:nrow(priors_range), 
						function(j) {

							# ** Priors -------------------------------------------------
							
							beta_mu = priors_range$mu[j] / 3.71
							beta_sigma = priors_range$sigma[j] / 3.71
							max_mu = priors_range[j, ]$max_mu
						

							# ** Counter for while loop ---------------------------------

							## NB: Use "<<-" to assign value to global workspace for next iteration
							yrs[[i]] <<- ifelse(max_mu, 
													 			  yrs_start[[i]], ## Major reset for next sigma round
																	yrs[[i]])      ## Minor reset for next mu round


							# ** While loop for TCR convergence -------------------------

							tcr_mean = 0 ## Place holder value
							
							while(round(tcr_mean, 1) < thresh & yrs[[i]] < yr_max) {
								
								## NB: Use "<<-" to assign value to global workspace for next iteration
								yrs[[i]] <<- yrs[[i]] + 1 
								
								clim_df =	climate[year <= yrs[[i]], !"had"]
								

								# *** Specify data list and params ------------------------

								N = nrow(clim_df)
								
								## Standard deviation for weakly informative priors. Here I follow
								## the advice of the Stan core team, using 2.5 * sd(y) / sd(x). The
								## prior means of the centered data are ofc zero.
								wi_sigma = function(gmst_var = 'gmst_sim') { ## CHANGED FOR THIS RUN
									sd_cols = setdiff(names(clim_df), c('year', gmst_var))
									sd_cols = setdiff(sd_cols, names(which(!sapply(clim_df, is.numeric))))
									clim_df[year<=2005 , 
													lapply(.SD, function(x) round(2.5*sd(get(gmst_var))/sd(x), 2)), 
													.SDcols = sd_cols]
									}
								prior_sigmas = wi_sigma('gmst_sim')
								
								data_list = 
									list(
										'N' = N, 
										'gmst' = clim_df$gmst_sim, ## CHANGED
										'trf' = clim_df$trf, 
										'volc' = clim_df$volc_mean, 'soi' = clim_df$soi_mean,
										'amo' = clim_df$amo_mean,
										'beta_mu' = beta_mu, 'beta_sigma' = beta_sigma,
										'gamma_sigma' = prior_sigmas$volc_mean,
										'delta_sigma' = prior_sigmas$soi_mean,
										'eta_sigma' = prior_sigmas$amo_mean
										)
								

								# ** Fit the CmdStanR model -------------------------------
								
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
								

								# ** Posterior draws --------------------------------------
								
								## TCR (mean only)
								beta_draws = as.matrix(fit$draws("beta"))[, 1]
								tcr_mean = mean(rf2x * beta_draws)
								
								# rm(fit, beta_draws, data_list, N); gc()
								
							} ### end of while loop
							
							## Progress bar
							if (i==2) {
								pb()
							}
							
							
							# ** Save results and restart loop --------------------------

							## "Restart" the iteration each time we get back to the "max mu" 
							## value for that sigma loop. (Max mu will always be at mu = 1 for 
							## the historic data, albeit not necessarily for simulated future 
							## data.) Being conservative, we'll start this new loop at the 
							## convergence year for the previous sigma value.
	
							## Again, use "<<-" to assign value to global workspace for the next 
							## iteration.
							if (max_mu) yrs_start[[i]] <<- yrs[[i]]
							
							evid_dt =
								data.table(
									mu = priors_range$mu[j],
									sigma = priors_range$sigma[j],
									tcr_mean = tcr_mean,
									yr_converge = yrs[[i]],
									yr_start = yrs_start[[i]], ## Mostly as a sanity check / efficiency measure
									thresh = thresh
									)
							
							return(evid_dt)
					})
				
				priors_loop = rbindlist(priors_loop)
				
				return(priors_loop)
			},
			future.seed = 123L)
	}


# Run Stan models (with progress bar) -------------------------------------

system.time(with_progress({res = evid_func()}))

## Bind internal data tables
res = rbindlist(res)

res$run = modrun


# Export results ----------------------------------------------------------

res_dir = 'results/evidence'

fwrite(res, here(res_dir, 'evid.csv'))

## Performance
ptime = proc.time() - ptime
mem_gb = as.numeric(system(
	"awk '/MemTotal/ { print $2/1024/1024 }' /proc/meminfo", 	
	intern=TRUE))
perf = data.table(run = modrun, sec = sprintf('%.3f', ptime[[3]]), 
									cores = future::availableCores(),
									mem_gb = round(mem_gb),	os = sessionInfo()$running, 
									arch = sessionInfo()$platform)
fwrite(perf, here('performance', paste0('perf-', modrun, '.csv')))
