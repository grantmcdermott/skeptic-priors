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
	.(year, had_sim = had, trf, volc_sim = volc_mean, soi_sim = soi_mean, 
		amo_sim = amo_mean)]

## Load priors data frame
priors_df = fread(here('data/priors.csv'))[prior_type!='ni']

## Model-derived standard errors (deviations) to needed to simulate "true" 
## future values.
had_dev = fread(here('results/main/had-dev.csv'))

## Regression results from main run. As above, will use the noninformative 
## posterior parameters to simulate "true" future values
pparams = fread(here('results/main/params.csv'))[prior=='ni']

## Recursive TCR estimates for four main sceptic types. Will use the results 
## from the weakest sceptic (i.e. moderate lukewarmer) to set the starting 
## point(s) of the while iteration below. 
# tcr_rec = fread(here('results/recursive/tcr-rec.csv'))[prior=='lukemod']
# yr_start_1.3C = 2005 - tcr_rec[round(tcr_mean, 1)==1.3, max(year_to)] + 1866
# yr_start_1.5C = 2005 - tcr_rec[round(tcr_mean, 1)==1.5, max(year_to)] + 1866
yr_start_1.3C = min(climate$year) + 60
yr_start_1.5C = min(climate$year) + 90

# * Priors range ----------------------------------------------------------

mu = seq(max(priors_df$mu), min(priors_df$mu), length = 11)
sigma = round(seq(max(priors_df$sigma), min(priors_df$sigma), length = 11), 3)
# threshold = c(1.3, 1.5)

priors_range = 
	expand.grid(mu = mu, sigma = sigma) %>%
	as.data.table() %>%
	.[, max_mu := mu == max(mu), by = sigma] ## For 'resetting the clock' below

rm(mu, sigma)


# * Simulate future climate data ------------------------------------------

## Add model noise
climate$noise = sample(had_dev$had_dev, nrow(climate), replace = TRUE)

## Also add noise to future covariate values; otherwise will cause collinearity 
## problems for future regressions (b/c values stay constant).
sd_cols = c('volc_sim', 'soi_sim', 'amo_sim')
climate[,n := .I
				][year > 2005,
					(sd_cols) := lapply(.SD, `+`, rnorm(1, mean = 0, sd = .1)),
					.SDcols = sd_cols,
					by = n
				][]
climate[, n:=NULL]; rm(sd_cols)

## Now simulate future climate data based on noninformative parameters
## Note that we need a shift parameter to adjust the "intercept".
# had05 = climate[year==2005, had_sim]
climate[year > 2005, 
				had_sim := 
					pparams[param=='alpha', mean] + #had05 +
					pparams[param=='beta', mean] * trf + 
					pparams[param=='gamma', mean] * volc_sim + 
					pparams[param=='delta', mean] * soi_sim + 
					pparams[param=='eta', mean] * amo_sim +
					noise
				]

# plot(climate$year, climate$had_sim, type = 'l')

# Evid function -----------------------------------------------------------

plan(multisession, workers = 2)

# yrs_j = list(thresh_1.3C = start_1.3C, thresh_1.5C = start_1.5C)
yrs_start = list(yr_start_1.3C, yr_start_1.5C)
# yrs = list(NULL, NULL)
yrs = yrs_start
threshold = c(1.3, 1.5)
yr_max = 2100

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

							# if (max_mu) {
							# 	yrs[[i]] <<- yrs_start[[i]]  ## Major reset for next sigma round
							# } else {
							# 	yrs[[i]] <<- yrs[[i]]       ## Minor reset for next mu round
							# }
							yrs[[i]] <<- ifelse(max_mu, 
													 			  yrs_start[[i]], ## Major reset for next sigma round
																	yrs[[i]])      ## Minor reset for next mu round


							# ** While loop for TCR convergence -------------------------

							tcr_mean = 0 ## Place holder value
							
							while(round(tcr_mean, 1) < thresh & yrs[[i]] < yr_max) {
								
								## NB: Use "<<-" to assign value to global workspace for next iteration
								yrs[[i]] <<- yrs[[i]] + 1 
								
								clim_df =	climate[year <= yrs[[i]]]
								
								## Center the vars according to pre 2005 values
								clim_df = center_pre2005(clim_df)
								

								# *** Specify data list and params ------------------------

								N = nrow(clim_df)
								
								## Standard deviation for weakly informative priors. Here I follow
								## the advice of the Stan core team, using 2.5 * sd(y) / sd(x). The
								## prior means of the centered data are ofc zero.
								wi_sigma = function(gmst_var = 'had_sim') { ## CHANGED FOR THIS RUN
									sd_cols = setdiff(names(clim_df), c('year', gmst_var))
									sd_cols = setdiff(sd_cols, names(which(!sapply(clim_df, is.numeric))))
									clim_df[year<=2005 , 
													lapply(.SD, function(x) round(2.5*sd(get(gmst_var))/sd(x), 2)), 
													.SDcols = sd_cols]
									}
								prior_sigmas = wi_sigma('had_sim')
								
								data_list = 
									list(
										'N' = N, 
										'gmst' = clim_df$had_sim,                           ## CHANGED
										'trf' = clim_df$trf, 
										'volc' = clim_df$volc_sim, 'soi' = clim_df$soi_sim, ## CHANGED 
										'amo' = clim_df$amo_sim,                            ## CHANGED
										'beta_mu' = beta_mu, 'beta_sigma' = beta_sigma,
										'gamma_sigma' = prior_sigmas$volc_sim,              ## CHANGED
										'delta_sigma' = prior_sigmas$soi_sim,               ## CHANGED
										'eta_sigma' = prior_sigmas$amo_sim                  ## CHANGED
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
