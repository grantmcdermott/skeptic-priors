modrun = 'main'

# Libraries ---------------------------------------------------------------

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(hrbrthemes)
library(data.table)
library(magrittr)
library(fst)
library(future.apply)
library(progressr)
library(here)

theme_set(theme_ipsum())

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
mod = cmdstan_model(here('stan/mod-pred.stan'))

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

## Load climate data
climate = fread(here("data/climate.csv"))

## Load priors data frame
priors_df = fread(here("data/priors.csv"))

## RCPS (might as well name them here)
rcps = c('rcp26', 'rcp45', 'rcp60', 'rcp85')


# Set up loop function ----------------------------------------------------

plan(multisession, workers = floor(availableCores()/n_chains))

priors_loop = function() {
  pb = progressor(along = 1:nrow(priors_df))
  future_lapply(
    1:nrow(priors_df),
    function(j) {

      # * Priors ----------------------------------------------------------      
      
      beta_mu = priors_df$mu[j]
      beta_sigma = round(priors_df$sigma[j], 3)
      prior_type = priors_df$prior_type[j]
      convic_type = priors_df$convic_type[j]
      prior_convic = paste0(prior_type, convic_type)
      
      ## Progress bar
      pb(sprintf("Prior = %s", prior_convic), class = "sticky")
      
      ## adjust subjective TCR priors to climate resistance scale 
      if (j!=1) {
        beta_mu = beta_mu/3.71
        beta_sigma = beta_sigma/3.71
      }
      

      # * Inner RCP loop --------------------------------------------------

      rcp_loop =
        lapply(
          rcps,
          function(i) {
              
            ## Which RCP
            clim_df = climate[rcp==i, .(year, had, trf, volc_mean, soi_mean, amo_mean)]
            
            ## Center the vars according to pre 2005 values
            clim_df = center_pre2005(clim_df)
            

            # ** Specify data list and params -----------------------------
            
            N1 = nrow(clim_df[year <= 2005]) 
            N2 = nrow(clim_df)
            
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
                  'N1' = N1, 'N2' = N2, 
                  'gmst' = clim_df[year <= 2005, had], 'trf' = clim_df$trf, 
                  'volc' = clim_df$volc_mean, 'soi' = clim_df$soi_mean, 
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
            
            ## Params of interest
            params = c("alpha", "beta", "gamma", "delta", "eta", "phi")
    
            
            # ** Posterior draws and diagnostics --------------------------
            
            tcr = NULL
            params_tab = NULL
            
            ## Only needed for a single RCP run, since we're using historic info.
            ## Obviously doesn't make a difference which, but I'll go with the 
            ## upper-middle-of-the-road RCP 6.0 for consistency elsewhere in the
            ## paper.
            if (i=='rcp60') {
              
              ## Diagnostics
              # fit$cmdstan_diagnose() ## helpful English summaries, but takes a while to run...
              sdg = summary(fit$sampler_diagnostics())
              fwrite(sdg, here('stan/diagnostics/main', paste0(prior_convic, '.csv')))
              
              color_scheme_set('blue'); theme_set(theme_ipsum())
              mcmc_trace(fit$draws(params), facet_args = list(labeller = label_parsed)) +
                ggsave(here('stan/diagnostics/main', paste0('traceplot-', prior_convic, '.png')),
                       width = 8, height = 5)

              ## Posterior densities
              color_scheme_set("darkgray");theme_set(theme_ipsum())
              mcmc_dens(fit$draws(params),
                        pars = params, alpha = 0.3,
                        facet_args = list(labeller = label_parsed, ncol = 2)) +
                ggsave(here('figs/main', paste0('params-', prior_convic, '.png')),
                       width = 8, height = 10)
              
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
              # num_cols = c('mean', 'q025', 'q975')
              # params_tab[ , (num_cols) := lapply(.SD, function(x) sprintf('%.3f', x)), .SDcols = num_cols]
            
            }
            
            ## Predicted temps
            pred_params = paste0("y_pred[", 1:nrow(clim_df), "]")
            gmst_pred = 
              data.table(cbind(
                year = seq(min(clim_df$year), length.out = nrow(clim_df)),
                fit$summary(pred_params, "mean", ~quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE))
                ))[, variable := NULL][]
            setnames(gmst_pred, old = c('2.5%', '97.5%'), new = c('q025', 'q975'))
            gmst_pred$prior = prior_convic
            
            ## 2100 temps (full distributions) 
            gmst2100 = data.table(gmst = as.matrix(fit$draws("y_pred[235]"))[, 1])
            gmst2100$rcp = i
            gmst2100$prior = prior_convic
        
            return(list(tcr = tcr, 
                        params_tab = params_tab, 
                        gmst_pred = gmst_pred,
                        gmst2100 = gmst2100))
            }
          )
      
      # Recombine the sub-elements of the list based on their common indexes
      rcp_loop =
        do.call(function(...) {
          mapply(rbind, ..., SIMPLIFY = FALSE)
          },
          args = rcp_loop)
      
      # ## Progress bar
      pb(sprintf("Prior = %s", prior_convic), class = "sticky")
      
      return(rcp_loop)
      })
}


# Run Stan models (with progress bar) -------------------------------------

system.time(with_progress({res = priors_loop()}))

## Recombine the sub-elements of the list based on their common indexes
res =
  do.call(function(...) mapply(rbind, ..., SIMPLIFY = FALSE), args = res)

## Get deviations of mean fitted values with actual HadCRUT obs. Will use this
## series as noise when simulating the future "true" temperatures in the
## "evidence" section of the code/paper.
had_dev =
  merge(
    climate[rcp=='rcp60' & !is.na(had_full), .(year, had_full)],
    res$gmst_pred[, .(year, mean)]
  ) %>%
  .[, .(year, had_dev = mean - had_full)]

# Export results ----------------------------------------------------------

## Note: Will compress large files (>5k rows) using FST format

res_dir = 'results/main'

res$tcr$run = modrun
res$gmst2100$run = modrun
res$params_tab$run = modrun
res$gmst_pred$run = modrun

write_fst(res$tcr, here(res_dir, 'tcr.fst'))
write_fst(res$gmst2100, here(res_dir, 'gmst2100.fst'))
fwrite(res$gmst_pred, here(res_dir, 'gmst-pred.csv'))
fwrite(res$params_tab, here(res_dir, 'params.csv'))
fwrite(had_dev, here(res_dir, 'had-dev.csv'))
