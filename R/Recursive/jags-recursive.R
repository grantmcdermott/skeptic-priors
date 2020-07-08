N <- nrow(clim_df)

##------------------------------------------------------------------------------
## THE BUGS/JAGS MODEL.

bugs_file <- bugs_model_func


##------------------------------------------------------------------------------
## SPECIFY THE DATA, INITIALIZATION VALUES AND PARAMETERS OF INTEREST.

## Tell JAGS where the data are coming from
## Notes: Using "_sim" versions in case of future recursive type. 
data_list <- 
  list(
    "N" = N, "had" = clim_df$had_sim, "trf" = clim_df$trf, 
    "volc" = clim_df$volc_sim, "soi" = clim_df$soi_sim, "amo" = clim_df$amo_sim,
    "mu_beta" = mu_beta, "sigma_beta" = sigma_beta
    )

## Give JAGS some initialization values for the model parameters
inits_list <- 
  function() {
    list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1, phi = 0.5)
    }

## Which parameters should R keep track of (i.e. return the posterior distributions for)?
parameters <- c("beta") ## Only need beta to calculate TCR

##------------------------------------------------------------------------------
## RUN THE PARALLEL JAGS MODEL.

mod_samples <- 
  jags_par_model(
    bugs_file=bugs_file, data_list=data_list, inits_list=inits_list, parameters=parameters
    )

##------------------------------------------------------------------------------
## Get TCR coefficient (i.e. beta) and then use to get TCR.

tcr <- 
  (as.matrix(mod_samples[, "beta"], iters = F) * rf2x) %>%
  as_tibble() %>%
  magrittr::set_colnames("tcr") %>%
  mutate(
    prior = paste0(prior_type, convic_type),
    year_to = ifelse(recurse_type == "historic", yr_min, max(clim_df$year))
    )

return(tcr = tcr)