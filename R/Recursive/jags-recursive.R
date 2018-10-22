N <- nrow(clim_df)

##------------------------------------------------------------------------------
## THE BUGS/JAGS MODEL.

bugs_file <- bugs_model_func

##------------------------------------------------------------------------------
## PARALLEL SETUP.

cl <- parallel::makeCluster(n_chains, type = cl_type) # no. of clusters (i.e. MCMC chains)
parLoadModule(cl, "lecuyer", quiet = T)
clusterSetRNGStream(cl, 123)

##------------------------------------------------------------------------------
## SPECIFY THE DATA AND INITIALIZE THE CHAINS.
## Notes: Using "_sim" versions in case of future recursive type. 

data_list <- 
  list(
    "N" = N, "had" = clim_df$had_sim, "trf" = clim_df$trf, 
    "volc" = clim_df$volc_sim, "soi" = clim_df$soi_sim, "amo" = clim_df$amo_sim,
    "mu_beta" = mu_beta, "sigma_beta" = sigma_beta
    )
inits_list <- 
  function() {
    list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1, phi = 0)
    }

##------------------------------------------------------------------------------
## RUN THE CHAINS/MCMC SAMPLES.

## Which parameters should R keep track of?
# parameters <- c("alpha", "beta", "gamma", "delta", "eta", "sigma", "y_pred")
parameters <- c("beta") ## Only need beta to calculate TCR
## Initialisation                
par_inits <- parallel.inits(inits_list, n.chains = n_chains) 
## Create the JAGS model object
parJagsModel(
  cl, name = "jags_mod", file = bugs_file, 
  data = data_list, inits = par_inits, n.chains = n_chains, n.adapt = 1000,
  quiet = T
  )
## Burn-in
parUpdate(cl, "jags_mod", n.iter = 1000, progress.bar = "none")
## Now we run the full model samples
mod_iters <- chain_length/n_chains
mod_samples <- 
  parCodaSamples(
    cl, "jags_mod", variable.names = parameters, 
    n.iter = mod_iters, n.chain = n_chains,
    progress.bar = "none"
    ) 

parallel::stopCluster(cl)

##------------------------------------------------------------------------------
## Get TCR coefficient (i.e. beta) and then use to get TCR.

tcr <- 
  (as.matrix(mod_samples[, "beta"], iters = F) * rf2x) %>%
  as_data_frame() %>%
  magrittr::set_colnames("tcr") %>%
  mutate(
    prior = paste0(prior_type, convic_type),
    year_to = ifelse(recurse_type == "historic", yr_min, max(clim_df$year))
    )

return(tcr = tcr)