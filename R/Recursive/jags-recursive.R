##------------------------------------------------------------------------------
## THE BUGS/JAGS MODEL. 
## Notes: Using "_sim" versions in case of future recursive type. Also dropping 
## y_pred variable b/c don't need predictions for recursive estimates.

N <- nrow(clim_df)

mod_string <- paste(
  "model{
    
    for(t in 1:N) {
      mu[t] <- alpha + beta*trf[t] + gamma*volc_sim[t] + delta*soi_sim[t] + eta*amo_sim[t]
      had_sim[t]  ~ dnorm(mu[t], tau)
    }
    ", 
  paste0(
    "
            mu_beta <- ", mu_beta),
  paste0(
    "
            sigma_beta <- ", sigma_beta),
  
    "

    ## Priors for all parameters   
    alpha ~ dnorm(0, 0.0001)            ## intercept
    beta ~ dnorm(mu_beta, tau_beta)     ## trf coef
      tau_beta <- pow(sigma_beta, -2)
    gamma ~ dnorm(0, 0.0001)            ## volc coef
    delta ~ dnorm(0, 0.0001)            ## soi coef
    eta ~ dnorm(0, 0.0001)              ## amo coef
    sigma ~ dunif(0, 100)               ## Residual std dev
      tau <- pow(sigma, -2)  	          
    had0 ~ dnorm(0.0, 1.0E-6)           ## Initialising value
    }" 
) 
bugs_file <- paste0("BUGSFiles/Recursive/", prior_type, "-", convic_type, recurse_type, ".txt")
if(prior_type == "ni"){bugs_file <- gsub("--","-",bugs_file)}
writeLines(mod_string, con = bugs_file)

load.module("lecuyer") ## JAGS module uses lecuyer random number generator (to avoid overlap/correlation in a parallel format)

cl <- parallel::makeCluster(n_chains, type = cl_type) # no. of clusters (i.e. MCMC chains)
parLoadModule(cl, "lecuyer", quiet = T)

##------------------------------------------------------------------------------
## INTIALIZE THE CHAINS.

data_list <- list("N" = N, "had_sim" = clim_df$had_sim, "trf" = clim_df$trf, 
                  "volc_sim" = clim_df$volc_sim, 
                  "soi_sim" = clim_df$soi_sim, "amo_sim" = clim_df$amo_sim)

inits_list <- function() {
  list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1)
}

##------------------------------------------------------------------------------
## RUN THE CHAINS.

parameters <- c("alpha", "beta", "gamma", "delta", "eta", "sigma")
par_inits <- parallel.inits(inits_list, n.chains = n_chains) # Initialisation

parJagsModel(cl, name = "jags_mod", file = bugs_file, 
             data = data_list, inits = par_inits, n.chains = n_chains, n.adapt = 1000)
parUpdate(cl, "jags_mod", n.iter = 1000) # burn-in
mod_iters <- chain_length/n_chains
mod_samples <- parCodaSamples(cl, "jags_mod", variable.names = parameters, 
                              n.iter = mod_iters, n.chain = n_chains)
parallel::stopCluster(cl)

##------------------------------------------------------------------------------
## Get TCR coefficient (i.e. beta) and then use to get TCR.

## Convert coefficients MCMC list into data frame for later. First combines
## all chains into one matrix.
coefs_df <-
  as.matrix(mod_samples[, c(1:6)], iters = F) %>%
  data.frame() %>% 
  tbl_df() %>%
  gather(coef, values)

## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
tcr <- 
  data_frame(
    beta = filter(coefs_df, coef == "beta")$values,
    prior = paste0(prior_type, convic_type)
    ) %>%
  mutate(year_to = ifelse(recurse_type == "historic", yr_min, max(clim_df$year)))

tcr$tcr <- tcr$beta * rf2x

tcr <- tcr %>% select_("-beta")

return(tcr = tcr)