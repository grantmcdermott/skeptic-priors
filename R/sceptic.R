## Choose run type
run_type <- "main"

# library(here)
## Load all packages, as well as some helper functions that will be used for plotting and tables
source(here::here("R/sceptic_funcs.R"))

## Decide on total length of MCMC chains (i.e. summed parallel chains JAGS model)
## Each individual chain will thus be chain_length/n_chains.
chain_length <- 30000
## The below below tries to optimize the number of parallel MCMC chains given
## available CPUs, but balanced against the diminishing returns brought on by
## repeating the burn-in period for each parallel worker. 
n_chains <- n_chains_func(chain_length)

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 

# Run the nested loop (takes about 1min on my laptop)
## Outer: Loop over priors ##
priors_loop <-
  pblapply(1:nrow(priors_df), function(j){
    
    m <- priors_df[j, ]$mu
    s <- priors_df[j, ]$sigma
    prior_type <- priors_df[j, ]$prior_type
    convic_type <- priors_df[j, ]$convic_type 
    
    mu_beta <- m/3.71
    sigma_beta <- s/3.71
    
    ## Inner: Loop over climate scenarios
    if(prior_type == "ni")  {
        # source("R/jags-loop.R", local = T) ## For vague noninformative riors using the rjags package
        source(here("R/noninf-loop.R"), local = T) ## For "proportional" noninformative prors using the LearnBayes package
      }
      else{
        source(here("R/jags-loop.R"), local = T)
      }
    
  })

priors_loop

## Extract only the .$value elements (drop surplus .$visibility element)
priors_loop <- 
  lapply(1:length(priors_loop), function(x) priors_loop[[x]]$value)

priors_loop <- 
  do.call(function(...) mapply(bind_rows, ..., SIMPLIFY = F), args = priors_loop)

## Extract/copy the data frames within the priors_loop list to the (local) global 
# environment and then delete the list itself.
list2env(priors_loop, .GlobalEnv) ## will take extract all the data frames
rm(priors_loop)


##################################
### COMBINED TABLES AND GRAPHS ###
##################################
pref <- here("TablesFigures/")
suff <- ""

source(here("R/sceptic_tablesfigures.R"))
