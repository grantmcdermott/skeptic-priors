rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("Data/climate.csv")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
## Total chain length will thus be chain_length * n_chains
chain_length <- 10000
n_chains <- detectCores() - 1 

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 

## Priors data frame
priors_df <- 
  data_frame(mu = c(0, 1, 1, 0, 0),
             sigma = c(100, 0.25, 0.065, 0.25, 0.065),
             prior_type = c("ni", "luke", "luke", "den", "den"),
             convic_type = c("", "mod", "strong", "mod", "strong")
             )
priors_df

# Run the nested loop (takes about 2min on my laptop)
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
        # source("jags-loop.R", local = T) ## For vague noninformative riors using the rjags package
        source("noninf-loop.R", local = T) ## For "proportional" noninformative prors using the LearnBayes package
      }
      else{
        source("jags-loop.R", local = T)
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
pref <- "TablesFigures/"
suff <- ""

source("sceptic_tablesfigures.R")
