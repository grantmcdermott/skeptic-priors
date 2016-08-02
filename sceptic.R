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


## Preallocate coefficients, tcr and temperature in 2100 lists for loop
coefs_tab <- list()
l <- 0 ## count variable for coefs_tab list
tcr <- list()
all_2100 <- list()


ptm <- proc.time()
## Loop over priors ##
for (k in 1:3)  {  
  prior_type <- c("ni", "luke", "den")[k] 

  ## Loop over conviction strength ##
  if(prior_type == "ni")  {
    convic_type <- ""
    # source("jags-loop.R") ## For vague noninformative riors using the rjags package
    source("noninf-loop.R") ## For "proportional" noninformative prors using the LearnBayes package
  }
  else{for (j in 1:2)  {   
    convic_type <- c("mod", "strong")[j]
    source("jags-loop.R")   
  } } ## End of conviction loop

} ## End of prior loop
proc.time() - ptm
# user  system elapsed 
# 28.140   2.912 118.166 


##################################
### COMBINED TABLES AND GRAPHS ###
##################################
pref <- "./TablesFigures/"
suff <- ""

source("sceptic_tablesfigures.R")
