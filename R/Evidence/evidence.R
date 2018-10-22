## NB: This script is best run from the shell ("~$ Rscript evidence.R"), *not* from
## within the RStudio IDE. See description of `evid_func` below (+/- line 55).

## Choose run type
run_type <- "evidence"

# library(here)
## Load all packages, as well as some helper functions that will be used for plotting and tables
source(here::here("R/sceptic_funcs.R"))

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv(here("Data/climate.csv"))

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
chain_length <- 9000
n_chains <- max(sapply(1:detectCores(), function(x) gcd(x, chain_length)))

## Only need to compare forcings and (relevant) predicted temps, so which RCP  
## doesn't really matter.
rcp_type <- "rcp26" 

## Limit climate DF to one RCP (doesn't matter which one) and then add model-derived 
## standard errors (deviations) to simulate "true" future values.
climate <- 
  climate %>%
  filter(rcp == rcp_type) 
y_dev <- read_csv(here("Results/Evidence/y-dev.csv"))
climate$noise <- sample(y_dev$dev, nrow(climate), replace = T)

## Also add noise to future "mean" covariate values; otherwise collinearity 
## problems for future regressions b/c values stay constant.
climate$volc_noise <- rnorm(nrow(climate), sd = .1)
climate$soi_noise <- rnorm(nrow(climate), sd = .1)
climate$amo_noise <- rnorm(nrow(climate), sd = .1)

## Simulate climate data based on noninformative parameters
climate <- 
  climate %>% 
  mutate(
    volc_sim = ifelse(year <= 2005, volc_mean, volc_mean + volc_noise),
    soi_sim = ifelse(year <= 2005, soi_mean, soi_mean + soi_noise),
    amo_sim = ifelse(year <= 2005, amo_mean, amo_mean + amo_noise)
    ) %>%
  mutate(had_sim = -0.102 + 0.418*trf + 0.051*volc_sim + -0.028*soi_sim + 0.473*amo_sim + noise)


## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 

## EVIDENCE FUNCTION
## DESCRIPTION:
## The function loops along a sequence of (decreasing) sigma values for a given 
## mu. Once it reaches the end of the sequence for that mu, it moves on to the  
## next, lower mu value and repeats the procedure. It does this for two TCR
## threshold levels: 1.3 °C and 1.5°C. (Note that the function will need to be 
## adjusted  if different thresholds are chosen.) The mean posterior TCR values  
## must be at least as big as the relevant threshold to qualify as fully converged  
## with mainstream. 
## NOTE: This function, and hence the whole script, is best run from the shell
## ("~$ Rscript evidence.R") and not from within the RStudio IDE when the cluster 
## type (cl_type) is FORK. Otherwise the function will hang at random points in 
## the loop, requiring manual stopping. See ?mcfork for more details. I still 
## recommend FORK over SOCK when it is available (i.e. Linux or MacOS) because it 
## is considerably (3x) faster and less memory intensive. That all being said, the
## function also allows caching for recovering completed runs if a manual override
## or exit does occur.
evid_func <- 
  function(d){
    lapply(1:nrow(d), function(x){
    # pblapply(1:nrow(d), function(x){
      
      
      # 1. Try to load cached data, if already generated
      key <- list(x)
      df_evid <- loadCache(key)
      if (!is.null(df_evid)) {
        cat(paste0("\nLoaded cached results ", x, " of ", nrow(d),"\n"))
        return(df_evid);
        yrs <<- df_evid$yrs[nrow(df_evid)]
        yrs_j <<- df_evid$yrs_j[nrow(df_evid)]
      }
      
      # 2. If not available, generate it.
      cat(paste0("Generating new results ", x, " of ", nrow(d), "..."))
      
      m <- d[x, ]$mu
      s <- d[x, ]$sigma
      max_m <- d[x, ]$max_mu
      thresh <- d[x, ]$threshold
      
      yr_max <- 2100
      
      yrs_j <<- ifelse(!exists("yrs_j"), 65, yrs_j)

      yrs <<- 
        ifelse(
          !exists("yrs"),
          65, ## If yrs variable doesn't exist, start from 65, else...
          ifelse(
            thresh==1.5 & round(m, 1)==1.0 & round(s, 3)==0.250, #x == 122, ## i.e. first row where thresh == 1.5
            100,
            ifelse(
              max_m!=1,
              yrs - 1,
              yrs_j
              )
            )
          )

      
      tcr_mean <- 0 ## Place holder value
      
      while((round(tcr_mean, digits = 1) < thresh) & 
            (yrs < nrow(filter(climate, year <= yr_max)))) {
        
        ## NB: Use "<<-" to assign value to global workspace for next iteration
        yrs <<- yrs + 1 
        
        clim_df <-
          climate %>%
          filter(year <= min(year) + yrs)
        
        N <- nrow(clim_df)
        
        mu_beta <- m/3.71
        sigma_beta <- s/3.71
        
        ##------------------------------------------------------------------------------
        ## THE BUGS/JAGS MODEL.
        
        bugs_file <- bugs_model_func
        
        ##------------------------------------------------------------------------------
        ## PARALLEL SETUP.
        
        cl <- parallel::makeCluster(n_chains, type = cl_type, #"PSOCK"#"
                                    outfile = here("BUGSFiles/Evidence/debug.txt")) # no. of clusters (i.e. MCMC chains)
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
        
        tcr <- as.matrix(mod_samples[, "beta"], iters = F) * rf2x
        tcr_mean <- mean(tcr)
        
        rm(bugs_file, cl, data_list, mod_samples, par_inits, parameters, tcr)
        gc()
        
      } ### end of while loop
      
      ## "Restart" the iteration each time we get back to the "max mu" value for that
      ## sigma loop. (Max mu will always be at mu = 1 for historic data frame, but not 
      ## necessarily for future data frame.) Will begin at convergence year for the 
      ## previous sigma value.
      ## Again, use "<<-" to assign value to global workspace for next iteration
      if(max_m == 1){yrs_j <<- yrs - 1}
      
      df_evid <-
        data_frame(
          mu = m,
          sigma = s,
          tcr_mean = tcr_mean,
          yrs = yrs,
          yrs_j = yrs_j,
          thresh = thresh
          )
      
      cat("done\n")
      saveCache(df_evid, key=key, comment="evidence()")
      
      return(df_evid)
      on.exit(stopCluster(cl)) ## In case function crashes
      
    } ## end of function
    ) %>% ## end of pbapply
      bind_rows()
  } ## end of entire evid_func function


###################################
## Sceptic prior range (data frame)

mu <- round(seq(1, 0, length = 11), 1)
sigma <- round( seq(.25, .065, length = 11), 4)
threshold <- c(1.3, 1.5)

df <- 
  expand.grid(mu = mu, sigma = sigma, threshold = threshold) %>%
  group_by(threshold, sigma) %>% 
  mutate(max_mu = ifelse(mu == max(mu), 1, 0)) %>%
  ungroup() 

####
## Run the function over the full range of sceptic priors
## NOTE: Takes about 20 min to run on my laptop (12 core, 32 GB RAM). You may wish to skip 
## diectly to the already-exported file ("Data/Evidence/tcr-evidence.csv") to avoid this wait.
tic()
evid <- evid_func(df)
toc()
rm(yrs, yrs_j)

## Write data
write_csv(evid, here("Results/Evidence/tcr-evidence.csv"))
