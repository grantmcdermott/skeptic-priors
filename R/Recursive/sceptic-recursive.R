## Choose run type
run_type <- "recursive"

## Load all packages, as well as some helper functions that will be used for plotting and tables
source(here::here("R/sceptic_funcs.R"))

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
chain_length <- 15000
## The below below tries to optimize the number of parallel MCMC chains given
## available CPUs, but balanced against the diminishing returns brought on by
## repeating the burn-in period for each parallel worker. 
n_chains <- n_chains_func(chain_length)

## NB: Decide whether using historical data only (and work *backwards* from most recent 
## date), or whether using simulated future data (and work *fowards* from earliest date).
recurse_type <- c("historic", "future")[1]

## Decide on max year for the recursive regressions
if (recurse_type == "historic") {yr_max <- 2005} 
if (recurse_type == "future") {yr_max <- 2075}

## Only need to compare forcings and (relevant) predicted temps, so which RCP doesn't 
## really matter.
rcp_type <- "rcp26" 

## Limit climate DF to one RCP (doesn't matter which one) and then add model-derived 
## standard errors (deviations) to simulate "true" future values.
climate <- 
  climate %>%
  filter(rcp == rcp_type) %>%
  filter(year <= yr_max) 

y_dev <- read_csv(here("Results/Evidence/y-dev.csv"))

climate$noise <- sample(y_dev$dev, nrow(climate), replace = T)

## Also add noise to future "mean" covariate values; otherwise collinearity 
## problems for future regressions b/c values stay constant.
climate$volc_noise <- rnorm(nrow(climate), sd = .1)
climate$soi_noise <- rnorm(nrow(climate), sd = .1)
climate$amo_noise <- rnorm(nrow(climate), sd = .1)

climate <- 
  climate %>% 
  mutate(
    volc_sim = ifelse(year <= 2005, volc_mean, volc_mean + volc_noise),
    soi_sim = ifelse(year <= 2005, soi_mean, soi_mean + soi_noise),
    amo_sim = ifelse(year <= 2005, amo_mean, amo_mean + amo_noise)
    ) %>%
  mutate(
    had_sim = 
      ifelse(year <= 2005, 
             had,
             -0.102 + 0.418*trf + 0.051*volc_sim + -0.028*soi_sim + 0.473*amo_sim + noise
             )
    )


## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 

## Loop over specified year intervals for recursive estimation
if (recurse_type == "historic") {recurse_seq <- seq(from = yr_max - 15, to = 1866, by = -1)} 
if (recurse_type == "future") {recurse_seq <- seq(from = 1880, to = yr_max, by = 1)}

## Run the recursive nested regressions.
## NOTE: This takes about 25min to complete on my system. If you don't feel like
## waiting that long (or even longer for older machines), skip to the bottom of
## the script to see where my saved results are stored.
tic()
tcr_rec <-
  ## Outer: Loop over recursive years
  # pblapply(seq_along(recurse_seq), function(n) {
  lapply(seq_along(recurse_seq), function(n) {
    
    cat(paste0("Generating new results ", n, " of ", length(recurse_seq), "..."))
    
    if(recurse_type == "historic"){
      clim_df <- 
        climate %>%
        filter(year >=  recurse_seq[n])
    }
    
    if(recurse_type == "future"){
      clim_df <- 
        climate %>%
        filter(year <=  recurse_seq[n]) 
    }
    
    yr_min <- min(clim_df$year)
    
    ## Inner: Loop over priors ##
    tcr_df <-
      lapply(1:nrow(priors_df), function(j){
        
        m <- priors_df[j, ]$mu
        s <- priors_df[j, ]$sigma
        prior_type <- priors_df[j, ]$prior_type
        convic_type <- priors_df[j, ]$convic_type 
        
        mu_beta <- m/3.71
        sigma_beta <- s/3.71
        
        if(prior_type == "ni")  {
          # source(here("R/Recursive/jags-recursive.R"), local = T) ## For vague noninformative riors using the rjags package
          source(here("R/Recursive/noninf-recursive.R"), local = T) ## For "proportional" noninformative prors using the LearnBayes package
        }
        else{
          source(here("R/Recursive/jags-recursive.R"), local = T)
        }

      })

    ## Extract only the .$value elements (drop surplus .$visibility element)
    tcr_df <-
      lapply(1:length(tcr_df), function(x) tcr_df[[x]]$value) %>%
      bind_rows()

    
    ##################################################################################
    ### Animated / recursive version of Figure 1 (TCR densities) for presentations ###
    ##################################################################################
    
    # ## Year tracker annotation label
    # year_tracker_lab <-
    #   ifelse(
    #     recurse_type == "historic",
    #     paste0("Looking back to ", tcr_df$year_to[1],"\n(Sample size = ", 2005-tcr_df$year_to[1]+1, " years)"),
    #     paste0("Looking forward to ", tcr_df$year_to[1],"\n(Sample size = ", tcr_df$year_to[1]-2005+1, " years)")
    #     )
    # 
    # fig_1 <- 
    #   tcr_plot(tcr_df) +
    #   annotate("text", x = 1.75, y = 1, label = year_tracker_lab, size = 3.5)
    # 
    # fig_1 +
    #   ggsave(
    #     file = here(paste0("TablesFigures/Untracked/Recursive/Animation/", str_to_title(recurse_type), "/rec-tcr-", 1000 + n, ".png")),
    #     width = 8, height = 4.5
    #     )
    
    
    ##################################################################
    ### Save the summary data frame (95 CI) for the recursive plot ###
    ##################################################################
    
    tcr_rec <-
      tcr_df %>% 
      group_by(prior, year_to) %>%
      summarise(
        tcr_mean = mean(tcr),
        tcr_q025 = quantile(tcr, p = 0.025),
        tcr_q975 = quantile(tcr, p = 0.975)
        )
    
    cat("done\n")
    
    return(tcr_rec)
  
  }) %>%
  bind_rows()
toc()
  
tcr_rec

## Write to disk for future use
write_csv(tcr_rec, here(paste0("Results/Recursive/tcr-rec-", recurse_type, ".csv")))

