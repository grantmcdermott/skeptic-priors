rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("Data/climate.csv")

## The 95% measurement error bounds are not quite symmetrical, but very close. 
# We can also compare the measurement error with model uncertainty i.e. sigma 
# from noninformative model. Note: dividing 95% ME bounds by 2 to get 1 std dev 
# (i.e. fair comparison to model sigma).
climate %>% 
  filter(rcp == "rcp26") %>%
  filter(!is.na(had_full)) %>%
  mutate(
    omega_low = (had_full - had_025)/2,
    omega_up = (had_975 - had_full)/2
    ) %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = omega_low), col = "blue") +
  geom_line(aes(y = omega_up), col = "red") +
  geom_hline(yintercept = 0.075) + ## Sigma mean = 0.075, sigma s.d. = 0.0045
  geom_hline(yintercept = 0.075 + 0.0045*2, lty = 2) + 
  geom_hline(yintercept = 0.075 - 0.0045*2, lty = 2) +
  labs(y = expression(paste("Measurement Error (", omega, ") vs Sigma (", sigma, ")")))

## Decide on total length of MCMC chains (i.e. summed parallel chains JAGS model)
## Each individual chain will thus be chain_length/n_chains.
chain_length <- 30000
## Function below ensures that the individual chains sum exactly to the desired
## total chain length, whilst still making full use of the available CPUs for
## for parallel processing power. (Note: If you want to use less than your full
## CPU allotment, use e.g. "...sapply(1:(detectCores-1)), ...)". The extra 
## parentheses is important.)
n_chains <- max(sapply(1:detectCores(), function(x) gcd(x, chain_length)))

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 

## Priors data frame
priors_df <- 
  data_frame(mu = c(0, 1, 1, 0, 0),
             sigma = c(100, 0.25, 0.065, 0.25, 0.065),
             prior_type = c("ni", "luke", "luke", "den", "den"),
             convic_type = c("", "mod", "strong", "mod", "strong")
             )
priors_df

# Run the nested loop (takes about two and a half minutes on my laptop)
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
    source("Robustness/MeasError/jags-loop-me.R", local = T)
    
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
run_type <- "me"
pref <- "TablesFigures/Untracked/Robustness/"
suff <- "-me"

source("sceptic_tablesfigures.R")

tcr %>%
  filter(prior == "ni") %>%
  mutate(series = "me") %>%
  group_by(series) %>%
  summarise(
    mean = mean(tcr),
    q025 = quantile(tcr, .025),
    q975 = quantile(tcr, .975)
    ) %>%
  write_csv("Data/Robustness/tcr-me.csv")
