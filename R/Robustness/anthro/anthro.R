run_type <- "anthro"

## Load all packages, as well as some helper functions that will be used for plotting and tables
source(here::here("R/sceptic_funcs.R"))

## Create new "trf_less_anthro" variable since we will be breaking the anthropogenic component out on its own
climate <-
  climate %>%
  mutate(trf_less_anthro = trf - anthro)

## Decide on total length of MCMC chains (i.e. summed parallel chains of the JAGS model)
## Each individual chain will thus be chain_length/n_chains.
chain_length <- 30000

## Now set the JAGS parallel MCMC parameters
n_chains <- detectCores() ## no. of parallel chains. Use `detectCores()-1` if you are worried about CPU resources
n_adapt = 5000 ## no. of tuning or adaptation steps
burn_in = 1000 ## no. of burn-in steps

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 

# Run the nested loop (takes about 10min on my laptop)
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
    source(here("R/Robustness/anthro/jags-loop-anthro.R"), local = T) 
    
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
run_type <- "anthro"
pref <- here("TablesFigures/Untracked/Robustness/")
suff <- "-anthro"

source(here("R/sceptic_tablesfigures.R"))

########################################
### EXPORT TCR SUMMARY FOR LATER USE ###
########################################
tcr %>%
  filter(prior == "ni") %>%
  mutate(series = "anthro") %>%
  group_by(series) %>%
  summarise(
    mean = mean(tcr),
    q025 = quantile(tcr, .025),
    q975 = quantile(tcr, .975)
    ) %>%
  write_csv("Results/Robustness/tcr-anthro.csv")


###################################################
### Table 3: Posterior regressions coefficients ###
###################################################

## Add TCR summary info to coefficents table
coefs_tab <-
  bind_rows(
    coefs_tab,
    tcr %>%
      group_by(prior) %>%
      summarise(
        mean = mean(tcr),
        q025 = quantile(tcr, .025),
        q975 = quantile(tcr, .975)
        ) %>%
      mutate(coef = "tcr") %>%
      select(coef, mean, q025, q975, prior)
    )

## Convert table into nicer format for the paper. Requires some minor tinkering in
## in LaTeX afterwards (e.g. not all vars have 3 decimals), but close enough.
coefs_tab <- 
  coefs_tab %>%
  filter(coef != "sigma", coef != "phi") %>%
  mutate(
    mean = ifelse(coef=="tcr", decimals(mean, 1), decimals(mean, 3)),
    ci = ifelse(coef=="tcr",
                paste0("[", decimals(q025, 1), ", ", decimals(q975, 1), "]"),
                paste0("[", decimals(q025, 3), ", ", decimals(q975, 3), "]")
                )
    ) %>%
  select(prior, coef, mean, ci) %>%
  gather(key, value, -c(prior, coef)) %>% ## This step causes some values to drop the (zero) 3rd decimal
  mutate(
    prior = factor(prior, levels = c("ni", "lukemod", "lukestrong", "denmod", "denstrong")),
    key = factor(key, levels = c("mean", "ci"))
  ) %>%
  arrange(prior, coef) %>%
  spread(prior, value) %>% 
  mutate(coef = factor(coef, levels = c("beta_anthro", "beta", "gamma", "delta", "eta", "alpha", "tcr"))) %>%
  arrange(coef) %>% 
  mutate(coef = ifelse(key == "mean", paste(coef), "")) %>%
  select(-key) %>%
  as.matrix()

rownames(coefs_tab) <- match_coefs(coefs_tab[, "coef"])
colnames(coefs_tab) <- match_priors(colnames(coefs_tab))

coefs_tab[, 2:ncol(coefs_tab)] %>%
  stargazer(
    align = T, header = F, rownames = T,
    title = "Posterior regression results",
    label = paste0("tab:reg-coefs", suff),
    # notes.align = "l",
    # notes.append = T,
    # notes = c("Dependent variable: Global Mean Surface Temperature (GMST)."),
    out = paste0(pref, "tab-3", suff, ".tex")
    )

rm(coefs_tab)
