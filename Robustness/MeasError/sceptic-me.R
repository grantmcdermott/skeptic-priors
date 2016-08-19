rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## The 95% measurement error bounds are not quite symmetrical, but very close. 
# We can also compare the measurement error with model uncertainty i.e. sigma 
# from noninformative model. Note: dividing 95% ME bounds by 2 to get 1 std dev 
# (i.e. fair comparison to model sigma).
ggplot(climate %>% 
         filter(rcp == "rcp26") %>%
         filter(!is.na(had_full)) %>%
         mutate(omega_low = (had_full - had_025)/2,
                omega_up = (had_975 - had_full)/2), 
       aes(x = year)) +
  geom_line(aes(y = omega_low), col = "blue") +
  geom_line(aes(y = omega_up), col = "red") +
  geom_hline(yintercept = 0.075) + ## Sigma mean = 0.075, sigma s.d. = 0.0045
  geom_hline(yintercept = 0.075 + 0.0045*2, lty = 2) + 
  geom_hline(yintercept = 0.075 - 0.0045*2, lty = 2) +
  labs(y = expression(paste("Measurement Error (", omega, ") vs Sigma (", sigma, ")")))

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
## Total chain length will thus be chain_length * n_chains
chain_length <- 5000
n_chains <- detectCores() - 1 

## Preallocate coefficients, tcr and temperature in 2100 lists for loop
coefs_tab <- list()
l <- 0 ## count variable for coefs_tab list
tcr <- list()
all_2100 <- list()

ptm <- proc.time()
## Loop over prior ##
for (k in 1:3)  {  
  prior_type <- c("ni", "luke", "den")[k] 
  
  ## Loop over conviction strength ##
  if(prior_type == "ni")  {
    convic_type <- ""
    source("./Robustness/MeasError/jags-loop-me.R") ## For vague noninformative riors using the rjags package
    # source("./Robustness/MeasError/noninf-loop-me.R") ## For "proportional" noninformative prors using the LearnBayes package
  }
  else{for (j in 1:2)  {   
    convic_type <- c("mod", "strong")[j]
    source("./Robustness/MeasError/jags-loop-me.R")   
  } } ## End of conviction loop
  
} ## End of prior loop
ptm <- proc.time() - ptm
# user  system elapsed 
# 19.44    1.95  122.73 

tcr %>%
  bind_rows() %>%
  
  
  mutate(series = "me") %>%
  summarise(mean = round(mean(tcr), 2),
            q025 = round(quantile(tcr, .025), 2),
            q975 = round(quantile(tcr, .975), 2)) %>%
  arrange(mean)
tcr_secondary_tab


##################################
### COMBINED TABLES AND GRAPHS ###
##################################
pref <- "./Robustness/TablesFigures/"
suff <- "-me"

source("sceptic_tablesfigures.R")

tcr %>%
  filter(prior == "ni") %>%
  mutate(series = "me") %>%
  group_by(series) %>%
  summarise(mean = mean(tcr),
            q025 = quantile(tcr, .025),
            q975 = quantile(tcr, .975)) %>%
  write_csv("Robustness/Data/tcr-me.csv")
