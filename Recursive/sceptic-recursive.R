rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
chain_length <- 5000
n_chains <- detectCores() - 1 

## NB: Decide whether using historical data only (and work *backwards* from most recent 
## date), or whether using simulated future data (and work *fowards* from earliest data).
recurse_type <- c("historic", "future")[1]

## Decide on max year for the recursive regressions
if (recurse_type == "historic") {yr_max <- 2005} 
if (recurse_type == "future") {yr_max <- 2075}

## Only need to compare forcings and (relevant) predicted temps, so which RCP doesn't 
## really matter.
rcp_type <- "rcp26" 

## Load climate data and deviation DF for simulating "true" future values
climate <- read_csv("./Data/climate.csv") %>%
  filter(rcp == rcp_type) %>%
  filter(year <= yr_max) 

y_dev <- read_csv("./Evidence/Data/y-dev.csv")

climate$noise <- sample(y_dev$dev, nrow(climate), replace = T)

## Also add noise to future "mean" covariate values; otherwise collinearity 
## problems for future regressions b/c values stay constant.
climate$volc_noise <- rnorm(nrow(climate), sd = .1)
climate$soi_noise <- rnorm(nrow(climate), sd = .1)
climate$amo_noise <- rnorm(nrow(climate), sd = .1)

climate <- 
  climate %>% 
  mutate(volc_sim = ifelse(year <= 2005, volc_mean, volc_mean + volc_noise),
         soi_sim = ifelse(year <= 2005, soi_mean, soi_mean + soi_noise),
         amo_sim = ifelse(year <= 2005, amo_mean, amo_mean + amo_noise)
         ) %>%
  mutate(had_sim = 
           ifelse(year <= 2005, 
                  had,
                  -0.102 + 0.418*trf + 0.051*volc_sim + -0.028*soi_sim + 0.473*amo_sim + noise
                  )
         )


## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 

tcr_rec <- list()
a <- 0 ## Count variable for combining TCR lists later
A <- 0 ## Count variable for (combined TCR) animation figures


## Loop over specified year intervals for recursive estimation
if (recurse_type == "historic") {recurse_seq <- seq(from = yr_max - 15, to = 1866, by = -1)} 
if (recurse_type == "future") {recurse_seq <- seq(from = 1880, to = yr_max, by = 1)}

## Run the recursive regression for each prior type
## NOTE: This takes about 25min to complete on my system. If you don't feel like
## waiting that long (or even longer for older machines) then skip to line 105
## and read in the saved results from a previous session.
ptm <- proc.time()
for (n in recurse_seq) {
  ## Loop over prior ##
  for (k in 1:3)  {  
    prior_type <- c("ni", "luke", "den")[k] 
    ## Loop over conviction strength ##
    if(prior_type == "ni")  {
      convic_type <- ""
      # source("./Recursive/jags-recursive.R") ## For vague noninformative riors using the rjags package
      source("./Recursive/noninf-recursive.R") ## For "proportional" noninformative prors using the LearnBayes package
    }
    else{for (j in 1:2)  {   
      convic_type <- c("mod", "strong")[j]
      source("./Recursive/jags-recursive.R")   
    } } ## End of conviction loop
  } ## End of prior loop
} ## End of recursive loop
proc.time() - ptm
## Historic
#   user   system  elapsed 
# 42.004   28.320 1438.644
## Same for Future set until 2005. But then, e.g., going then from 2006 to 2050...
# user  system elapsed 
# 18.55    5.07  752.94 

tcr_rec <- dplyr::bind_rows(tcr_rec)

write_csv(tcr_rec, paste0("./Recursive/Data/tcr-rec-", recurse_type, ".csv"))


###################################
## Read data
tcr_rec <- read_csv(paste0("./Recursive/Data/tcr-rec-", recurse_type,".csv"))

tcr_rec <- 
  tcr_rec %>%
  arrange(series) %>% 
  mutate(serieslab = series)

## Next "rebind" the ni series data in the way that can be easily faceted by the four sceptic priors
## Basically involves duplicating the ni series for each sceptic prior in the data frame
tcr_rec <-
  bind_rows(tcr_rec %>% filter(series != "ni"),
            do.call("bind_rows",
                    lapply(c("lukemod", "lukestrong", "denmod", "denstrong"),
                           function(x) {
                             tcr_rec %>%
                               filter(series == "ni") %>%
                               mutate(serieslab = x)
                             }
                           )))

tcr_rec$serieslab <- 
  ifelse(tcr_rec$serieslab == "lukemod",
         "(A) Moderate Lukewarmer", 
         ifelse(tcr_rec$serieslab == "lukestrong",
                "(B) Strong Lukewarmer",
                ifelse(tcr_rec$serieslab == "denmod",
                       "(C) Moderate Denier",
                       "(D) Strong Denier"
                       )))

tcr_plot <- 
  ggplot(tcr_rec, aes(x = year_to, #x = samp_size,
                     y = mean, col = series)) +
  geom_line(lwd = .85) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = series), lty = 0, alpha = 0.4) +
  ylab(expression(~degree~C)) + xlab("Year") + #xlab("Sample size (years)") +
  ## scale_y_continuous(limits = c(-1, 3)) +
  scale_colour_manual(values = prior_cols, 
                      limits = c("denstrong","denmod","lukestrong","lukemod","ni"),
                      labels = prior_names) +
  scale_fill_manual(values = c(prior_cols[1:4], "gray60"), 
                    limits = c("denstrong","denmod","lukestrong","lukemod","ni"),
                    labels = prior_names) +
  background_grid(major = "y", minor = "none", colour.major = "gray90") +
  facet_wrap(~ serieslab, ncol = 2) +
  theme_recursive 

if(recurse_type == "historic"){
  tcr_plot <- 
    tcr_plot + 
    scale_x_reverse(breaks = seq(max(tcr_rec$year_to), min(tcr_rec$year_to), by = -30)) +
    xlab("Year")
  }

tcr_plot +
  ggsave(file = paste0("./Recursive/TablesFigures/rec-tcr-all-", recurse_type, ".pdf"),
         width = 8, height = 7, 
         device = cairo_pdf)


## As above, but using stippled lines for 95% CI instead of geom_ribbon
tcr_plot_lines <- 
  ggplot(tcr_rec, aes(x = year_to, col = series)) +
  geom_line(aes(y = mean), lwd = .85) +
  geom_line(aes(y = q025), lty = 2) +
  geom_line(aes(y = q975), lty = 2) +
  ylab(expression(~degree~C)) + xlab("Year") + #xlab("Sample size (years)") +
  ## scale_y_continuous(limits = c(-1, 3)) +
  scale_colour_manual(values = prior_cols, 
                      limits = c("denstrong","denmod","lukestrong","lukemod","ni"),
                      labels = prior_names) +
  scale_fill_manual(values = c(prior_cols[1:4], "gray60"), 
                    limits = c("denstrong","denmod","lukestrong","lukemod","ni"),
                    labels = prior_names) +
  background_grid(major = "y", minor = "none", colour.major = "gray90") +
  facet_wrap(~ serieslab, ncol = 2) +
  theme_recursive 

if(recurse_type == "historic"){
  tcr_plot_lines <- 
    tcr_plot_lines + 
    scale_x_reverse(breaks = seq(max(tcr_rec$year_to), min(tcr_rec$year_to), by = -30)) +
    xlab("Year")
}

tcr_plot_lines +
  ggsave(file = paste0("./Recursive/TablesFigures/rec-tcr-all-", recurse_type, "-lines.pdf"),
         width = 8, height = 7, 
         device = cairo_pdf)