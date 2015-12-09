rm(list = ls()) # Clear data

## NB: Decide whether using historical data only (and work *backwards* from most recent 
## date), or whether using simulated future data (and work *fowards* from earliest data).
recurse_type <- c("historic", "future")[2]

## Load packages used in the below loop (or further down in the code)##
library(readr) ## For reading in data files
library(dplyr) ## For manipulating and munging data frames
library(tidyr) ## For tidying data frames
library(purrr) ## For manipulating vectors and functions (complements dplyr)
library(LearnBayes) ## Mostly for simulating noninformative prior (using random multivarite normal command)
# library(arm)
library(rjags) ## For running the MCMC (Gibbs) sampler
# library(coda) ## For converting MCMC objects and diagnostics. Loads as rjags dependency.
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) # For extracting summary statistics from MCMC chain
library(ggplot2)
library(cowplot) ## For cowplot ggplot theme
library(ggthemes) ## For additional (e.g. "few") ggplot2 themes
require(RColorBrewer)
library(grid) ## To adjust legend key width and size in ggplot2 themes that don't naturally support a grid
library(gridExtra) ## Facilitates easier labelling in ggplot2
library(extrafont) ## For additional fonts in ggplot2
library(stargazer) ## For nice LaTeX tables

## NB: Check working directory. Should be main "C:/Users/.../sceptic-priors", even though
## now working with files from the "C:/Users/.../sceptic-priors/Recursive" sub-directory.

## Optional for replication
set.seed(123) 

## Decide on max year for the recursive regressions
if (recurse_type == "historic") {yr_max <- 2005} 
if (recurse_type == "future") {
  yr_max <- 2050#2005
  }

## Only need to compare forcings and (relevant) predicted temps. Which RCP doesn't 
## really matter.
rcp_type <- "rcp26" 

## Load climate data and deviation DF for simulating "true" future values
climate <- read_csv("./Data/climate.csv") %>%
  filter(rcp == rcp_type) %>%
  filter(year <= yr_max) 

y_dev <- read_csv("./Recursive/Data/y-dev.csv")

climate$noise <- sample(y_dev$dev, nrow(climate), replace = T)

## Also add noise to future "mean" covariate values; otherwise collinearity problems for some
## regressions.
climate$volc_noise <- rnorm(nrow(climate))
climate$soi_noise <- rnorm(nrow(climate))
climate$amo_noise <- rnorm(nrow(climate))

climate <- 
  climate %>% 
  mutate(volc_sim = ifelse(year <= 2005, volc_mean, volc_mean + volc_noise),
         soi_sim = ifelse(year <= 2005, soi_mean, soi_mean + soi_noise),
         amo_sim = ifelse(year <= 2005, amo_mean, amo_mean + amo_noise)
         ) %>%
  mutate(had_sim = 
           ifelse(year <= 2005, 
                  had,
                  -0.110 + 0.415*trf + 0.047*volc_sim + -0.028*soi_sim + 0.481*amo_sim + noise
                  )
         )

## Load some misc functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
chain_length <- 30000
n_chains <- 3

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 

tcr_rec <- list()
a <- 0 ## Count variable for combining TCR lists later

A <- 0 ## Count variable for (combined TCR) animation figures


## Loop over specified year intervals for recursive estimation
if (recurse_type == "historic") {recurse_seq <- seq(from = yr_max - 15, to = 1866, by = -1)} 
if (recurse_type == "future") {recurse_seq <- seq(from = 1880, to = yr_max, by = 1)}

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
# user   system  elapsed 
# 49.39   15.50  1472.46 
## Same for Future set until 2005. But then, e.g., going then from 2006 to 2050...
# user  system elapsed 
# 18.55    5.07  752.94 

tcr_rec <- dplyr::bind_rows(tcr_rec)

write_csv(tcr_rec, paste0("./Recursive/Data/tcr-rec-", recurse_type, ".csv"))


###################################
## Read data
tcr_rec <- read_csv(paste0("./Recursive/Data/tcr-rec-", recurse_type,".csv"))

tcr_df <- 
  tcr_rec %>%
  arrange(series) %>% 
  mutate(serieslab = factor(series))

## Next "rebind" the ni series data in the way that can be easily faceted by the four sceptic priors
## Basically involves duplicating the ni series for each sceptic prior in the data frame
tcr_df <-
  bind_rows(tcr_df %>% filter(series != "ni"),
            do.call("bind_rows",
                    sapply(c("denmod", "denstrong", "lukemod", "lukestrong"),
                           function(x) {
                             tcr_df %>%
                               filter(series == "ni") %>%
                               mutate(serieslab = x)
                             },
                           simplify = F
                           )
                    )
            )

tcr_df$serieslab <- 
  ifelse(tcr_df$serieslab == "lukemod",
         "(A) Moderate Lukewarmer", 
         ifelse(tcr_df$serieslab == "lukestrong",
                "(B) Strong Lukewarmer",
                ifelse(tcr_df$serieslab == "denmod",
                       "(C) Moderate Denier",
                       "(D) Strong Denier"
                       )
                )
         )

tcr_df$serieslab <- as.factor(tcr_df$serieslab)

## Reorder colour scheme slightly (to match the way ggplot facets)
prior_cols <- c(brewer.pal(12, "Paired")[c(4, 2, 6, 8#7
                                           )], "#000000")
prior_names <- c("Strong Denier", "Moderate Denier", 
                 "Strong Lukewarmer", "Moderate Lukewarmer", 
                 "Noninformative")

tcr_plot <- 
  ggplot(tcr_df, aes(x = year_to, #x = samp_size,
                     y = mean, col = series)) +
  theme_cowplot() +
  geom_line(lwd = .85) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = series), lty = 0, alpha = 0.4) +
  ylab(expression(bold(~degree~C))) + xlab("Year") + #xlab("Sample size (years)") +
  ## scale_y_continuous(limits = c(-1, 3)) +
  scale_colour_manual(values = prior_cols, 
                      breaks = c("denstrong","denmod","lukestrong","lukemod","ni"),
                      labels = prior_names) +
  scale_fill_manual(values = c(prior_cols[1:4], "gray60"), 
                    breaks = c("denstrong","denmod","lukestrong","lukemod","ni"),
                    labels = prior_names) +
  background_grid(major = "y", minor = "none", colour.major = "gray90") +
  facet_wrap(~ serieslab, ncol = 2) +
  theme(
    text = element_text(family = "Palatino Linotype"),
    axis.title.x = element_text(face="bold", size=18),
    axis.title.y = element_text(face="bold", size=18),#, angle = 0),
    axis.text  = element_text(size=17),
    legend.position = "none",
    strip.text = element_text(size = 17, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.margin = unit(2, "lines") ## Increase gap between facet panels
  ) 

if(recurse_type == "historic"){
  tcr_plot <- 
    tcr_plot + 
    scale_x_reverse() +
    xlab("Year")
  }

tcr_plot +
  ggsave(file = paste0("./Recursive/TablesFigures/rec-tcr-all-", recurse_type, ".pdf"),
         width = 8, height = 7, 
         device = cairo_pdf)