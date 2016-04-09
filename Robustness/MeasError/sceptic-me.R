rm(list = ls()) # Clear data

## Load packages used in the below loop (or further down in the code)##
library(readr) ## For reading in data files
library(dplyr) ## For manipulating and munging data frames
library(tidyr) ## For tidying data frames
library(purrr) ## For manipulating vectors and functions (complements dplyr)
# library(LearnBayes) ## Mostly for simulating noninformative prior (using random multivarite normal command)
library(rjags) ## For running the MCMC (Gibbs) sampler
# library(coda) ## For converting MCMC objects and diagnostics. Loads as rjags dependency.
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) # For extracting summary statistics from MCMC chain
library(ggplot2)
library(cowplot) ## For cowplot ggplot theme
# library(ggthemes) ## For additional (e.g. "few") ggplot2 themes
library(RColorBrewer)
# library(grid) ## To adjust legend key width and size in ggplot2 themes that don't naturally support a grid
# library(gridExtra) ## Facilitates easier labelling in ggplot2
library(extrafont) ## For additional fonts in ggplot2
library(stargazer) ## For nice LaTeX tables

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## Load some misc functions that will be used for plotting and tables
source("sceptic_funcs.R")
## Choose font type for graphs (note extrafont package installation instructions)
font_type <- c("Palatino Linotype", "Lato")[1]

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
## Total chain length will thus be chain_length * n_chains
chain_length <- 3000#10000
n_chains <- 3

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 


## Loop over prior ##
for (k in 1:3)  {  
  prior_type <- c("ni", "luke", "den")[k] 
  
  ## Loop over conviction strength ##
  if(prior_type == "ni")  {
    convic_type <- ""
    source("./Robustness/MeasError/jags-loop.R") ## For vague noninformative riors using the rjags package
    # source("./Robustness/MeasError/noninf-loop.R") ## For "proportional" noninformative prors using the LearnBayes package
  }
  else{for (j in 1:2)  {   
    convic_type <- c("mod", "strong")[j]
    source("./Robustness/MeasError/jags-loop.R")   
  } } ## End of conviction loop
  
} ## End of prior loop


##################################
### COMBINED TABLES AND GRAPHS ###
##################################

## Remove data unnecessary to further analysis ##
rm(climate, n_chains, chain_length, prior_type, convic_type,
   coefs_tab, predictions, temp_2100, rf2x, tcr)


## Set colours and names for consistency among graphs ##
# prior_cols <- c("dodgerblue2", "limegreen", "orange", "red2", "gray20")
prior_cols <- c(brewer.pal(12, "Paired")[c(2, 4, 8, 6)], "#000000")
prior_names <- c("Strong Denier", "Moderate Denier", 
                 "Strong Lukewarmer", "Moderate Lukewarmer", 
                 "Noninformative")


## Combined regression coefficient table 
### Clunky, but I don't feel like changing earlier code
coefs_mean_tab_all <- 
  cbind(coefs_ni[, 1], 
        coefs_luke_mod[, 1], 
        coefs_luke_strong[, 1],
        coefs_den_mod[, 1], coefs_den_strong[, 1]
        )
coefs_95_tab_all <- 
  cbind(coefs_ni[, 2], 
        coefs_luke_mod[, 2], 
        coefs_luke_strong[, 2],
        coefs_den_mod[, 2], 
        coefs_den_strong[, 2]
        )

coefs_tab_all <- rbind(coefs_mean_tab_all[1, ], coefs_95_tab_all[1, ])
for(i in 2:5){
  coefs_tab_all <- rbind(coefs_tab_all, coefs_mean_tab_all[i, ], coefs_95_tab_all[i, ])
}

rm(coefs_ni, coefs_luke_mod, coefs_luke_strong, coefs_den_mod, coefs_den_strong,
   coefs_mean_tab_all, coefs_95_tab_all)

colnames(coefs_tab_all) <- rev(prior_names)
rownames(coefs_tab_all) <- c("Total radiative forcing", "", 
                             "Stratospheric aerosols", "",
                             "SOI", "",
                             "AMO", "",
                             "Constant", "")

stargazer(coefs_tab_all, align = T, header = F, 
          title = "Posterior regression results",
          label = "tab:coefs-all-me",
          #             notes.align = "l",
          #             notes.append = T,
          #             notes = c("Dependent variable: Global Mean Surface Temperature (GMST).")
          out = "./Robustness/TablesFigures/coefs-all-me.tex"
)