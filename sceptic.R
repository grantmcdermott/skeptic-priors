rm(list = ls()) # Clear data

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

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## Load some misc functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
chain_length <- 30000
n_chains <- 3

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 


## Loop over prior ##
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
          label = "tab:coefs-all",
          #             notes.align = "l",
          #             notes.append = T,
          #             notes = c("Dependent variable: Global Mean Surface Temperature (GMST).")
          out = "./TablesFigures/coefs-all.tex"
          )


## Combined TCR density plots with shaded IPCC "likely" range ##

## First, mostly for presentations, priors only ###
pdf(file = "./TablesFigures/tcr-combined-prior.pdf", 
    width = 16, height = 10, family = "Palatino")
par(las = 1, 
    mar = c(5.5, 5.7, 4, 2) + 0.1, 
    mgp = c(4, 1.5, 0)) # Change margin and tick labels to give space for different font
plot(1, type = "n", 
     ylab = "Density", xlab = expression(~degree~C), 
     main = "", cex.lab = 2.5, cex.axis = 2.25, cex.main = 2.25, cex.sub = 2.25,
     xlim = c(-1, 3),
     ylim = c(0, 8),
     font.lab = 2)
## IPCC "likely" range
x <- c(1.0, 2.5, 2.5, 1.0)
y <- c(8, 8, 0, 0) # y <- c(1/1.4, 1/1.4, 0, 0) ## If uniform over 1.0-2.5 deg C
polygon(x, y, col = "gray90", border = NA)
## Moderate lukewarmer
curve(dnorm(x, mean = 1, sd = 0.25), 
      from=0, to = 2, lwd = 4, lty = 2, col = prior_cols[4], add = T, yaxt = "n")
## Strong lukewarmer
curve(dnorm(x, mean = 1, sd = 0.065), 
      from = 0.5, to = 1.5, lwd = 4, lty = 2, col = prior_cols[3], add = T, yaxt = "n")
## Moderate Denier
curve(dnorm(x, mean = 0, sd = 0.25), 
      from = -1, to = 1, lwd = 4, lty = 2, col = prior_cols[2], add = T, yaxt = "n")
## Strong denier
curve(dnorm(x, mean = 0, sd = 0.065), 
      from = -0.5, to = 0.5, lwd = 4, lty = 2, col = prior_cols[1], add = T, yaxt = "n")
## Legend
legend("topleft", 
       prior_names[c(1:4)], 
       col = prior_cols[c(1:4)], 
       cex = 2, lty = 1, lwd = 4, bty = "n")
# legend("topright", 
#        c("IPCC ''likely''"
#        cex = 1.25, lty = 1, col = "gray90", lwd = 10, bty = "n")
dev.off()

## Now with both priors and posterior TCR densities ##
pdf(file = "./TablesFigures/tcr-combined.pdf", width = 16, height = 10, family = "Palatino")
par(las = 1, 
    mar = c(5.5, 5.7, 4, 2) + 0.1, 
    mgp = c(4, 1.5, 0)) # Change margin and tick labels to give space for different font
plot(1, type = "n", 
     ylab = "Density", xlab = expression(~degree~C), 
     main = "", cex.lab = 2.5, cex.axis = 2.25, cex.main = 2.25, cex.sub = 2.25,
     xlim = c(-1, 3),
     ylim = c(0, 8),
     font.lab = 2)
## IPCC "likely" range
x <- c(1.0, 2.5, 2.5, 1.0)
y <- c(8, 8, 0, 0) # y <- c(1/1.4, 1/1.4, 0, 0) ## If uniform over 1.0-2.5 deg C
polygon(x, y, col = "gray90", border = NA)
## Noninformative
lines(density(tcr_ni), lwd = 4, col = prior_cols[5])
## Moderate lukewarmer
curve(dnorm(x, mean = 1, sd = 0.25), 
      from = 0, to = 2, lwd = 4, lty = 2, col = prior_cols[4], add = T, yaxt = "n")
lines(density(tcr_luke_mod), lwd = 4, col = prior_cols[4])
## Strong lukewarmer
curve(dnorm(x, mean = 1, sd = 0.065), 
      from = 0.5, to = 1.5, lwd = 4, lty = 2, col = prior_cols[3], add = T, yaxt = "n")
lines(density(tcr_luke_strong), lwd = 4, col = prior_cols[3])
## Moderate Denier
curve(dnorm(x, mean = 0, sd = 0.25), 
      from = -1, to = 1, lwd = 4, lty = 2, col = prior_cols[2], add = T, yaxt = "n")
lines(density(tcr_den_mod), lwd = 4, col = prior_cols[2])
## Strong denier
curve(dnorm(x, mean = 0, sd = 0.065),
      from = -0.5, to = 0.5, lwd = 4, lty = 2, col = prior_cols[1], add = T, yaxt = "n")
lines(density(tcr_den_strong), lwd = 4, col = prior_cols[1])
## Legend
legend("topleft", 
       prior_names, 
       col = prior_cols,
       cex = 2, lty = 1, lwd = 4, bty = "n")
# legend("topright", 
#        "IPCC ''likely''", 
#        cex = 1.25, lty = 1, col = "gray90"), lwd = 10, bty = "n")
dev.off()


## Combined temperature prediction in 2100 by scenario ##

for (i in 1:4)  {  
  combined_rcp <- c("rcp26", "rcp45", "rcp60", "rcp85")[i] 
  pdf(file = paste0("./TablesFigures/all", combined_rcp, ".pdf"), 
      width = 10, height = 8, family = "Palatino")
  par(las = 1, mar = c(5.5, 5, 4, 2) + 0.1) # Change margin and tick labels to give space for different font
  plot(density(filter(ni_2100, rcp == combined_rcp)$temp),    
       xlab = expression(~degree~C), xlim = c(0, 4), ylim = c(0, 5), lwd = 3.5, #xlim was c(-1,1.5)
       font.lab = 2,
       # main = rcp_names[i],
      main = "",
      cex.lab = 3, cex.axis = 2.5, cex.main = 2.5, cex.sub = 2.25)
  lines(density(filter(den_strong_2100, rcp == combined_rcp)$temp), col = prior_cols[1], lty = 1, lwd = 3.5)
  lines(density(filter(den_mod_2100, rcp == combined_rcp)$temp), col = prior_cols[2], lty = 1, lwd = 3.5)
  lines(density(filter(luke_strong_2100, rcp == combined_rcp)$temp), col = prior_cols[3], lty = 1, lwd = 3.5)
  lines(density(filter(luke_mod_2100, rcp == combined_rcp)$temp), col = prior_cols[4], lty = 1, lwd = 3.5)
#   legend("topright", legend = prior_names, 
#          lty = 1:1,
#          lwd = 3.5:3.5,
#          col = prior_cols,
#          bty = "n",
#          cex = 2)
  dev.off()
}

pdf(file = "./TablesFigures/all-combined-legend.pdf", 
    width = 10, height = .8, family = "Palatino")
par(mar = c(0, 0, 0, 0))
plot.new()
legend(x = "top",inset = 0, ncol = 3,
       legend = prior_names, 
       col = prior_cols, 
       bty = "n", 
       lwd = 3, 
       seg.len = 2, 
       cex = 1.5#, 
       # horiz = TRUE
       )
dev.off()


######################################################
### 95% Conf. Intervals combined in single table ###
######################################################

## Use own function to help pull desired data in summary form for table ##
clean_func <- function(x) {
  rbind(decimals(mean(x), 1),
        paste0("[",
               decimals(quantile(x, p = 0.025, names = F), 1), 
               ", ", 
               decimals(quantile(x, p = 0.975, names = F), 1),
               "]"))
  }

## TCRs ##
# qnorm(c(.025, .975), mean = 0, sd = 100)
# [1] -191.8876  195.9964
tab_tcr_all <- 
  rbind(cbind(rbind("0", "[-196.0  196.0]"), clean_func(tcr_ni)),
        cbind(rbind("1.0","[0.5, 1.5]"), clean_func(tcr_luke_mod)),
        cbind(rbind("1.0","[0.9, 1.1]"), clean_func(tcr_luke_strong)),
        cbind(rbind("0.0","[-0.5, 0.5]"), clean_func(tcr_den_mod)),
        cbind(rbind("0.0","[-0.1, 0.1]"), clean_func(tcr_den_strong))
        )

rm(tcr_ni, tcr_luke_mod, tcr_luke_strong, tcr_den_mod, tcr_den_strong)

colnames(tab_tcr_all) <- c("Prior", "Posterior")
rownames(tab_tcr_all) <- c("Noninformative", "", 
                           "Moderate Lukewarmer", "", 
                           "Strong Lukewarmer", "", 
                           "Moderate Denier", "", 
                           "Strong Denier", "")

tab_tcr_all

stargazer(tab_tcr_all,
          title ="Transient climate response (TCR), $^\\circ$C",
          header = F,
          label = "tab:tcr-all",
          # align = T, 
          notes.align = "l",
          notes = c(#"\\footnotesize TCR is the temperature change that results from a steady doubling of CO$^2$.",
                    "\\footnotesize Mean estimates are given, with 95\\% probability interval in brackets."),
          out = "./TablesFigures/tcr-all.tex"
          )

## Temperatures at 2100 ##
tab_2100_all <- c("Noninformative", "",
                  "Moderate Lukewarmer", "", 
                  "Strong Lukewarmer", "", 
                  "Moderate Denier", "", 
                  "Strong Denier", "")

for (i in 1:4)  {
  combined_rcp <- c("rcp26", "rcp45", "rcp60", "rcp85")[i]
  tab_2100 <- NULL
  tab_2100 <- rbind(tab_2100,
                    clean_func(filter(ni_2100, rcp == combined_rcp)$temp),
                    clean_func(filter(luke_mod_2100, rcp == combined_rcp)$temp),
                    clean_func(filter(luke_strong_2100, rcp == combined_rcp)$temp),
                    clean_func(filter(den_mod_2100, rcp == combined_rcp)$temp),
                    clean_func(filter(den_strong_2100, rcp == combined_rcp)$temp)
                    )
  tab_2100_all <- cbind(tab_2100_all, tab_2100)
}

rm(combined_rcp, tab_2100, 
   ni_2100, luke_mod_2100, luke_strong_2100, den_mod_2100, den_strong_2100)

colnames(tab_2100_all) <- c("Prior", "RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5")

tab_2100_all

stargazer(tab_2100_all,
          title = "Temperature anomaly by 2100: Prediction using various priors",
          label = "tab:predictions",
          header = F,
          align = T, 
          notes.align = "l",
          notes = c("\\footnotesize Relative to pre-industrial (1851-1880) average.",
                    "\\footnotesize Mean estimates are given, with 95\\% credible interval in brackets."),
          out = "./TablesFigures/predictions2100.tex"
          )