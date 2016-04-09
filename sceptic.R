rm(list = ls()) # Clear data

## Load packages used in the below loop (or further down in the code)##
library(readr) ## For reading in data files
library(LearnBayes) ## Mostly for simulating noninformative prior (using random multivarite normal command)
# library(arm)
library(rjags) ## For running the MCMC (Gibbs) sampler
# library(coda) ## For converting MCMC objects and diagnostics. Loads as rjags dependency.
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) ## devtools::install_github("johnbaums/jagstools") For extracting summary statistics from MCMC chain
library(ggplot2)
library(cowplot) ## For cowplot ggplot theme
library(ggthemes) ## For additional ggplot2 themes (e.g. "few") 
library(RColorBrewer)
library(grid) ## To adjust legend key width and size in ggplot2 themes that don't naturally support a grid
library(gridExtra) ## Facilitates easier labelling in ggplot2 (e.g. annote with extrafont fonts)
library(extrafont) ## For additional fonts in ggplot2
library(stargazer) ## For nice LaTeX tables
library(dplyr) ## For manipulating and munging data frames
library(tidyr) ## For tidying data frames
library(purrr) ## For manipulating vectors and functions (complements dplyr)

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## Choose font type for graphs (note extrafont package installation instructions)
font_type <- c("Palatino Linotype", "Lato", "Arial")[1]
## Assign colours and names for later graphs ##
rcp_names <- c(expression("RCP 2.6 (420 ppmv CO"[2]*")"),
               expression("RCP 4.5 (540 ppmv CO"[2]*")"),
               expression("RCP 6.0 (670 ppmv CO"[2]*")"),
               expression("RCP 8.5 (940 ppmv CO"[2]*")"))  
# rcp_cols <- c("limegreen", "orchid", "orange", "red2")
rcp_cols <- c("darkgreen", "darkorchid", "darkorange2", "darkred")
rcp_fills <- c("lightgreen", "orchid", "orange", "red")

## Load some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
## Total chain length will thus be chain_length * n_chains
chain_length <- 10000
n_chains <- detectCores() - 1 


## Preallocate coefficients, tcr and temperature in 2100 lists for loop
coefs_tab <- list()
tcr <- list()
all_2100 <- list()
l <- 0 ## count variable for coefs_tab list


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
# 38.19    5.25  133.73 


##################################
### COMBINED TABLES AND GRAPHS ###
##################################

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 

## Remove data unnecessary to further analysis ##
rm(climate, n_chains, chain_length, prior_type, convic_type,
   theme_coefs, theme_pred)


## Set colours and names for consistency among graphs ##
# prior_cols <- c("dodgerblue2", "limegreen", "orange", "red2", "gray20")
prior_cols <- c(brewer.pal(12, "Paired")[c(2, 4, 8, 6)], "#000000")
prior_names <- c("Strong Denier", "Moderate Denier", 
                 "Strong Lukewarmer", "Moderate Lukewarmer", 
                 "Noninformative")


###################
### COEFICIENTS ###
###################

## Unlist coefficents list
coefs_tab <- bind_rows(coefs_tab) 
## Convert into nice table format for the paper
coefs_tab <- 
  coefs_tab %>%
  filter(coef != "sigma") %>%
  mutate(mean = decimals(mean, 3),
         ci = paste0("[", decimals(q025, 3), ", ", decimals(q975, 3), "]")
         ) %>%
  select(prior, coef, mean, ci) %>%
  gather(key, value, -c(prior, coef)) %>%
  mutate(prior = factor(prior, levels = c("ni", "lukemod", 
                                          "lukestrong", "denmod", "denstrong")),
         key = factor(key, levels = c("mean", "ci"))) %>%
  arrange(prior, coef) %>%
  spread(prior, value) %>%
  mutate(coef = ifelse(key == "mean", paste(coef), "")) %>%
  select(-key) %>%
  as.matrix()
rownames(coefs_tab) <- match_coefs(coefs_tab[, "coef"])
colnames(coefs_tab) <- match_priors(colnames(coefs_tab))

coefs_tab[,2:ncol(coefs_tab)] %>%
  stargazer(align = T, header = F, rownames = T,
            title = "Posterior regression results",
            label = "tab:coefs-all",
            #             notes.align = "l",
            #             notes.append = T,
            #             notes = c("Dependent variable: Global Mean Surface Temperature (GMST).")
            out = "./TablesFigures/coefs-all.tex"
            )
  

###################
###     TCR     ###
###################

## Unlist tcr list and generate values
tcr <-
  bind_rows(
          mapply(
            function(x) {
              tcr[[x]]$tcr <- tcr[[x]]$beta * rf2x
              return(tcr[[x]])
            },
            seq(1:5), 
            SIMPLIFY = F
            )
          ) 

## TCR Density plot
tcr %>%
  mutate(prior = factor(match_priors(prior),
                        levels = prior_names)) %>%
  ggplot(aes(x = tcr, col = prior)) +
  geom_line(stat = "density") +
  labs(x = expression(~degree*C), y = "Density") +
  xlim(-1, 3) + 
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf,
           alpha = .2) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = .065), 
                lty=2, col=prior_cols[1]) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = .25), 
                lty=2, col=prior_cols[2]) + 
  stat_function(fun = dnorm, args = list(mean = 1, sd = .065), 
                lty=2, col=prior_cols[3]) +
  stat_function(fun = dnorm, args = list(mean = 1, sd = .25), 
                lty=2, col=prior_cols[4]) + 
  scale_colour_manual(values = prior_cols) +
  guides(col = guide_legend(nrow = 2)) +
  theme(
    axis.line.x = element_line(linetype = 1), ## Temporary bug(?) in cowplot theme: missing axis line
    axis.line.y = element_line(linetype = 1), ## Ditto
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank()
    ) +
  ggsave(file = "./TablesFigures/tcr-combined.pdf",
         width = 6, height = 4,
         device = cairo_pdf)

## Just the priors this time (plus ni posterior for comparison) for presentations
tcr %>%
  filter(prior == "ni") %>%
  mutate(prior = factor(match_priors(prior),
                        levels = prior_names)) %>%
  ggplot(aes(x = tcr, col = prior)) +
  geom_line(stat = "density") +
  labs(x = expression(~degree*C), y = "Density") +
  xlim(-1, 3) +
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf,
           alpha = .2) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = .065),
                lty=2, aes(col=prior_names[1])) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = .25),
                lty=2, aes(col=prior_names[2])) +
  stat_function(fun = dnorm, args = list(mean = 1, sd = .065),
                lty=2, aes(col=prior_names[3])) +
  stat_function(fun = dnorm, args = list(mean = 1, sd = .25),
                lty=2, aes(col=prior_names[4])) + 
  scale_colour_manual(values = prior_cols,
                      limits = prior_names) +
  guides(col = guide_legend(nrow = 2)) +
  theme(
    axis.line.x = element_line(linetype = 1), ## Temporary bug(?) in cowplot theme: missing axis line
    axis.line.y = element_line(linetype = 1), ## Ditto
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank()
    ) +
  ggsave(file = "./TablesFigures/tcr-combined-prior.pdf",
         width = 6, height = 4,
         device = cairo_pdf)

## Table summary of the above ##
tcr_tab <-
  tcr %>%
  group_by(prior) %>%
  summarise(mean = mean(tcr),
            q025 = quantile(tcr, .025),
            q975 = quantile(tcr, .975)) 
## Format for nice looking table in LaTeX
tcr_tab <-
  tcr_tab %>%
  mutate(mean = round(mean, 1),
         ci = paste0("[", round(q025, 1), ", ", round(q975, 1), "]")
         ) %>%
  select(prior, mean, ci) %>%
  gather(key, Posterior, -prior) %>%
  mutate(prior = factor(prior, levels = c("ni", "lukemod", 
                                          "lukestrong", "denmod", "denstrong")),
         key = factor(key, levels = c("mean", "ci"))) %>%
  arrange(prior) %>%
  mutate(prior = ifelse(key == "mean", paste(prior), "")) %>%
  select(-key) %>%
  rename(Prior = prior) %>%
  as.matrix()
rownames(tcr_tab) <- match_priors(tcr_tab[, "Prior"])
tcr_tab[, "Prior"] <- c("-", "",
                        "1.0","[0.5, 1.5]",
                        "1.0","[0.9, 1.1]",
                        "0.0","[-0.5, 0.5]",
                        "0.0","[-0.1, 0.1]")

stargazer(tcr_tab,
          title ="Transient climate response (TCR), $^\\circ$C",
          header = F,
          label = "tab:tcr-tab",
          # align = T,
          notes.align = "l",
          notes = c(paste0("\\footnotesize Mean estimates are given, with square ",
                           "brackets denoting 95\\% probability interval.")),
          out = "./TablesFigures/tcr-tab.tex"
          )



####################
### TEMP IN 2100 ###
####################

## Combined temperature prediction in 2100 by scenario ###
## Unlist temperature in 2100 list
all_2100 <-
  bind_rows(all_2100) 

## All 2100 density plot
all_2100 %>%
  mutate(prior = factor(match_priors(prior),
                        levels = prior_names)) %>%
  mutate(rcp = match_rcps(rcp)) %>%
  ggplot(aes(x = temp, col = prior)) +
  geom_line(stat = "density") +
  labs(x = expression(~degree*C), y = "Density") +
  scale_colour_manual(values = prior_cols) +
  guides(col = guide_legend(nrow = 2)) +
  facet_wrap(~ rcp, labeller = label_parsed) +
  theme(
    axis.line.x = element_line(linetype = 1), ## Temporary bug(?) in cowplot theme: missing axis line
    axis.line.y = element_line(linetype = 1), ## Ditto
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank(),
    # strip.text = element_text(size = 18, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.margin = unit(2, "lines") ## Increase gap between facet panels
    ) +
  ggsave(file = "./TablesFigures/all-2100.pdf",
         width = 6, height = 5.5,
         device = cairo_pdf)

## Summarise in tabular form
all_2100 %>%
  group_by(rcp, prior) %>%
  summarise(mean = decimals(mean(temp), 1),
            q025 = decimals(quantile(temp, .025), 1),
            q975 = decimals(quantile(temp, .975), 1)) %>%
  mutate(prior = factor(prior, levels = c("ni", "lukemod", "lukestrong",
                                          "denmod", "denstrong"))) %>%
  arrange(prior)

#      rcp      prior  mean  q025  q975
#    (chr)     (fctr) (chr) (chr) (chr)
# 1  rcp26         ni   1.0   0.8   1.2
# 2  rcp26    lukemod   1.0   0.8   1.1
# 3  rcp26 lukestrong   0.9   0.7   1.1
# 4  rcp26     denmod   1.0   0.8   1.1
# 5  rcp26  denstrong   0.4   0.0   0.7
# 6  rcp45         ni   1.7   1.5   1.9
# 7  rcp45    lukemod   1.7   1.5   1.8
# 8  rcp45 lukestrong   1.5   1.3   1.7
# 9  rcp45     denmod   1.6   1.5   1.8
# 10 rcp45  denstrong   0.6   0.2   0.9
# 11 rcp60         ni   2.2   2.0   2.4
# 12 rcp60    lukemod   2.2   2.0   2.4
# 13 rcp60 lukestrong   1.9   1.7   2.1
# 14 rcp60     denmod   2.1   1.9   2.3
# 15 rcp60  denstrong   0.7   0.3   1.1
# 16 rcp85         ni   3.4   3.1   3.6
# 17 rcp85    lukemod   3.3   3.1   3.6
# 18 rcp85 lukestrong   2.9   2.7   3.2
# 19 rcp85     denmod   3.3   3.0   3.5
# 20 rcp85  denstrong   1.0   0.6   1.5
