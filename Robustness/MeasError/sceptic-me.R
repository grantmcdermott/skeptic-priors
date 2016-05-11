rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## The 95% measurement error bounds are not quite equal, but very close. We can also compare 
## measurement error with model uncertainty i.e. sigma from noninformative model.
## Note, divide 95% ME bounds by 2 to get 1 std dev (i.e. fair comparison to model sigma).
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

## Load some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

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


##################################
### COMBINED TABLES AND GRAPHS ###
##################################
pref <- "./Robustness/TablesFigures/"
suff <- "-me"

source("sceptic_tablesfigures.R")

# ## Set radiative forcing distribution used for calulating TCRs later in code.
# ## Centered around 3.71 °C +/- 10% (within 95% CI). 
# ## Length of disbn equals length of MCMC chain for consistency
# rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 
# 
# ## Remove data unnecessary to further analysis ##
# rm(climate, n_chains, chain_length, prior_type, convic_type,
#    theme_coefs, theme_pred)
# 
# 
# ####################
# ### COEFFICIENTS ###
# ####################
# 
# ## Unlist coefficents list
# coefs_tab <- bind_rows(coefs_tab) 
# ## Convert into nice table format for the paper
# coefs_tab <- 
#   coefs_tab %>%
#   filter(coef != "sigma") %>%
#   mutate(mean = decimals(mean, 3),
#          ci = paste0("[", decimals(q025, 3), ", ", decimals(q975, 3), "]")
#          ) %>%
#   select(prior, coef, mean, ci) %>%
#   gather(key, value, -c(prior, coef)) %>%
#   mutate(prior = factor(prior, levels = c("ni", "lukemod", 
#                                           "lukestrong", "denmod", "denstrong")),
#          key = factor(key, levels = c("mean", "ci"))) %>%
#   arrange(prior, coef) %>%
#   spread(prior, value) %>%
#   mutate(coef = ifelse(key == "mean", paste(coef), "")) %>%
#   select(-key) %>%
#   as.matrix()
# rownames(coefs_tab) <- match_coefs(coefs_tab[, "coef"])
# colnames(coefs_tab) <- match_priors(colnames(coefs_tab))
# 
# coefs_tab[, 2:ncol(coefs_tab)] %>%
#   stargazer(align = T, header = F, rownames = T,
#             title = "Posterior regression results",
#             label = "tab:coefs-all-me",
#             #             notes.align = "l",
#             #             notes.append = T,
#             #             notes = c("Dependent variable: Global Mean Surface Temperature (GMST).")
#             out = "./Robustness/TablesFigures/coefs-all-me.tex"
#             )
# 
# 
# ###################
# ###     TCR     ###
# ###################
# 
# ## Unlist tcr list and generate values
# tcr <-
#   bind_rows(
#     mapply(
#       function(x) {
#         tcr[[x]]$tcr <- tcr[[x]]$beta * rf2x
#         return(tcr[[x]])
#         },
#       seq(1:5), 
#       SIMPLIFY = F
#       )
#     ) 
# 
# ## TCR Density plot
# tcr %>%
#   mutate(prior = factor(match_priors(prior),
#                         levels = prior_names)) %>%
#   ggplot(aes(x = tcr, col = prior)) +
#   geom_line(stat = "density") +
#   labs(x = expression(~degree*C), y = "Density") +
#   xlim(-1, 3) + 
#   annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf,
#            alpha = .2) +
#   stat_function(fun = dnorm, args = list(mean = 0, sd = .065), 
#                 lty=2, col=prior_cols[1]) +
#   stat_function(fun = dnorm, args = list(mean = 0, sd = .25), 
#                 lty=2, col=prior_cols[2]) + 
#   stat_function(fun = dnorm, args = list(mean = 1, sd = .065), 
#                 lty=2, col=prior_cols[3]) +
#   stat_function(fun = dnorm, args = list(mean = 1, sd = .25), 
#                 lty=2, col=prior_cols[4]) + 
#   scale_colour_manual(values = prior_cols) +
#   guides(col = guide_legend(nrow = 2)) +
#   theme_tr +
#   ggsave(file = "./Robustness/TablesFigures/tcr-combined-me.pdf",
#          width = 6, height = 4,
#          device = cairo_pdf)
# 
# ## Table summary of the above ##
# tcr_tab <-
#   tcr %>%
#   group_by(prior) %>%
#   summarise(mean = mean(tcr),
#             q025 = quantile(tcr, .025),
#             q975 = quantile(tcr, .975)) 
# ## Format for nice looking table in LaTeX
# tcr_tab <-
#   tcr_tab %>%
#   mutate(mean = round(mean, 1),
#          ci = paste0("[", round(q025, 1), ", ", round(q975, 1), "]")
#          ) %>%
#   select(prior, mean, ci) %>%
#   gather(key, Posterior, -prior) %>%
#   mutate(prior = factor(prior, levels = c("ni", "lukemod", 
#                                           "lukestrong", "denmod", "denstrong")),
#          key = factor(key, levels = c("mean", "ci"))) %>%
#   arrange(prior) %>%
#   mutate(prior = ifelse(key == "mean", paste(prior), "")) %>%
#   select(-key) %>%
#   rename(Prior = prior) %>%
#   as.matrix()
# rownames(tcr_tab) <- match_priors(tcr_tab[, "Prior"])
# tcr_tab[, "Prior"] <- c("-", "",
#                         "1.0","[0.5, 1.5]",
#                         "1.0","[0.9, 1.1]",
#                         "0.0","[-0.5, 0.5]",
#                         "0.0","[-0.1, 0.1]")
# 
# stargazer(tcr_tab,
#           title ="Transient climate response (TCR), $^\\circ$C",
#           header = F,
#           label = "tab:tcr-tab",
#           # align = T,
#           notes.align = "l",
#           notes = c(paste0("\\footnotesize Mean estimates are given, with square ",
#                            "brackets denoting 95\\% probability interval.")),
#           out = "./Robustness/TablesFigures/tcr-tab-me.tex"
#           )
# 
# 
# ####################
# ### TEMP IN 2100 ###
# ####################
# 
# ## Combined temperature prediction in 2100 by scenario ###
# ## Unlist temperature in 2100 list
# all_2100 <-
#   bind_rows(all_2100) 
# 
# ## All 2100 density plot
# all_2100 %>%
#   mutate(prior = factor(match_priors(prior),
#                         levels = prior_names)) %>%
#   mutate(rcp = match_rcps(rcp)) %>%
#   ggplot(aes(x = temp, col = prior)) +
#   geom_line(stat = "density") +
#   labs(x = expression(~degree*C), y = "Density") +
#   scale_colour_manual(values = prior_cols) +
#   guides(col = guide_legend(nrow = 2)) +
#   facet_wrap(~ rcp, labeller = label_parsed) +
#   theme_2100 +
#   ggsave(file = "./Robustness/TablesFigures/all-2100-me.pdf",
#          width = 6, height = 5.5,
#          device = cairo_pdf)
# 
# ## Summarise in tabular form
# all_2100 %>%
#   group_by(rcp, prior) %>%
#   summarise(mean = decimals(mean(temp), 1),
#             q025 = decimals(quantile(temp, .025), 1),
#             q975 = decimals(quantile(temp, .975), 1)) %>%
#   mutate(prior = factor(prior, levels = c("ni", "lukemod", "lukestrong",
#                                           "denmod", "denstrong"))) %>%
#   arrange(prior)

#      rcp      prior  mean  q025  q975
#    (chr)     (fctr) (dbl) (dbl) (dbl)
# 1  rcp26         ni   1.0   0.8   1.1
# 2  rcp26    lukemod   1.0   0.8   1.1
# 3  rcp26 lukestrong   0.9   0.7   1.0
# 4  rcp26     denmod   1.0   0.8   1.1
# 5  rcp26  denstrong   0.4   0.0   0.7
# 6  rcp45         ni   1.7   1.5   1.8
# 7  rcp45    lukemod   1.7   1.5   1.8
# 8  rcp45 lukestrong   1.5   1.3   1.6
# 9  rcp45     denmod   1.6   1.5   1.8
# 10 rcp45  denstrong   0.6   0.2   0.9
# 11 rcp60         ni   2.2   2.0   2.4
# 12 rcp60    lukemod   2.2   2.0   2.3
# 13 rcp60 lukestrong   1.9   1.7   2.1
# 14 rcp60     denmod   2.1   1.9   2.3
# 15 rcp60  denstrong   0.7   0.3   1.1
# 16 rcp85         ni   3.4   3.1   3.6
# 17 rcp85    lukemod   3.3   3.1   3.5
# 18 rcp85 lukestrong   2.9   2.7   3.2
# 19 rcp85     denmod   3.2   3.0   3.5
# 20 rcp85  denstrong   1.0   0.6   1.5
