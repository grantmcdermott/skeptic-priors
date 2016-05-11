## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 

## Remove data unnecessary to further analysis ##
rm(climate, n_chains, chain_length, prior_type, convic_type,
   theme_coefs, theme_pred)


####################
### COEFFICIENTS ###
####################

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

coefs_tab[, 2:ncol(coefs_tab)] %>%
  stargazer(align = T, header = F, rownames = T,
            title = "Posterior regression results",
            label = paste0("tab:coefs-all", suff),
            #             notes.align = "l",
            #             notes.append = T,
            #             notes = c("Dependent variable: Global Mean Surface Temperature (GMST).")
            out = paste0(pref, "coefs-all", suff, ".tex")
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
  theme_tcr +
  ggsave(file = paste0(pref, "tcr-combined", suff, ".pdf"),
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
  theme_tcr +
  ggsave(file = paste0(pref, "tcr-combined-prior", suff, ".pdf"),
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
          label = paste0("tab:tcr-tab", suff),
          # align = T,
          notes.align = "l",
          notes = c(paste0("\\footnotesize Mean estimates are given, with square ",
                           "brackets denoting 95\\% probability interval.")),
          out = paste0(pref, "tcr-tab", suff, ".tex")
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
  theme_2100 +
  ggsave(file = paste0(pref, "all-2100", suff, ".pdf"),
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