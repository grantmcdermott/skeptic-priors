##########################
### COEFFICIENTS TABLE ###
##########################

## Add TCR summary info to coefficents table
coefs_tab <-
  bind_rows(coefs_tab,
            tcr %>%
              group_by(prior) %>%
              summarise(mean = mean(tcr),
                        q025 = quantile(tcr, .025),
                        q975 = quantile(tcr, .975)) %>%
              mutate(coef = "tcr") %>%
              select(coef, mean, q025, q975, prior)
            )

## Convert table into nicer format for the paper. Requires some minor tinkering in
## in LaTeX afterwards (e.g. not all vars have 3 decimals), but close enough.
coefs_tab <- 
  coefs_tab %>%
  filter(coef != "sigma") %>%
  mutate(mean = ifelse(coef=="tcr", decimals(mean, 1), decimals(mean, 3)),
         ci = ifelse(coef=="tcr",  
                     paste0("[", decimals(q025, 1), ", ", decimals(q975, 1), "]"),
                     paste0("[", decimals(q025, 3), ", ", decimals(q975, 3), "]")
                     )
         ) %>%
  select(prior, coef, mean, ci) %>%
  gather(key, value, -c(prior, coef)) %>% ## This step causes some values to drop the (zero) 3rd decimal
  mutate(prior = factor(prior, levels = c("ni", "lukemod", 
                                          "lukestrong", "denmod", "denstrong")),
         key = factor(key, levels = c("mean", "ci"))) %>%
  arrange(prior, coef) %>%
  spread(prior, value) %>% 
  mutate(coef = factor(coef, levels = c("beta", "gamma", 
                                        "delta", "eta", "alpha", "tcr"))) %>%
  arrange(coef) %>% 
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

## As above, but fill only
tcr %>%
  mutate(prior = factor(match_priors(prior),
                        levels = prior_names)) %>%
  ggplot(aes(x = tcr, fill = prior)) +
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
  geom_area(stat = "density", position = "dodge", alpha = .5) +
  scale_colour_manual(values = prior_cols) +
  scale_fill_manual(values = prior_cols) +
  guides(col = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  theme_tcr +
  ggsave(file = paste0(pref, "tcr-combined", suff, "-fill.pdf"),
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

## Summarise in tabular form
tcr %>%
  group_by(prior) %>%
  summarise(tcr_mean = decimals(mean(tcr), 1),
            q025 = decimals(quantile(tcr, .025), 1),
            q975 = decimals(quantile(tcr, .975), 1)) %>%
  mutate(prior =
           factor(prior,
                  levels = c("ni", "lukemod", "lukestrong", "denmod", "denstrong"))
         ) %>%
  arrange(prior)

####################
### TEMP IN 2100 ###
####################

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

## As above, but fill only
all_2100 %>%
  mutate(prior = factor(match_priors(prior),
                        levels = prior_names)) %>%
  mutate(rcp = match_rcps(rcp)) %>%
  ggplot(aes(x = temp, fill = prior)) +
  geom_area(stat = "density", position = "dodge", alpha = .5) +
  labs(x = expression(~degree*C), y = "Density") +
  scale_fill_manual(values = prior_cols) +
  guides(col = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  facet_wrap(~ rcp, labeller = label_parsed) +
  theme_2100 +
  ggsave(file = paste0(pref, "all-2100", suff, "-fill.pdf"),
         width = 6, height = 5.5,
         device = cairo_pdf)

## Summarise in tabular form
all_2100 %>%
  group_by(rcp, prior) %>%
  summarise(mean_2100 = decimals(mean(temp), 1),
            q025 = decimals(quantile(temp, .025), 1),
            q975 = decimals(quantile(temp, .975), 1)) %>%
  mutate(prior = factor(prior, levels = c("ni", "lukemod", "lukestrong",
                                          "denmod", "denstrong"))) %>%
  arrange(prior)