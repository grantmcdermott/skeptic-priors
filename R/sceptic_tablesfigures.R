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


###############################
### Figure 1: TCR densities ###
###############################

## Priors and posteriors
fig_1 <- tcr_plot(tcr)
fig_1 +
  ggsave(
    file = paste0(pref, "PNGs/fig-1", suff, ".png"),
    width = 8, height = 4.5
    )
fig_1 +
  ggsave(
    file = paste0(pref, "fig-1", suff, ".pdf"),
    width = 6, height = 4,
    device = cairo_pdf
    )
rm(fig_1)

## Just the priors this time (for presentations)
tcr_plot_priors <- tcr_plot_priors(tcr)
tcr_plot_priors +
  ggsave(
    file = paste0(pref, "PNGs/tcr-prior", suff, ".png"),
    width = 6, height = 4
    )
tcr_plot_priors +
  ggsave(
    file = paste0(pref, "tcr-prior", suff, ".pdf"),
    width = 6, height = 4,
    device = cairo_pdf
    )
rm(tcr_plot_priors)

## Summarise in tabular form
tcr %>%
  group_by(prior) %>%
  summarise(tcr_mean = decimals(mean(tcr), 1),
            q025 = decimals(quantile(tcr, .025), 1),
            q975 = decimals(quantile(tcr, .975), 1)) %>%
  mutate(prior = match_priors(prior)) %>% 
  arrange(desc(tcr_mean))

####################
### TEMP IN 2100 ###
####################

## All 2100 pointrange plot

all_2100_plot <- all_2100_plot_func(all_2100)
all_2100_plot +
  ggsave(
    file = paste0(pref, "PNGs/all-2100", suff, ".png"),
    width = 9, height = 6
    )
all_2100_plot +
  ggsave(
    file = paste0(pref, "all-2100", suff, ".pdf"),
    width = 9, height = 6,
    device = cairo_pdf
    )
rm(all_2100_plot)

## Summarise in tabular form
all_2100 %>%
  group_by(rcp, prior) %>%
  summarise(mean_2100 = decimals(mean(temp), 1),
            q025 = decimals(quantile(temp, .025), 1),
            q975 = decimals(quantile(temp, .975), 1)) %>%
  mutate(prior = match_priors(prior)) %>% 
  group_by(prior) %>%
  arrange(desc(mean_2100))
