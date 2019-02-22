## Default TablesFigures directory prefix if not already set by run type
pref <- ifelse(!exists("pref"), here("TablesFigures/"), pref)

##########################
##########################
####    * TABLES *    ####
##########################
##########################

## Only export tables for the main run
if (run_type == "main") {


  ###################################################
  ### Table 3: Posterior regressions coefficients ###
  ###################################################
  
  ## Add TCR summary info to coefficents table
  coefs_tab <-
    bind_rows(
      coefs_tab,
      tcr %>%
        group_by(prior) %>%
        summarise(
          mean = mean(tcr),
          q025 = quantile(tcr, .025),
          q975 = quantile(tcr, .975)
          ) %>%
        mutate(coef = "tcr") %>%
        select(coef, mean, q025, q975, prior)
      )
  
  ## Convert table into nicer format for the paper. Requires some minor tinkering in
  ## in LaTeX afterwards (e.g. not all vars have 3 decimals), but close enough.
  coefs_tab <- 
    coefs_tab %>%
    filter(coef != "sigma") %>%
    mutate(
      mean = ifelse(coef=="tcr", decimals(mean, 1), decimals(mean, 3)),
      ci = ifelse(coef=="tcr",
                  paste0("[", decimals(q025, 1), ", ", decimals(q975, 1), "]"),
                  paste0("[", decimals(q025, 3), ", ", decimals(q975, 3), "]")
                  )
      ) %>%
    select(prior, coef, mean, ci) %>%
    gather(key, value, -c(prior, coef)) %>% ## This step causes some values to drop the (zero) 3rd decimal
    mutate(
      prior = factor(prior, levels = c("ni", "lukemod", "lukestrong", "denmod", "denstrong")),
      key = factor(key, levels = c("mean", "ci"))
      ) %>%
    arrange(prior, coef) %>%
    spread(prior, value) %>% 
    mutate(coef = factor(coef, levels = c("beta", "gamma", "delta", "eta", "alpha", "tcr"))) %>%
    arrange(coef) %>% 
    mutate(coef = ifelse(key == "mean", paste(coef), "")) %>%
    select(-key) %>%
    as.matrix()
  
  rownames(coefs_tab) <- match_coefs(coefs_tab[, "coef"])
  colnames(coefs_tab) <- match_priors(colnames(coefs_tab))
  
  coefs_tab[, 2:ncol(coefs_tab)] %>%
    stargazer(
      align = T, header = F, rownames = T,
      title = "Posterior regression results",
      label = paste0("tab:reg-coefs", suff),
      # notes.align = "l",
      # notes.append = T,
      # notes = c("Dependent variable: Global Mean Surface Temperature (GMST)."),
      out = paste0(pref, "tab-3", suff, ".tex")
      )
  
  rm(coefs_tab)
  
  
  ########################################
  ### Table 4: Robustness checks (TCR) ###
  ########################################
  
  ## Read and combine the main (noninformative) TCR outputs in a single table
  tcr_robust <- 
    lapply(list.files(here("Results/Robustness"), full.names = T), read_csv) %>%
    bind_rows()
  
  tcr_robust <-
    tcr_robust %>%
    mutate(series = factor(series, levels=c("had","cw","giss","me","marvel_a","marvel_b"))) %>%
    arrange(series)
  
  ## Format for LaTeX export
  tcr_robust_tab <-
    tcr_robust %>%
    rename("Mean" = "mean", "Series" = "series") %>%
    mutate(
      Mean = 
        decimals(Mean, 1), "95% C.I." = paste0("[", decimals(q025, 1), ", ", decimals(q975, 1), "]")
      ) %>%
    mutate(
      Series = gsub("had", "HadCRUT4", Series),
      Series = gsub("cw", "CW2014", Series),
      Series = gsub("giss", "GISTEMP", Series),
      Series = gsub("me", "Measurement Error", Series),
      Series = gsub("marvel", "Marvel", Series)
      ) %>%
    mutate(Comment = "") %>%
    # magrittr::set_rownames(.$series) %>%
    select(-c(q025, q975)) %>%
    as.matrix() 
  
  tcr_robust_tab[, "Comment"] <-
    c("Primary GMST series. For comparison.",
      "Alternative GMST series.",
      "Alternative GMST series.",
      "Specifying measurement error in HadCRUT4.",
      "Adjusted forcing efficacies (means).",
      "Adjusted forcing efficacies (distributions)."
      )
  
  ## Requires some addtional tinkering in LaTeX (threeparttable package, left-align, etc)
  tcr_robust_tab %>%
    stargazer(
      header = F,
      summary = F,
      rownames = F,
      title ="TCR: Alternative specifications and robustness checks",
      label = "tab:tcr-robust",
      # align = T, 
      notes.align = "l",
      notes = c("\\footnotesize All estimates are computed using noninformative priors. See text for details."),
      # type = "text"
      out = paste0(pref, "tab-4.tex")
      )
  
  rm(tcr_robust, tcr_robust_tab)
  
  
  ####################
  ### Table 6: SCC ###
  ####################
  
  scc <- read_csv(here("Results/PAGE09/scc.csv"))
  
  scc_tab <-
    scc %>%
    gather(prior, scc) %>%
    group_by(prior) %>%
    summarise(
      Mean = decimals(mean(scc), 2),
      Median = decimals(quantile(scc, .5), 2),
      q025 = decimals(quantile(scc, .025), 2),
      q975 = decimals(quantile(scc, .975), 2)
      ) %>%
    mutate(prior = factor(match_priors(prior), levels = rev(prior_names))) %>%
    arrange(prior)
  
  scc_tab <-
    scc_tab %>%
    mutate("95% probability interval" = 
             paste0("[", sprintf("%.2f", q025), ", ", sprintf("%.2f", q975), "]")
           ) %>%
    select(-c(q025, q975)) %>%
    magrittr::set_colnames(c("", "Mean", "Median", "95% Probability Interval"))
  
  ## Export table to LaTeX. Still requires some manual tinkering to get ideal 
  ## formatting and alignment, as well as include table notes.
  scc_tab %>%
    xtable(
      #align = c("l", "l","c","c","c"), ## Note: extra col align. char. (yet to exlude row names)
      caption = "Social cost of carbon (US\\$2005 per tonne)",
      label = "tab:scc"
      ) %>%
    print(
      booktabs = T, caption.placement = "top", 
      table.placement = "t", include.rownames = F,
      file = paste0(pref, "tab-6.tex")
      )
  
  rm(scc, scc_tab)

}


###########################
###########################
####    * FIGURES *    ####
###########################
###########################


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
## Only plot this figure for the main run
if (run_type == "main") {
  fig_1_priors <- tcr_plot_priors(tcr)
  fig_1_priors +
    ggsave(
      file = paste0(pref, "Untracked/PNGs/fig-1-priors", suff, ".png"),
      width = 6, height = 4
      )
  fig_1_priors +
    ggsave(
      file = paste0(pref, "Untracked/fig-1-priors", suff, ".pdf"),
      width = 6, height = 4,
      device = cairo_pdf
      )
  rm(fig_1_priors)
}

## Summarise in tabular form
tcr %>%
  group_by(prior) %>%
  summarise(
    tcr_mean = decimals(mean(tcr), 1),
    q025 = decimals(quantile(tcr, .025), 1),
    q975 = decimals(quantile(tcr, .975), 1)
    ) %>%
  mutate(prior = match_priors(prior)) %>% 
  arrange(desc(tcr_mean)) %>% 
  mutate(run_type = run_type)


#########################################
### Figure 2: Recursive TCR estimates ###
#########################################

## Only plot this figure for the main run
if (run_type == "main") {
  lapply(c("historic", "future"), function(recurse_type) {
    tcr_rec <- read_csv(here(paste0("Results/Recursive/tcr-rec-", recurse_type, ".csv")))
    fig_2 <- recursive_plot(tcr_rec)
    fig_2_pref <- here("TablesFigures/Untracked/")
    fig_2_suff <- paste0("-", recurse_type)
    if(recurse_type == "historic"){
      fig_2 <- 
        fig_2 + 
        scale_x_reverse(breaks = seq(max(tcr_rec$year_to), min(tcr_rec$year_to), by = -30))
      fig_2_pref <- here("TablesFigures/")
      fig_2_suff <- ""
    }
    fig_2 +
      ggsave(
        file = paste0(fig_2_pref, "PNGs/fig-2", fig_2_suff, ".png"),
        width = 8, height = 7
        )
    fig_2 +
      ggsave(
        file = paste0(fig_2_pref, "fig-2", fig_2_suff, ".pdf"),
        width = 8, height = 7,
        device = cairo_pdf
        )
    rm(fig_2)
  })
}


###################################
### Figure 3: Evidence required ###
###################################  

if (run_type=="main") {
  ## Read data
  evid <- read_csv(here("Results/Evidence/tcr-evidence.csv"))
  
  ## Plot the data
  ## Years with red-white-blue colour scheme
  fig_3 <- evid_plot(evid)
  fig_3 +
    ggsave(
      file = paste0(pref, "PNGs/fig-3.png"),
      width = 8, height = 4
      )
  fig_3 +
    ggsave(
      file = paste0(pref, "fig-3.pdf"),
      width = 8, height = 4,
      device = cairo_pdf
      )
  rm(fig_3)
  
  ## Lines instead of grid
  fig_3_lines <- evid_plot_lines(evid) 
  fig_3_lines +
    ggsave(
      file = paste0(pref, "Untracked/PNGs/fig-3-lines.png"),
      width = 8, height = 4
    )
  fig_3_lines +
    ggsave(
      file = paste0(pref, "Untracked/fig-3-lines.pdf"),
      width = 8, height = 4,
      device = cairo_pdf
    )
  rm(fig_3_lines)
}
  

###########################################
### Figure 4: Model fit and predictions ###
########################################### 

## Figure(s) already exported as part of the main Bayesian (RCP) loop.


###############################
### Figure 5: Temps in 2100 ###
###############################  

## All 2100 pointrange plot

fig_5 <- temp2100_plot(temp2100)
fig_5 +
  ggsave(
    file = paste0(pref, "PNGs/fig-5", suff, ".png"),
    width = 9, height = 6
    )
fig_5 +
  ggsave(
    file = paste0(pref, "fig-5", suff, ".pdf"),
    width = 9, height = 6,
    device = cairo_pdf
    )
rm(fig_5)

## Summarise in tabular form
temp2100 %>%
  group_by(rcp, prior) %>%
  summarise(
    mean_2100 = decimals(mean(temp), 1),
    q025 = decimals(quantile(temp, .025), 1),
    q975 = decimals(quantile(temp, .975), 1)
    ) %>%
  mutate(prior = match_priors(prior)) %>% 
  # group_by(rcp) %>%
  arrange(desc(rcp), desc(mean_2100)) %>%
  mutate(run_type = run_type)



################################################
### Figure S1: Parameter posterior densities ###
################################################ 

## Figure(s) already exported as part of the main Bayesian (RCP) loop.


################################
### Figure S2: SCC densities ###
################################

if (run_type=="main") {
  
  scc <- read_csv(here("Results/PAGE09/scc.csv"))
  
  fig_s2 <- scc_plot(scc)
  fig_s2 +
    ggsave(
      file = paste0(pref, "PNGs/fig-s2.png"),
      width = 6, height = 4.5
      )
  fig_s2 +
    ggsave(
      file = paste0(pref, "fig-s2.pdf"),
      width = 6, height = 4.5,
      device = cairo_pdf
      )
  rm(scc, fig_s2)
  
  }
