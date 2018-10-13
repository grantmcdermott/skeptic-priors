rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("R/sceptic_funcs.R")

## Uncomment to run the individual robustness files one by one
# source("Robustness/CW2014-GISTEMP/cw-2014-gsitemp.R")
# source("Robustness/MeasError/sceptic-me.R")
# source("Robustness/Marvel/marvel.R")

## Combine the main (noninformative) TCR outputs in a single table
tcr_robust <- 
  lapply(list.files("Data/Robustness", full.names = T), read_csv) %>%
  bind_rows()

tcr_robust <-
  tcr_robust %>%
  mutate(series = 
           factor(series, 
                  levels = c("had", "cw", "giss", "me", "marvel_a", "marvel_b"))) %>%
  arrange(series)

## Format for LaTeX export
tcr_robust_tab <-
  tcr_robust %>%
  rename_("Mean" = "mean", "Series" = "series") %>%
  mutate(Mean = decimals(Mean, 1),
         "95% C.I." = paste0("[", 
                           decimals(q025, 1), ", ", 
                           decimals(q975, 1),
                           "]")
         ) %>%
  mutate(Series = gsub("had", "HadCRUT4", Series),
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
    "Adjusted forcing efficacies (distributions).")

## Requires some addtional tinkering in LaTeX (threeparttable package, left-align, etc)
stargazer(tcr_robust_tab,
          header = F,
          summary = F,
          rownames = F,
          title ="TCR: Alternative specifications and robustness checks",
          label = "tab:tcr-robust",
          # align = T, 
          notes.align = "l",
          notes = c("\\footnotesize All estimates are computed using noninformative priors. See text for details."),
          # type = "text" 
          out = "Robustness/TablesFigures/tcr-robust.tex")