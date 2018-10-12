## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

scc <- read_csv("Data/PAGE09/scc.csv")

fig_s2 <- scc_plot(scc)
fig_s2 +
  ggsave(
    file = "TablesFigures/PNGs/scc.png",
    width = 6, height = 4.5
    )
fig_s2 +
  ggsave(
    file = "TablesFigures/scc.pdf",
    width = 6, height = 4.5,
    device = cairo_pdf
    )
rm(fig_s2)

scc_tab <-
  scc %>%
  gather(prior, scc) %>%
  group_by(prior) %>%
  summarise(Mean = decimals(mean(scc), 2),
            Median = decimals(quantile(scc, .5), 2),
            q025 = decimals(quantile(scc, .025), 2),
            q975 = decimals(quantile(scc, .975), 2)) %>%
  mutate(prior = factor(match_priors(prior), levels = rev(prior_names))) %>%
  arrange(prior)
scc_tab <-
  scc_tab %>%
  mutate("95% probability interval" = paste0("[", sprintf("%.2f", q025), 
                                             ", ", 
                                             sprintf("%.2f", q975), "]")
         ) %>%
  select(-c(q025, q975)) %>%
  magrittr::set_colnames(c("", "Mean", "Median", "95% Probability Interval"))

## Export table to LaTeX. Still requires some manual tinkering to get ideal 
## formatting and alignment, as well as include table notes.
scc_tab %>%
  xtable(#align = c("l", "l","c","c","c"), ## Note: extra col align. char. (yet to exlude row names)
         caption = "Social cost of carbon (US\\$2005 per tonne)",
         label = "tab:scc") %>%
  print(booktabs = T, caption.placement = "top", 
        table.placement = "t", include.rownames = F,
        file = "TablesFigures/scc.tex"
        )

## Using xable.decimal function for alignment (see scep_funcs.R file)
## Not quite right either
# scc_tab %>%
#   xtable.decimal(caption = "Social cost of carbon (US\\$2005 per tonne)",
#                  label = "tab:scc",
#                  booktabs = T, caption.placement = "top",
#                  table.placement = "t"#,
#                  #colAlignment = c("95% probability interval" = "c"), 
#                  )


## Alternatively, use Stargazer package... but still requires additional tinkering
## in LaTeX and also doesn't have booktabs option.
# library(stargazer) 
scc_tab %>%
  magrittr::set_colnames(c("prior", "Mean", "Median", "95% Probability Interval")) %>%
  mutate(prior = as.character(prior)) %>%
  magrittr::set_colnames(c("", "Mean", "Median", "95% Probability Interval")) %>%
  stargazer(summary = F,
          header = F,
          title = "Social cost of carbon (US\\$2005)",
          label = "tab:scc",
          rownames = F,
          align = T,  
          notes.align = "l"#,
          # out = "TablesFigures/scc.tex"
          )
