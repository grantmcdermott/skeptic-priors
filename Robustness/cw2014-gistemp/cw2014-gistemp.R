rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
## Total chain length will thus be chain_length * n_chains
chain_length <- 10000
n_chains <- 3

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 


## Choose prior, conviction and RCP types
## Only need noninformative priors and RCP 2.6, respectively for this robustness check
prior_type <- "ni"
convic_type <- ""
rcp_type <- "rcp26"

tcr <- list()

## Loop over the secondary temperature series ##
for (i in 1:3) {  
  
  temp_series <- c("had", "cw" , "giss")[i]
  
  clim_df <-
    climate %>%
    filter(rcp == rcp_type) %>%
    filter(year <= 2005) %>%
    select(had, cw, giss, trf, volc, soi, amo) %>%
    gather(series, temp, had:giss) %>%
    filter(series == temp_series)
  
  ## Using the automated LearnBayes commands
  theta_sample <- blinreg(clim_df$temp, 
                          cbind(alpha = 1, beta = clim_df$trf, gamma = clim_df$volc, 
                                delta = clim_df$soi, eta = clim_df$amo), 
                          chain_length * n_chains)
    
    ## Get coefficients MCMC list into separate matrix for later. Combines all chains into one matrix.##
    beta_coef <- theta_sample$beta[, "Xbeta"]
    
    rm(theta_sample)
    
    ## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
    tcr[[i]] <- data.frame(tcr = beta_coef * rf2x,
                         series = temp_series)
    rm(list = setdiff(ls(), 
                      c("climate", "rf2x", 
                        "match_coefs", "match_priors", "decimals",
                        "chain_length", "n_chains",
                        "prior_type", "rcp_type",
                        "tcr")))
  
} 

tcr <-
  bind_rows(tcr) %>%
  group_by(series)

tcr %>%
  ungroup() %>%
  mutate(series = gsub("cw", "CW2014", series),
         series = gsub("giss", "GISTEMP", series),
         series = gsub("had", "HadCRUT4", series)) %>% 
  mutate(series = factor(series, levels = c("HadCRUT4", "CW2014", "GISTEMP"))) %>%
  ggplot(aes(x = tcr, group = series, col = series)) +
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf,
           alpha = .2) +
  labs(x = expression(~degree*C), y = "Density") +
  xlim(-1, 3) +
  geom_line(stat = "density") +
  labs(x = expression(~degree*C), y = "Density") +
  theme(text = element_text(family = "Palatino Linotype"),
        legend.title=element_blank(),
        legend.position="bottom"
        ) +
  ggsave(file = "./Robustness/TablesFigures/tcr-other.pdf",
         width = 5, height = 4, 
         device = cairo_pdf) ## See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466

  

tcr %>%
  summarise(mean = round(mean(tcr), 2),
            q025 = round(quantile(tcr, .025), 2),
            q975 = round(quantile(tcr, .975), 2)) %>%
  arrange(mean)

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
tcr_other <- 
  tcr %>%
  split(.$series) %>%
  map(~ clean_func(.$tcr)) %>% 
  transpose() %>%
  do.call("cbind", .) %>%
  # `colnames<-`(c("Mean", "95% C.I.")) %>% ## Works, but next line more elegant
  magrittr::set_colnames(c("Mean", "95% credible interval")) %>%
  as.data.frame() %>%
  mutate(Series = row.names(.)) %>%
  mutate(Series = gsub("cw", "CW2014", Series),
         Series = gsub("giss", "GISTEMP", Series),
         Series = gsub("had", "HadCRUT4", Series)) %>%
  select(Series, everything()) %>% 
  slice(c(3,1,2))
tcr_other$"Effective sample period" <- c(rep("1866-2005", 2), "1880-2005")

## Requires some addtional; tinkering in LaTeX (threeparttable package, etc)
stargazer(tcr_other,
          header = F,
          summary = F,
          rownames = F,
          title ="Transient climate response (TCR), $^\\circ$C: Results from different temperature series",
          label = "tab:tcr-other",
          # align = T, 
          notes.align = "l",
          notes = c("\\footnotesize All estimates are computed using noninformative priors."),
          # type = "text" 
          out = "./Robustness/TablesFigures/tcr-other.tex"
          )

# tcr_other %>%
#   xtable::xtable() %>%
#   print(caption = "Transient climate response (TCR), $^\\circ$C: Results from different temperature series", 
#         caption.placement = 'top',
#         label = "tab:tcr-other",
#         include.rownames=F,
#         booktabs = T)

