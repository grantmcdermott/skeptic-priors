## Load all packages, as well as some helper functions that will be used for plotting and tables
source(here::here("R/sceptic_funcs.R"))

## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv(here("Data/climate.csv"))

## Decide on total length of MCMC chains (i.e. summed parallel chains of the JAGS model)
## Each individual chain will thus be chain_length/n_chains.
chain_length <- 30000

## Now set the JAGS parallel MCMC parameters
n_chains <- detectCores() ## no. of parallel chains. Use `detectCores()-1` if you are worried about CPU resources
n_adapt = 5000 ## no. of tuning or adaptation steps
burn_in = 1000 ## no. of burn-in steps

## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 


## Choose prior, conviction and RCP types
## Only need noninformative priors and RCP 2.6, respectively for this robustness check
prior_type <- "ni"
convic_type <- ""
rcp_type <- "rcp26"


tcr_secondary <-
  lapply(c("had", "cw" , "giss"), function(x){
    
    clim_df <-
      climate %>%
      filter(rcp == rcp_type) %>%
      filter(year <= 2005) %>%
      select(had, cw, giss, trf, volc, soi_mean, amo_mean) %>%
      gather(series, temp, had:giss) %>%
      filter(series == x) %>%
      filter(!is.na(temp))
    
    ## Get noninformative priors
    m <- filter(priors_df, prior_type=="ni")$mu
    s <- filter(priors_df, prior_type=="ni")$sigma
    
    mu_beta <- m/3.71
    sigma_beta <- s/3.71
    
    N <- nrow(clim_df)
    
    ##------------------------------------------------------------------------------
    ## THE BUGS/JAGS MODEL.
    bugs_file <- bugs_model_func
    
    ## Tell JAGS where the data are coming from
    ## Notes: Using "clim_df$temp" for the "had" variable, since this is actually changing
    ## depending on the data source. (Easier than re-writing the BUGS model func.)
    data_list <- 
      list(
        "N" = N, "had" = clim_df$temp, "trf" = clim_df$trf, 
        "volc" = clim_df$volc, "soi" = clim_df$soi_mean, "amo" = clim_df$amo_mean,
        "mu_beta" = mu_beta, "sigma_beta" = sigma_beta
        )
    
    ## Give JAGS some initialization values for the model parameters
    inits_list <- 
      function() {
        list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1, phi = 0.5)
        }
    
    ## Which parameters should R keep track of (i.e. return the posterior distributions for)?
    parameters <- c("beta")
    
    ## Run the parallel JAGS model
    mod_samples <- 
      jags_par_model(
        bugs_file=bugs_file, data_list=data_list, inits_list=inits_list, parameters=parameters
        )
    
    ## Get TCR
    tcr <- as.vector(as.matrix(mod_samples, iters = F)) * rf2x
    ## Put in a DF
    tcr_df <- tibble(tcr, series = x)
    
    return(tcr_df)
    
  }
  ) %>%
  bind_rows()

tcr_secondary %>%
  ungroup() %>%
  mutate(
    series = gsub("cw", "CW2014", series),
    series = gsub("giss", "GISTEMP", series),
    series = gsub("had", "HadCRUT4", series)
    ) %>% 
  mutate(series = factor(series, levels = c("HadCRUT4", "CW2014", "GISTEMP"))) %>%
  ggplot(aes(x = tcr, group = series, col = series)) +
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf, alpha = .2) +
  labs(x = expression(~degree*C), y = "Density") +
  xlim(-1, 3) +
  geom_line(stat = "density") +
  labs(x = expression(~degree*C), y = "Density") +
  scale_color_brewer(palette = "Set1") +
  theme(
    text = element_text(family = font_type),
    legend.title=element_blank(),
    legend.position="bottom"
    ) +
  ggsave(
    file = here("TablesFigures/Untracked/Robustness/tcr-secondary.pdf"),
    width = 5, height = 4, 
    device = cairo_pdf
    )

tcr_secondary %>%
  group_by(series) %>%
  summarise(mean = mean(tcr),
            q025 = quantile(tcr, .025),
            q975 = quantile(tcr, .975)) %>%
  arrange(mean) %>%
  write_csv(here("Results/Robustness/tcr-secondary.csv"))

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
tcr_secondary_tab <- 
  tcr_secondary %>%
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
tcr_secondary$"Effective sample period" <- c(rep("1866-2005", 2), "1880-2005")

## Requires some addtional tinkering in LaTeX (threeparttable package, etc)
stargazer(tcr_secondary_tab,
          header = F,
          summary = F,
          rownames = F,
          title ="Transient climate response (TCR), $^\\circ$C: Results from different temperature series",
          label = "tab:tcr-secondary",
          # align = T, 
          notes.align = "l",
          notes = c("\\footnotesize All estimates are computed using noninformative priors."),
          # type = "text" 
          out = here("TablesFigures/Untracked/Robustness/tcr-secondary.tex")
          )
