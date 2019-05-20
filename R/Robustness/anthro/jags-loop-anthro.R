## Loop over all four RCPs ##
rcp_loop <-
  
  pblapply(c("rcp26" , "rcp45" , "rcp60" , "rcp85"), function(i){

  ## Subset the data to the relevant RCP ##
  clim_df <- 
    climate %>%
    filter(rcp == i) 
  
  N <- nrow(clim_df)
  
  ##------------------------------------------------------------------------------
  ## THE BUGS/JAGS MODEL.
  
  bugs_file <- bugs_model_func_anthro ## CHANGED
  
  
  ##------------------------------------------------------------------------------
  ## SPECIFY THE DATA, INITIALIZATION VALUES AND PARAMETERS OF INTEREST.
  
  ## Tell JAGS where the data are coming from
  data_list <- 
    list(
      "N" = N, "had" = clim_df$had, "trf_less_anthro" = clim_df$trf_less_anthro, 
      "anthro" = clim_df$anthro, ## NEW
      "volc" = clim_df$volc_mean, "soi" = clim_df$soi_mean, "amo" = clim_df$amo_mean,
      "mu_beta" = mu_beta, "sigma_beta" = sigma_beta
      )
  
  ## Give JAGS some initialization values for the model parameters
  inits_list <- 
    function() {
      list(alpha = 0, beta_anthro = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1, phi = 0.5) ## CHANGED
      }
  
  ## Which parameters should R keep track of (i.e. return the posterior distributions for)?
  parameters <- c("alpha", "beta_anthro", "beta", "gamma", "delta", "eta", "sigma", "phi", "y_pred") ## CHANGED
  
  
  ##------------------------------------------------------------------------------
  ## RUN THE PARALLEL JAGS MODEL.

  mod_samples <- 
    jags_par_model(
      bugs_file=bugs_file, data_list=data_list, inits_list=inits_list, parameters=parameters
      )

  
  ##------------------------------------------------------------------------------
  ## SUMMARY AND DIAGNOSTICS
  
  # head(mod_samples)
  # summary(mod_samples)
  # quantile(mod_samples, probs = c(5, 50, 95)/100)
  #   
  # mod_samples[,1:6] %>% gelman.plot() 
  # mod_samples[,1:3] %>% plot(); mod_samples[,4:6] %>% plot()
  # geweke.diag(mod_samples)
  # heidel.diag(mod_samples)
  # raftery.diag(mod_samples)
  
  ## Extract regression coefficients ##
  ## (NOTE: Use only those from RCP 2.6 for consistency) ##
  if (i == "rcp26") {
    
    ## Convert coefficients MCMC list into data frame for later. First combines
    ## all chains into one matrix.
    coefs_df <-
      as.matrix(mod_samples[, c(1:8)], iters = F) %>% ## CHANGED
      # data.frame() %>% 
      as_tibble() %>%
      gather(coef, values)
    
    ## Get summary statistics for tables ##
    coefs_tab <-
      coefs_df %>% 
      group_by(coef) %>% 
      summarise(
        mean = mean(values), 
        q025 = quantile(values, 0.025),
        q975 = quantile(values, 0.975)
        ) %>% 
      mutate(coef = factor(coef, levels=c("beta_anthro", "beta","gamma","delta","eta","alpha","sigma","phi"))) %>% 
      mutate(prior = paste0(prior_type, convic_type)) %>%
      arrange(coef)
    
    ## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
    tcr <- 
      tibble(
        beta_anthro = filter(coefs_df, coef == "beta_anthro")$values, ## CHANGED
        prior = paste0(prior_type, convic_type)
        ) 
    tcr$tcr <- tcr$beta_anthro * rf2x ## CHANGED 
    
  
    # if (prior_type=="ni") {
    #   #############################################
    #   ### Figure S1: Coefficient densities plot ###
    #   #############################################
    #   fig_s1 <- coefs_plot(coefs_df)
    #   fig_s1 +
    #     ggsave(
    #       file = here("TablesFigures/Untracked/Robustness/PNGs/fig-s1-anthro.png"),
    #       width = 8, height = 10
    #       )
    #   fig_s1 +
    #     ggsave(
    #       file = here("TablesFigures/Untracked/Robustness/fig-s1-anthro.pdf"),
    #       width = 8, height = 10, 
    #       device = cairo_pdf
    #       )
    #   rm(fig_s1)
    # } ## End of mini NI "if" clause
    
    rm(coefs_df)
    
  } ## End of RCP 2.6 "if" clause
  
  if (i != "rcp26") {
    coefs_tab <- NULL
    tcr <- NULL
  }
  
  ## Summarise temperature predictions over 1866-2100 ##
  predictions <- jagsresults(mod_samples, params = "y_pred", regex = T)
  predictions <- as_tibble(predictions[, c("mean", "2.5%", "97.5%")])
  colnames(predictions) <- c("mean", "q025", "q975")
  
  predictions$series <- i #rcp_type
  predictions$year <- seq(from=1866, length.out=nrow(predictions))
  predictions <- 
    predictions %>% 
    gather(stat, temp, -c(year, series)) %>%
    select(year, everything())
  
  ## Full distribution of temps in 2100 by themselves 
  temp2100 <- as_tibble(as.matrix(mod_samples[, "y_pred[235]"]))
  colnames(temp2100) <- "temp"
  temp2100$rcp <- i #rcp_type 
  temp2100$prior <- paste0(prior_type, convic_type) 
  
  return(list(tcr=tcr, coefs_tab=coefs_tab,
              predictions=predictions, temp2100=temp2100, 
              N=data.frame(N))) 
  
  }) ## END OF RCP LOOP FOR JAGS SIMLUATIONS

## Recombine the sub-elements of the list based on their common indexes
rcp_loop <- 
  do.call(function(...) mapply(bind_rows, ..., SIMPLIFY = F), args = rcp_loop)

## Extract/copy some data frames within the rcp_loop list to the (local) global environment
predictions <- rcp_loop$predictions
N <- rcp_loop$N[1,"N"]

## Add historic temperature obs to predictions data frame
predictions <- 
  bind_rows(
    tibble(
      year = climate$year[1:N],
      series = "had_full",
      stat = "mean",
      temp = climate$had_full[1:N]
      ),
    predictions
    ) 

## Tidy the data:
## Get rid of duplicate historic (pre-2006) model fits from different RCPs. 
## Simultaneously rename historic RCP 2.6 series as "fitted".
## Lastly, order series as factor and define readable labels for plotting
predictions <- 
  predictions %>%
  mutate(series = ifelse(year <= 2005, gsub("rcp26", "fitted", series), series)) %>%
  filter(year >= 2005 | series %in% c("had_full", "fitted")) %>%
  filter(!(year > 2005 & series %in% c("fitted"))) %>% 
  spread(stat, temp) %>% 
  arrange(series) %>%
  mutate(series = factor(series, levels=c("had_full","fitted","rcp26","rcp45","rcp60","rcp85"))) 


##########################################
### Figure 4: Model fit and prediction ###
##########################################

fig_4 <- pred_plot(predictions)
fig_4_dir <- here("TablesFigures/Untracked/Robustness/")
fig_4_lab <- paste0("fig-4-", prior_type, convic_type, "-anthro")
fig_4 +
  ggsave(
    file = paste0(fig_4_dir, "PNGs/", fig_4_lab, ".png"),
    width = 8, height = 6
    )
fig_4 +
  ggsave(
    file = paste0(fig_4_dir, fig_4_lab, ".pdf"),
    width = 8, height = 6,
    device = cairo_pdf 
    )
rm(fig_4, fig_4_dir, fig_4_lab)

## Remove data frames no longer needed
rm(N, predictions)
## Similarly, subset rcp_loop list to relevant variables for outer (prior) loop
rcp_loop <- rcp_loop[c("coefs_tab", "tcr", "temp2100")]

return(rcp_loop=rcp_loop)
