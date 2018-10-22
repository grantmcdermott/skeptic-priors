## Loop over all four RCPs ##
rcp_loop <-
  
  pblapply(c("rcp26" , "rcp45" , "rcp60" , "rcp85"), function(i) {
  
  ## Subset the data to the relevant RCP ##
  clim_df <- 
    climate %>%
    filter(rcp == i) 
  
  ## Fill in measurement error for future values (needed for prediction) based on
  ## measurement error of last two decades
  omega_had_recent <- 
    (clim_df %>% 
       filter(!is.na(had_full)) %>% 
       tail(20))$had_omega %>% 
    mean()

  clim_df$had_omega <- ifelse(is.na(clim_df$had_omega), omega_had_recent, clim_df$had_omega)
  
  N <- nrow(clim_df)
  
  ##------------------------------------------------------------------------------
  ## THE BUGS/JAGS MODEL.
  
  bugs_file <- bugs_model_func_me
  
  ##------------------------------------------------------------------------------
  ## PARALLEL SETUP.
  
  cl <- parallel::makeCluster(n_chains, type = cl_type) # no. of clusters (i.e. MCMC chains)
  parLoadModule(cl, "lecuyer", quiet = T)
  clusterSetRNGStream(cl, 123)
  
  ##------------------------------------------------------------------------------
  ## SPECIFY THE DATA AND INITIALIZE THE CHAINS.
  ## Notes: Adding "had_omega" (or cw_omega) to account for measurement errors
  
  data_list <- 
    list(
      "N" = N, "had" = clim_df$had, "trf" = clim_df$trf, 
      "volc" = clim_df$volc_mean, "soi" = clim_df$soi_mean, "amo" = clim_df$amo_mean,
      "mu_beta" = mu_beta, "sigma_beta" = sigma_beta, 
      "had_omega" = clim_df$had_omega
    )
  inits_list <- 
    function() {list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1)}

  ##------------------------------------------------------------------------------
  ## RUN THE CHAINS/MCMC SAMPLES.
  
  ## Which parameters should R keep track of?
  parameters <- c("alpha", "beta", "gamma", "delta", "eta", "sigma", "y_pred")
  ## Initialisation
  par_inits <- parallel.inits(inits_list, n.chains = n_chains) 
  ## Create the JAGS model object
  parJagsModel(
    cl, name = "jags_mod", file = bugs_file, 
    data = data_list, inits = par_inits, n.chains = n_chains, n.adapt = 1000,
    quiet = T
    )
  ## Burn-in
  parUpdate(cl, "jags_mod", n.iter = 1000, progress.bar = "none")
  ## Now we run the full model samples
  mod_iters <- chain_length/n_chains
  mod_samples <- 
    parCodaSamples(
      cl, "jags_mod", variable.names = parameters,
      n.iter = mod_iters, n.chain = n_chains,
      progress.bar = "none"
      )
  
  parallel::stopCluster(cl)
  
  ##------------------------------------------------------------------------------
  ## SUMMARY AND DIAGNOSTICS
  
  # head(mod_samples)
  # summary(mod_samples)
  # quantile(mod_samples, probs = c(5, 50, 95)/100)
  #   
  # geweke.diag(mod_samples)
  # heidel.diag(mod_samples)
  # raftery.diag(mod_samples)
  
  ## Extract regression coefficients ##
  ## (NOTE: Use only those from RCP 2.6 for consistency) ##
  if (i == "rcp26") {
    
    ## Convert coefficients MCMC list into data frame for later. First combines
    ## all chains into one matrix.
    coefs_df <-
      as.matrix(mod_samples[, c(1:6)], iters = F) %>%
      # data.frame() %>% 
      as_data_frame() %>%
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
      mutate(coef = factor(coef, levels=c("beta","gamma","delta","eta","alpha","sigma"))) %>% 
      mutate(prior = paste0(prior_type, convic_type)) %>%
      arrange(coef)
    
    ## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
    tcr <- 
      data_frame(
        beta = filter(coefs_df, coef == "beta")$values,
        prior = paste0(prior_type, convic_type)
        ) 
    tcr$tcr <- tcr$beta * rf2x
    
  
    if (prior_type=="ni") {
      #############################################
      ### Figure S1: Coefficient densities plot ###
      #############################################
      fig_s1 <- coefs_plot(coefs_df)
      fig_s1 +
        ggsave(
          file = here("TablesFigures/Untracked/Robustness/PNGs/fig-s1-me.png"),
          width = 8, height = 10
          )
      fig_s1 +
        ggsave(
          file = here("TablesFigures/Untracked/Robustness/fig-s1-me.pdf"),
          width = 8, height = 10, 
          device = cairo_pdf
          )
      rm(fig_s1)
    } ## End of mini NI "if" clause
    
    rm(coefs_df)
    
  } ## End of RCP 2.6 "if" clause
  
  if (i != "rcp26") {
    coefs_tab <- NULL
    tcr <- NULL
  }
  
  ## Summarise temperature predictions over 1866-2100 ##
  
  predictions <- jagsresults(mod_samples, params = "y_pred", regex = T)
  predictions <- as_data_frame(predictions[, c("mean", "2.5%", "97.5%")])
  colnames(predictions) <- c("mean", "q025", "q975")
  
  predictions$series <- i #rcp_type
  predictions$year <- seq(from=1866, length.out=nrow(predictions))
  predictions <- 
    predictions %>% 
    gather(stat, temp, -c(year, series)) %>%
    select(year, everything())
  
  ## Full distribution of temps in 2100 by themselves 
  temp2100 <- as_data_frame(as.matrix(mod_samples[, "y_pred[235]"]))
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
    data_frame(
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
fig_4_lab <- paste0("fig-4-", prior_type, convic_type, "-me")
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


## Zoom in on historic record with comparison between model and measurement error
bind_rows(
  data_frame(
    year = climate$year[c(1:N)],
    series = "had_full",
    mean = climate$had_full[c(1:N)],
    q025 = climate$had_025[c(1:N)],
    q975 = climate$had_975[c(1:N)]
    ),
  predictions %>% 
    filter(series == "fitted") %>%
    mutate(series = as.character(series))
  ) %>%
  filter(!is.na(mean)) %>%
  mutate(series = factor(series, levels = c("had_full", "fitted"))) %>%
  ggplot(aes(x = year, y = mean, ymin = q025, ymax = q975, col = series, fill = series)) +
  geom_ribbon(lwd=0.25, alpha = .3) +
  # geom_line() +
  labs(
    title = "Comparing measurement error and model error",
    subtitle = paste0("Prior: ", match_priors(paste0(prior_type, convic_type))), 
    y = "Temperature anomaly (°C)\n"
    ) + 
  scale_colour_manual(
    values = c("had_full"="black", "fitted"="#377EB8"),
    labels = c("HadCRUT4 (measurement error)   ", "Model fit (95% CI)")
    ) +
  scale_fill_manual(
    values = c("had_full"="black", "fitted"="#377EB8"),
    labels = c("HadCRUT4 (measurement error)   ", "Model fit (95% CI)")
    ) +
  ## Historic vs Forecast period
  geom_vline(xintercept = 2005, colour = "gray35", linetype = 6) +
  theme(
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
    ) +
  ggsave(
    file = here(paste0("TablesFigures/Untracked/Robustness/PNGs/compare-", prior_type, convic_type, "-fit-me.png")),
    width = 8, height = 6
    )

## Remove data frames no longer needed
rm(N, predictions)
## Similarly, subset rcp_loop list to relevant variables for outer (prior) loop
rcp_loop <- rcp_loop[c("coefs_tab", "tcr", "temp2100")]

return(rcp_loop=rcp_loop)
