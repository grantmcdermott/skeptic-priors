## Loop over all four RCPs ##
rcp_loop <-
  
  pblapply(c("rcp26" , "rcp45" , "rcp60" , "rcp85"), function(i){
  
  ## Subset the data to the relevant RCP ##
  ## NEW: Add omega columns for measurement error
  clim_df <- 
    climate %>%
    filter(rcp == i) %>%
    mutate(had_omega = (had_975 - had_full)/2, ## Currently 95% bound, i.e. 2*sigma.
           cw_omega = cw_1sigma)
  
  ## Fill in measurement error for future values (needed for prediction) based on
  ## measurement error of last two decades
  omega_had_recent <- 
    (clim_df %>% 
       filter(!is.na(had_full)) %>% 
       tail(20))$had_omega %>% 
    mean()

  clim_df$had_omega <- ifelse(is.na(clim_df$had_omega), omega_had_recent, clim_df$had_omega)
  
  
  ##------------------------------------------------------------------------------
  ## THE BUGS/JAGS MODEL.
  
  N <- nrow(clim_df)
  
  mod_string <- paste(
    "model{
    
    for(t in 1:N) {
    mu[t] <- alpha + beta*trf[t] + gamma*volc_mean[t] + delta*soi_mean[t] + eta*amo_mean[t]
    
    ## Meas. error in variance of y, i.e.: y[t]* = y[t] + v[t], where Var(v[t]) = omega[t]
    ## Note change in variance and tau formula because of combined normal distributions
    ## See: https://en.wikipedia.org/wiki/Sum_of_normally_distributed_random_variables
    
    sigmasq_tot[t] <- pow(sigma, 2) + pow(had_omega[t], 2)
    tau_tot[t] <- pow(sigmasq_tot[t], -1)
    
    had[t]  ~ dnorm(mu[t], tau_tot[t])
    y_pred[t] ~ dnorm(mu[t], tau_tot[t]) 
    }
    ", 
    paste0(
      "
            mu_beta <- ", mu_beta),
    paste0(
      "
            sigma_beta <- ", sigma_beta),
    "

    ## Priors for all parameters   
    alpha ~ dnorm(0, 0.0001)            ## intercept
    tau_beta <- pow(sigma_beta, -2)
    beta ~ dnorm(mu_beta, tau_beta)     ## trf coef
    gamma ~ dnorm(0, 0.0001)            ## volc coef
    delta ~ dnorm(0, 0.0001)            ## soi coef
    eta ~ dnorm(0, 0.0001)              ## amo coef
    sigma ~ dunif(0, 100)               ## Residual std dev
    #tau <- pow(sigma, -2) ## Changed to commented out
    had0 ~ dnorm(0.0, 1.0E-6)           ## Initialising value
    }" 
    ) 
  
  bugs_file <- paste0("./Robustness/MeasError/BUGSFiles/", 
                      prior_type, "-", convic_type, "-", i, #rcp_type, 
                      "-me.txt")
  if(prior_type == "ni"){bugs_file <- gsub("--","-",bugs_file)}
  writeLines(mod_string, con = bugs_file)
  
  load.module("lecuyer") ## JAGS module uses lecuyer random number generator (to avoid overlap/correlation in a parallel format)
  
  cl <- makeCluster(n_chains, type = "SOCK") # no. of clusters (i.e. MCMC chains), SOCK is simplest cluster
  parLoadModule(cl, "lecuyer", quiet = T)
  
  ##------------------------------------------------------------------------------
  ## INTIALIZE THE CHAINS.
  
  data_list <- list("N" = N, "had" = clim_df$had, "trf" = clim_df$trf, 
                    "volc_mean" = clim_df$volc_mean, 
                    "soi_mean" = clim_df$soi_mean, "amo_mean" = clim_df$amo_mean,
                    "had_omega" = clim_df$had_omega ## variance of ME
                    ) 
  
  inits_list <- function() {
    list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1)
  }
  
  ##------------------------------------------------------------------------------
  ## RUN THE CHAINS.
  
  parameters <- c("alpha", "beta", "gamma", "delta", "eta", "sigma", "y_pred")
  par_inits <- parallel.inits(inits_list, n.chains = n_chains) # Initialisation
  
  parJagsModel(cl, name = "jags_mod", file = bugs_file, 
               data = data_list, inits = par_inits, n.chains = n_chains, n.adapt = 1000)
  parUpdate(cl, "jags_mod", n.iter = 1000) # burn-in
  mod_iters <- chain_length/n_chains
  mod_samples <- parCodaSamples(cl, "jags_mod", variable.names = parameters, 
                                n.iter = mod_iters, n.chain = n_chains)
  stopCluster(cl)
  
  
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
      data.frame() %>% 
      tbl_df() %>% 
      gather(coef, values)
    
    ## Get summary statistics for tables ##
    coefs_tab <-
      coefs_df %>% 
      group_by(coef) %>% 
      summarise(mean = mean(values), 
                q025 = quantile(values, 0.025), 
                q975 = quantile(values, 0.975)) %>% 
      mutate(coef = factor(coef, levels = c("beta","gamma","delta",
                                            "eta","alpha","sigma"))
             ) %>% 
      mutate(prior = paste0(prior_type, convic_type)) %>%
      arrange(coef)
    
    ## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
    tcr <- 
      data_frame(beta = filter(coefs_df, coef == "beta")$values,
                 prior = paste0(prior_type, convic_type)
                 ) 
    tcr$tcr <- tcr$beta * rf2x
    
  
    if (prior_type=="ni") {
      #############################################
      ### Figure S1: Coefficient densities plot ###
      #############################################
      fig_s1 <- coef_plot(coefs_df)
      fig_s1 +
        ggsave(
          file = paste0("TablesFigures/Untracked/PNGs/fig-s1.png"),
          width = 8, height = 10
          )
      fig_s1 +
        ggsave(
          file = paste0("TablesFigures/Untracked/fig-s1.pdf"),
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
  # predictions <- tbl_df(as.data.frame(predictions[, c("mean", "2.5%", "97.5%")]))
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
    data_frame(year = climate$year[1:N],
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
  mutate(series = factor(series, 
                         levels = c("had_full", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"))
         ) 


##########################################
### Figure 4: Model fit and prediction ###
##########################################

fig_4 <- pred_plot(predictions)
fig_4_dir <- "TablesFigures/Untracked/Robustness/"
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
  data.frame(year = climate$year[c(1:N)],
             series = "had_full",
             mean = climate$had_full[c(1:N)],
             q025 = climate$had_025[c(1:N)],
             q975 = climate$had_975[c(1:N)]
             ),
  predictions %>% 
    filter(series == "fitted")
  ) %>%
  filter(!is.na(mean)) %>%
  mutate(series = factor(series, levels = c("had_full", "fitted"))) %>%
  ggplot(aes(x = year, col = series)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = series), alpha = .3) +
  geom_line(aes(y = mean), lwd = 1) +
  labs(x = "Year", y = expression(~degree*C)) +
  scale_colour_manual(values = c("black", "blue"),
                      labels = c("HadCRUT4 (plus measurement error)", 
                                 "Model fit (plus 95% CI)")) +
  scale_fill_manual(values = c("grey40", "blue"),
                    labels = c("HadCRUT4 (plus measurement error)", 
                               "Model fit (plus 95% CI)")) +
  geom_vline(xintercept = 2005, linetype = 4) +
  theme_pred +
  theme(legend.position = c(.325, .9)) +
  ggsave(
    file = paste0("TablesFigures/Untracked/Robustness/PNGS/compare-", prior_type, convic_type, "-fit-me.png"),
    width = 10, height = 6.75
    )

## Remove data frames no longer needed
rm(N, predictions)
## Similarly, subset rcp_loop list to relevant variables for outer (prior) loop
rcp_loop <- rcp_loop[c("coefs_tab", "tcr", "temp2100")]

return(rcp_loop=rcp_loop)
