## Preallocate lists for loop
predictions <- list()
temp_2100 <- list()


## Loop over all four RCPs ##
for (i in 1:4) {
  
  rcp_type <- c("rcp26" , "rcp45" , "rcp60" , "rcp85")[i]
  
  ## Subset the data to the relevant RCP ##
  clim_df <- 
    climate %>%
    filter(rcp == rcp_type)
  
  ##------------------------------------------------------------------------------
  ## THE BUGS/JAGS MODEL.
  
  N <- nrow(clim_df) #2100 - 1866 + 1
  
  mod_string <- paste(
    "model{
    
    for(t in 1:N) {
      mu[t] <- alpha + beta*trf[t] + gamma*volc_mean[t] + delta*soi_mean[t] + eta*amo_mean[t]
      had[t]  ~ dnorm(mu[t], tau)
      y_pred[t] ~ dnorm(mu[t], tau) ## For predictions into the future
    }
    ",
    
    if (prior_type == "ni") {
    "
    ## Noninformative prior on beta:         
    mu_beta <- 0
    sigma_beta <- 100
    "
    }else{
    "
    ## Prior on beta: One of four subjective prior, conviction combinations...
    # 1) Mod. lukewarmer: TCR ~ N(1, 0.25^2). Mean of 1 °C and uncertainty range of 1.0 °C (95% probability).
    #    Converting to beta (= TCR/3.71): B ~ N( 1/3.71, (0.25/3.71)^2 )
    # 2) Strong lukewarmer: TCR ~ N(1, 0.065^2). Mean of 1 °C and uncertainty range of 0.25 °C (95% probability).
    #    Converting to beta (= TCR/3.71): B ~ N( 1/3.71, (0.065/3.71)^2 )
    # 3) Mod. denier: TCR ~ N(0, 0.25^2). Mean of 0 °C and uncertainty range of 1.0 °C (95% probability).
    #    Converting to beta (= TCR/3.71): B ~ N( 0/3.71, (0.25/3.71)^2 )
    # 4) Strong denier: TCR ~ N(0, 0.065^2). Mean of 0 °C and uncertainty range of 0.25 °C (95% probability).
    #    Converting to beta (= TCR/3.71): B ~ N( 0/3.71, (0.065/3.71)^2 )
    "},
    
    if (prior_type == "luke") {
    "
    mu_beta <- 1/3.71"
    },
    if (prior_type == "den") {
    "
    mu_beta <- 0"
    },
    if (convic_type == "mod") {
    "
    sigma_beta <- 0.25/3.71"
    },
    if (convic_type == "strong") {
    "
    sigma_beta <- 0.065/3.71"
    }, 
    
    "

    ## Priors for all parameters   
    alpha ~ dnorm(0, 0.0001)            ## intercept
    beta ~ dnorm(mu_beta, tau_beta)     ## trf coef
      tau_beta <- pow(sigma_beta, -2)
    gamma ~ dnorm(0, 0.0001)            ## volc coef
    delta ~ dnorm(0, 0.0001)            ## soi coef
    eta ~ dnorm(0, 0.0001)              ## amo coef
    sigma ~ dunif(0, 100)               ## Residual std dev
      tau <- pow(sigma, -2)  	          
    had0 ~ dnorm(0.0, 1.0E-6)           ## Initialising value
    }" 
    ) 
  
  bugs_file <- paste0("./BUGSfiles/", 
                      prior_type, "-", convic_type, "-", rcp_type, ".txt")
  if(prior_type == "ni"){bugs_file <- gsub("--", "-", bugs_file)}
  writeLines(mod_string, con = bugs_file)
  
  load.module("lecuyer") ## JAGS module uses lecuyer random number generator (to avoid overlap/correlation in a parallel format)
  
  cl <- makeCluster(n_chains, type = "SOCK") # no. of clusters (i.e. MCMC chains), SOCK is simplest cluster
  parLoadModule(cl, "lecuyer", quiet = T)
  
  ##------------------------------------------------------------------------------
  ## INTIALIZE THE CHAINS.
  
  data_list <- list("N" = N, "had" = clim_df$had, "trf" = clim_df$trf, 
                    "volc_mean" = clim_df$volc_mean, 
                    "soi_mean" = clim_df$soi_mean, "amo_mean" = clim_df$amo_mean)
  
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
  mod_iters <- chain_length
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
  if (rcp_type == "rcp26") {
    
    ## Convert coefficients MCMC list into data frame for later. First combines
    ## all chains into one matrix.
    coefs_df <-
      as.matrix(mod_samples[, c(1:6)], iters = F) %>%
      data.frame() %>% 
      tbl_df() %>% 
      gather(coef, values)
    
    ## Get summary statistics for tables ##
    l <- l + 1
    coefs_tab[[l]] <-
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
    tcr[[l]] <- 
      data.frame(beta = filter(coefs_df, coef == "beta")$values,
                 prior = paste0(prior_type, convic_type)
                 ) %>%
      tbl_df()
    
    
    ### Density plot ###
    coefs_df %>%
      mutate(coef = gsub("alpha", "alpha[0]", coef),
             coef = gsub("beta", "beta[1]", coef),
             coef = gsub("gamma", "gamma[2]", coef),
             coef = gsub("delta", "delta[3]", coef),
             coef = ifelse(coef=="eta", "eta[4]", coef)) %>%
      mutate(coef = factor(coef, levels = c("alpha[0]", "beta[1]", "gamma[2]",
                                            "delta[3]", "eta[4]", "sigma"))
             ) %>%
      ggplot(aes(x = values, group = coef)) +
      geom_line(stat = "density") +
      theme_coefs + 
      facet_wrap(~coef, ncol = 2, scales = "free",
                 labeller = label_parsed) +
      ggsave(file = paste0("./TablesFigures/coefs-", 
                           prior_type, convic_type, ".pdf"),
             width = 8, height = 10, 
             device = cairo_pdf ## See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466
             )
    
    rm(coefs_df)
    
    } ## End of RCP 2.6 "if" clause
  
  ## Summarise temperature predictions over 1866-2100 ##
  
  pred <- jagsresults(mod_samples, params = "y_pred", regex = T)
  pred <- tbl_df(as.data.frame(pred[, c("mean", "2.5%", "97.5%")]))
  colnames(pred) <- c("mean", "q025", "q975")
  
  pred$series <- rcp_type
  pred$year <- seq(from=1866, length.out=nrow(pred))
  pred <- pred %>% 
    gather(stat, temp, -c(year, series)) %>%
    select(year, everything())
  
  ## Full distribution of temps in 2100 by themselves 
  df_2100 <- tbl_df(as.data.frame(as.matrix(mod_samples[, "y_pred[235]"])))
  colnames(df_2100) <- "temp"
  df_2100$rcp <- rcp_type 

  predictions[[i]] <- pred ## add it to the list
  temp_2100[[i]] <- df_2100 # add it to the list
  
  rm(bugs_file, cl, clim_df, data_list, df_2100, i, inits_list, mod_iters,
     mod_samples, par_inits, parameters, pred, rcp_type)
  
} ## END OF RCP LOOP FOR JAGS SIMLUATIONS


## Combine predictions from RCP loop into one data frame
predictions <- 
  bind_rows(
    data.frame(year = climate$year[c(1:N)],
               series = "had_full",
               stat = "mean",
               temp = climate$had_full[c(1:N)]
               ),
    bind_rows(predictions)
    )

## Ditto for temps in 2100
all_2100[[l]] <- 
  bind_rows(temp_2100) %>%
  mutate(prior = paste0(prior_type, convic_type))

## Tidy the data:
## Get rid of duplicate historic (pre-2006) model fits from different RCPs. 
## Simultaneously rename historic RCP 2.6 series as "fitted".
predictions <- 
  predictions %>%
  mutate(series = ifelse(year <= 2005, gsub("rcp26", "fitted", series), series)) %>%
  filter(year >= 2005 | series %in% c("had_full", "fitted")) %>%
  filter(!(year > 2005 & series %in% c("fitted"))) %>% 
  spread(stat, temp) %>% 
  arrange(series)

## Lastly, order series as factor and define readable labels for plotting
predictions <-
  predictions %>%
  mutate(series = factor(series, 
                         levels = c("had_full", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"))
         ) %>%
  arrange(series)
series_labs <- c("HadCRUT4", "Model fit", 
                 "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", 
                 "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")

## predictions plot
ggplot(data = predictions, 
       aes(x = year, col = series, fill = series, linetype = series)) +
  ylab(expression(~degree*C)) + xlab("Year") +
  geom_line(data = predictions %>% 
              filter(series %in% c("rcp26", "rcp45", "rcp60", "rcp85") ),
            aes(y = mean),
            lwd = 1) + 
  geom_ribbon(aes(ymin = q025, ymax = q975), 
              lty = 0, alpha = 0.3) +
  geom_line(data = predictions %>% 
              filter(series %in% c("had_full", "fitted")),
            aes(y = mean),
            lwd = 1) + 
  ## Historic vs Forecast period
  geom_vline(xintercept = 2005, colour = "gray50", linetype = "longdash") +
  annotate("text", x = 1985, y = -0.5, label = "Historic", size = 7, family = font_type) + 
  annotate("text", x = 2025, y = -0.5, label = "Forecast", size = 7, family = font_type) +
  scale_colour_manual(
    values = c("black", "blue", rcp_cols),
    labels = series_labs,
    limits = levels(predictions$series)
    ) +
  scale_fill_manual(
    values = c(NA, "blue", rcp_fills),
    labels = series_labs,
    limits = levels(predictions$series)
    ) +
  scale_linetype_manual(
    values = c(1, 1, 2, 2, 2, 2),
    labels = series_labs,
    limits = levels(predictions$series)
    ) +
  theme_pred +
  ggsave(file = paste("./TablesFigures/predictions-",
                      prior_type, convic_type, ".pdf", sep = ""),
         width = 10, height = 6.75,
         device = cairo_pdf) ## See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466

rm(predictions, temp_2100)
