## Loop over all four RCPs ##
rcp_loop <-
  
  pblapply(c("rcp26" , "rcp45" , "rcp60" , "rcp85"), function(i){
    
    clim_df <-
      climate %>%
      filter(rcp == i) %>%
      filter(!is.na(had)) %>% ## i.e. up to 2005
      select(had, trf, volc, soi, amo)
    
    ## Using the automated LearnBayes commands
    ## Assign to global environment to make available for later stages of the loop
    theta_sample <- blinreg(clim_df$had, 
                            cbind(alpha = 1, beta = clim_df$trf, gamma = clim_df$volc, 
                                  delta = clim_df$soi, eta = clim_df$amo), 
                            chain_length)
    
    if (i == "rcp26") {
    
    ## Convert coefficients MCMC list into data frame for later. 
    coefs_df <-
      cbind(as.matrix(theta_sample$beta),
            sigma = theta_sample$sigma) %>%
      magrittr::set_colnames(gsub("X", "", colnames(.))) %>%
      data.frame() %>% 
      as_data_frame() %>%
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
    
    ### Density plot ###
    coefs_plot <- coef_plot_func(coefs_df)
    coefs_plot +
      ggsave(
        file = paste0("TablesFigures/PNGs/coefs-", prior_type, convic_type, "-prop.png"),
        width = 8, height = 10
        )
    coefs_plot +
      ggsave(
        file = paste0("TablesFigures/coefs-", prior_type, convic_type, "-prop.pdf"),
        width = 8, height = 10, 
        device = cairo_pdf ## See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466
        )

    rm(coefs_df, coefs_plot)
    
    } ## End of RCP 2.6 "if" clause
    
    if (i != "rcp26") {
      coefs_tab <- NULL
      tcr <- NULL
    }
  
  
  ## Summarise temperature predictions over 1866-2100 ##
  N <- 2100 - 1866 + 1
  
  X1 <- (climate %>% 
           split(.$rcp))[[i]] %>%
    filter(year <= 2100) %>%
    select(trf, volc_mean, soi_mean, amo_mean) %>%
    as.matrix()
  X1 <- cbind(const = 1, X1)
  
  y_pred <- blinregpred(X1, theta_sample)
  colnames(y_pred) <- seq(from = 1866, to = 2100, by = 1)
  y_pred <- t(y_pred)

  predictions <-
    as_data_frame(
      cbind(mean = apply(y_pred, 1, mean),
            q025 = apply(y_pred, 1, quantile, probs = 0.025),
            q975 = apply(y_pred, 1, quantile, probs = 0.975)
            )
      ) 
  
  predictions$series <- i #rcp_type
  predictions$year <- seq(from = 1866, to = 2100, by = 1)
  
  predictions <- 
    predictions %>%
    gather(stat, temp, -c(year, series)) %>%
    select(year, everything())
  
  ## Full distribution of temps in 2100 by themselves 
  all_2100 <-  as_data_frame(y_pred["2100", ]) 
  colnames(all_2100) <- "temp"
  all_2100$rcp <- i #rcp_type 
  all_2100$prior <- paste0(prior_type, convic_type)
  
  return(list(tcr=tcr, coefs_tab=coefs_tab,
              predictions=predictions, all_2100=all_2100, 
              N=data.frame(N))) 
  
  }) ## END OF RCP LOOP FOR MCMC SIMLUATIONS

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

series_labs <- c("HadCRUT4", "Model fit", 
                 "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", 
                 "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")

## Predictions plot
pred_plot <- pred_plot_func(predictions)
pred_plot +
  ggsave(
    file = paste0("TablesFigures/PNGs/predictions-", prior_type, convic_type, ".png"),
    width = 8, height = 6
    )
pred_plot +
  ggsave(
    file = paste0("TablesFigures/predictions-", prior_type, convic_type, ".pdf"),
    width = 8, height = 6,
    device = cairo_pdf ## See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466
    )
rm(pred_plot)

## Lastly, export the mean, historic predicted temperature series (i.e. "fitted"),
## together with the had obs, to the Evidence data folder. We'll be using the
## difference between these series as noise when simulating future "true" temperatures
## in the evidence section of the paper.
y_dev <- 
  predictions %>% 
  filter(series %in% c("had_full", "fitted")) %>%
  select(year:mean) %>%
  spread(series, mean) %>%
  mutate(dev = fitted - had_full) %>%
  filter(!is.na(dev))

write_csv(y_dev, "Evidence/Data/y-dev.csv")

## Remove data frames no longer needed
rm(y_dev, N, predictions)
## Similarly, subset rcp_loop list to relevant variables for outer (prior) loop
rcp_loop <- rcp_loop[c("coefs_tab", "tcr", "all_2100")]

return(rcp_loop=rcp_loop)
