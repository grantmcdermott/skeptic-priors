## Preallocate lists for loop
predictions <- list()
temp_2100 <- list()


## Loop over all four RCPs ##
for (i in 1:4) {  
  
  rcp_type <- c("rcp26" , "rcp45" , "rcp60" , "rcp85")[i]
  
  ## Extract regression coefficients (NOTE: Use only those from RCP 2.6 for consistency) ##
  if (rcp_type == "rcp26") {
    
    clim_df <-
      climate %>%
      filter(rcp == rcp_type) %>%
      filter(!is.na(had)) %>% ## i.e. up to 2005
      select(had, trf, volc, soi, amo)
    
    ## Using the automated LearnBayes commands
    theta_sample <- blinreg(clim_df$had, 
                            cbind(alpha = 1, beta = clim_df$trf, gamma = clim_df$volc, 
                                  delta = clim_df$soi, eta = clim_df$amo), 
                            chain_length)
    
    ## Convert coefficients MCMC list into data frame for later. 
    coefs_df <-
      cbind(as.matrix(theta_sample$beta),
            sigma = theta_sample$sigma) %>%
      magrittr::set_colnames(gsub("X", "", colnames(.))) %>%
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
                           prior_type, convic_type, "-prop.pdf"),
             width = 8, height = 10, 
             device = cairo_pdf ## Need for Palatino fontspacing to work. See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466
             )
      
    ## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
    tcr[[l]] <- 
      data.frame(beta = filter(coefs_df, coef == "beta")$values,
                           prior = paste0(prior_type, convic_type)
                 ) %>%
      tbl_df()

    rm(coefs_df)
  
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
  
  ## Again, will plot the temperature prediction at 2100 for all priors and scenariors
  ## later. Uncomment the below only if specificaly want to check the results for the 
  ## noninf. priors by themselves.
#   plot(density(y_pred["2100", ]), col = "dodgerblue2", 
#        main = paste("2100", prior_type, convic_type, rcp_type), xlab = expression(~degree~C),
#        cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)

  pred <- 
    as.data.frame(
      cbind(mean = apply(y_pred, 1, mean),
            q025 = apply(y_pred, 1, quantile, probs = 0.025),
            q975 = apply(y_pred, 1, quantile, probs = 0.975)
            ),
      row.names = F
      ) 
  
  pred$series <- rcp_type
  pred$year <- seq(from = 1866, to = 2100, by = 1)
  
  pred <- 
    tbl_df(pred) %>% 
    gather(stat, temp, -c(year, series)) %>%
    select(year, everything())
  
  df_2100 <-  tbl_df(as.data.frame(y_pred["2100", ])) 
  colnames(df_2100) <- "temp"
  df_2100$rcp <- rcp_type 
  
  predictions[[i]] <- pred ## add it to the list
  temp_2100[[i]] <- df_2100 # add it to the list

  rm(list = setdiff(ls(), 
                    c("climate", 
                      "font_type",
                      "decimals", "match_coefs", "match_priors", "match_rcps",
                      "theme_coefs", "theme_pred",
                      "chain_length", "n_chains",
                      "prior_type", "convic_type", "N", "theta_sample",
                      "ptm", "l", "all_2100", "tcr",
                      "coefs_tab", "predictions", "temp_2100")
                    ))
  
  } ## END OF RCP LOOP FOR MCMC SIMLUATIONS

rm(theta_sample)

predictions <- 
  bind_rows(
    data.frame(year = climate$year[c(1:N)],
               series = "had_full",
               stat = "mean",
               temp = climate$had_full[c(1:N)]
               ),
    bind_rows(predictions)
    )

# temp_2100 <- bind_rows(temp_2100)
all_2100[[l]] <- 
  bind_rows(temp_2100) %>%
  mutate(prior = paste0(prior_type, convic_type))
rm(temp_2100)


## For later graphs ##
#     rcp_names <- c("RCP 2.6 (420 ppmv CO2)", "RCP 4.5 (540 ppmv CO2)", "RCP 6.0 (670 ppmv CO2)", "RCP 8.5 (940 ppmv CO2)")
rcp_names <- c(expression("RCP 2.6 (420 ppmv CO"[2]*")"),
              expression("RCP 4.5 (540 ppmv CO"[2]*")"),
              expression("RCP 6.0 (670 ppmv CO"[2]*")"),
              expression("RCP 8.5 (940 ppmv CO"[2]*")"))  
rcp_cols <- c("limegreen", "orchid", "orange", "red2")


# Tidy the data
## Get rid of duplicate historic (pre-2006) model fits from different RCPs 
## Simultaneously rename historic RCP 2.6 series as "fitted"
predictions <- 
  predictions %>%
  mutate(series = ifelse(year <= 2005, gsub("rcp26", "fitted", series), series)) %>%
  filter(year >= 2005 | series %in% c("had_full", "fitted")) %>%
  filter(!(year > 2005 & series %in% c("fitted"))) %>% 
  spread(stat, temp) %>% 
  arrange(series)

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
    values = c("blue", "black", "darkgreen", "darkorchid", "darkorange2", "darkred"),
    breaks = c("had_full", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"),
    labels = c("HadCRUT4", "Model fit", 
               "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", 
               "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")
    ) +
  scale_fill_manual(
    values = c("blue", NA, "lightgreen", "orchid", "orange", "red"),
    breaks = c("had_full", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"),
    labels = c("HadCRUT4", "Model fit", 
               "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", 
               "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")
    ) +
  scale_linetype_manual(
    values = c(1, 1, 2, 2, 2, 2),
    breaks = c("had_full", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"),
    labels = c("HadCRUT4", "Model fit", 
               "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", 
               "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")
    ) +
  theme_pred +
  ggsave(file = paste0("./TablesFigures/predictions-", 
                       prior_type, convic_type, "-prop.pdf"),
         width = 10, height = 6.75, 
         device = cairo_pdf) ## Need for Palatino font spacing to work. See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466

## Lastly, export the mean, historic predicted temperature series (i.e. "fitted"),
## together with the had obs, to the recursive data folder. We'll be using the
## difference between these series as noise when simulating future "true" temperatures
## in the recursive section of the paper.
y_dev <- 
  predictions %>% 
  filter(series %in% c("had_full", "fitted")) %>%
  select(year:mean) %>%
  spread(series, mean) %>%
  mutate(dev = fitted - had_full) %>%
  filter(!is.na(dev))

write_csv(y_dev, "./Recursive/Data/y-dev.csv")

rm(y_dev, predictions)
