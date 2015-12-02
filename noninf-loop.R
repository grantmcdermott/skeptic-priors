## Preallocate lists for loop
predictions <- list()
temp_2100 <- list()


## Loop over all four RCPs ##
for (i in 1:4) {  
  
  rcp_type <- c("rcp26" , "rcp45" , "rcp60" , "rcp85")[i]
  
  ## Extract regression coefficients (NOTE: Use only those from RCP 2.6 for consistency) ##
  if (rcp_type == "rcp26") {
    
    ## Using the automated LearnBayes commands
    theta_sample <- blinreg(climate$had[c(1:140)], 
                            cbind(1, climate$trf, climate$volc, climate$soi, climate$amo)[1:140, ], 
                            chain_length
                            )
    
    ## Get coefficients MCMC list into separate matrix for later. Combines all chains into one matrix.##
    coefs_mat <- as.matrix(theta_sample$beta)
    ## Get summary statistics for tables ##
    coefs_tab <- t(rbind(apply(coefs_mat, 2, mean),
                         apply(coefs_mat, 2, quantile, probs = c(0.025, 0.975)) 
                         ) 
                   )
    # rownames(coefs_tab) <- c("alpha", "beta", "gamma", "delta", "eta")
    coefs_tab <- rbind(coefs_tab["X2",], coefs_tab["X3",], coefs_tab["X4",],
                       coefs_tab["X5",], coefs_tab["X1",]) # i.e order as: beta, gamma, delta, eta, alpha
    coefs_tab <- cbind(decimals(coefs_tab[, 1], 3), 
                       paste("[",
                             decimals(coefs_tab[, 2], 3), 
                             ", ", 
                             decimals(coefs_tab[, 3], 3),
                             "]",
                             sep = ""))
    colnames(coefs_tab) <- c("Mean", "95% Credible Interval")
    rownames(coefs_tab) <- c("Total radiative forcing", "Stratospheric aerosols",
                            "SOI", "AMO", "Constant")
    
    ## The regression results for all prior types will be exported later. Uncomment 
    ## the below if specifically want the noninf. coefs by themselves for some reason.
#     stargazer(coefs_tab, align = T, header = F, 
#               title = "Posterior regression results: Noninformative priors",
#               label = "tab:ni-coefs",
#               out =  paste0("./TablesFigures/", prior_type, convic_type, "-coefs-prop.tex")
#               )
    
    ### Density plot ###
    coefs_df <- tbl_df(as.data.frame(cbind(coefs_mat, theta_sample$sigma)))
    colnames(coefs_df) <- c("alpha", "beta", "gamma", "delta", "eta", "sigma")
    coefs_df <- coefs_df %>% gather(coef, values)
    
#     coefs_df$coef <- factor(coefs_df$coef,
#                              levels = c("alpha", "beta", "gamma", "delta", "eta", "sigma")
#                             )
    coefs_plot <- 
      ggplot(coefs_df, aes(x = values)) + 
      cowplot::theme_cowplot() +
      theme(
        text = element_text(family = "Palatino Linotype"),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=14),
        legend.position = "none",
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "white"), ## Facet strip
        panel.margin = unit(2, "lines") ## Increase gap between facet panels
        ) +
      geom_line(aes(group = coef), stat = "density")  +
      facet_wrap( ~ coef, ncol = 2, scales = "free") 
    
      ggsave(file = paste0("./TablesFigures/coefs-", prior_type, convic_type, "-prop.pdf"),
             plot = facet_wrap_labeller(coefs_plot, 
                                        # labeller = label_parsed ## Works, but want to add no. subscripts
                                        expression(alpha[0],beta[1],gamma[2],delta[3],eta[4],sigma)
                                        ),
             width = 8, height = 10, 
             device = cairo_pdf ## Need for Palatino font spacing to work. See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466
             )
      
      rm(coefs_plot)
    
      
    ## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
    tcr <- coefs_mat[, 2] * rf2x
    
    tcr_ni <- tcr
    coefs_ni <- coefs_tab
    
    rm(coefs_mat)
  
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
  ## later. However, uncomment if specificaly want to check the results for the noninf.
  ## priors by themselves.
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
  
  pred <- tbl_df(pred) %>% 
    gather(stat, temp, -c(year, series)) %>%
    select(year, everything())
  
  # df_2100 <-  tbl_df(as.data.frame(y_pred[, N])) 
  df_2100 <-  tbl_df(as.data.frame(y_pred["2100", ])) 
  colnames(df_2100) <- "temp"
  df_2100$rcp <- rcp_type 
  
  predictions[[i]] <- pred ## add it to the list
  # predictions$rcp_type <- pred ## add it to the list
  temp_2100[[i]] <- df_2100 # add it to the list
  # temp_2100$rcp_type <- df_2100 # add it to the list

  rm(list = setdiff(ls(), 
                    c("climate", "rf2x",
                      "facet_wrap_labeller", "print.arrange", "decimals",
                      "chain_length", "n_chains",
                      "prior_type", "convic_type", "N", "theta_sample",
                      "coefs_tab", "predictions", "temp_2100",
                      "ni_2100", "luke_mod_2100", "luke_strong_2100", "den_mod_2100", "den_strong_2100",
                      "tcr_ni", "tcr_luke_mod", "tcr_luke_strong", "tcr_den_mod", "tcr_den_strong",
                      "coefs_ni", "coefs_luke_mod", "coefs_luke_strong", "coefs_den_mod", "coefs_den_strong"
                      )))
  
  } ## END OF RCP LOOP FOR MCMC SIMLUATIONS

rm(theta_sample)

predictions <- 
  bind_rows(
    data.frame(year = climate$year[c(1:N)],
               series = "had",
               stat = "mean",
               temp = climate$had[c(1:N)]
               ),
    dplyr::bind_rows(predictions)
    )

temp_2100 <- dplyr::bind_rows(temp_2100)



## For later graphs ##
#     rcp_names <- c("RCP 2.6 (420 ppmv CO2)", "RCP 4.5 (540 ppmv CO2)", "RCP 6.0 (670 ppmv CO2)", "RCP 8.5 (940 ppmv CO2)")
rcp_names <- c(expression("RCP 2.6 (420 ppmv CO"[2]*")"),
              expression("RCP 4.5 (540 ppmv CO"[2]*")"),
              expression("RCP 6.0 (670 ppmv CO"[2]*")"),
              expression("RCP 8.5 (940 ppmv CO"[2]*")"))  
rcp_cols <- c("limegreen", "orchid", "orange", "red2")

ni_2100 <- temp_2100

# Tidy the data
## Get rid of duplicate historic (pre-2006) model fits from different RCPs 
## Simultaneously rename historic RCP 2.6 series as "fitted"
predictions <- predictions %>%
  mutate(series = ifelse(year <= 2005, gsub("rcp26", "fitted", series), series)) %>%
  filter(year >= 2005 | series %in% c("had", "fitted")) %>%
  filter(!(year > 2005 & series %in% c("fitted"))) %>% 
  spread(stat, temp) %>% 
  arrange(series)

predictions_plot <- 
  ggplot(data = predictions, aes(x = year, col = series, linetype = series)) +
  theme_few() + # use few theme rather than grey (do this before font changes, or it overrides them)
  ylab(expression(~degree~C)) + xlab("Year") +
  geom_line(data = predictions %>% 
              filter(series %in% c("rcp26", "rcp45", "rcp60", "rcp85") ),
            aes(y = mean),
            lwd = 1
            ) + 
  geom_ribbon(data = predictions,
              aes(ymin = q025, ymax = q975, fill = series), lty = 0, alpha = 0.3
              ) +
  geom_line(data = predictions %>% 
              filter(series %in% c("had", "fitted")),
            aes(y = mean),
            lwd = 1
            ) + 
  ## Historic vs Forecas period
  geom_vline(xintercept = 2005, colour = "gray50", linetype = "longdash") +
  annotate("text", x = 1985, y = -0.5, label = "Historic", size = 7, family = "Palatino Linotype") + 
  annotate("text", x = 2025, y = -0.5, label = "Forecast", size = 7, family = "Palatino Linotype") +
  scale_colour_manual(
    values = c("blue", "black", "darkgreen", "darkorchid", "darkorange2", "darkred"),
    breaks = c("had", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"),
    labels = c("HadCRUT4", "Model fit", 
               "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")
    ) +
  scale_fill_manual(
    values = c("blue", NA, "lightgreen", "orchid", "orange", "red"),
    breaks = c("had", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"),
    labels = c("HadCRUT4", "Model fit", 
               "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")
    ) +
  scale_linetype_manual(
    values = c(1, 1, 2, 2, 2, 2),
    breaks = c("had", "fitted", "rcp26", "rcp45", "rcp60", "rcp85"),
    labels = c("HadCRUT4", "Model fit", 
               "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")
    ) 

predictions_plot +
  theme(
    text = element_text(family = "Palatino Linotype"),
    axis.title.x = element_text(face="bold", size=20),
    axis.title.y = element_text(face="bold", size=20, angle = 0),
    axis.text  = element_text(size=18),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", size = 1),
    # legend.position = c(.16, .81),
    legend.position = c(.2, .75),
    legend.title = element_blank(), # switch off the legend title
    legend.text = element_text(size=18),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.key.width = unit(3.75, "line"),
    legend.key.height = unit(2.25, "line"),
    legend.key.size = unit(2, "line")
  ) +
  ggsave(file = paste0("./TablesFigures/predictions-", prior_type, convic_type, "-prop.pdf"),
         width = 10, height = 6.75, 
         device = cairo_pdf) ## Need for Palatino font spacing to work. See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466
