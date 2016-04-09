
###################################
###################################
## Measurement error file begins ##
###################################
###################################

ggplot(climate %>% filter(year < 2006), aes(x = year)) +
  geom_line(aes(y = had_full), lwd = 1) +
  geom_ribbon(aes(ymax = had_975,
                  ymin = had_025), alpha = 0.3, fill = "blue")

## 95 bounds are not quite equal, but very close
## Also: Compare measurement error with model uncertainty i.e. sigma from noninformative model.
## Note, divide 95% ME bounds by 2 to get 1 std dev (i.e. fair comparison to model sigma).
ggplot(climate %>% 
         filter(rcp == "rcp26") %>%
         filter(!is.na(had_full)) %>%
         mutate(omega_low = (had_full - had_025)/2,
                omega_up = (had_975 - had_full)/2), 
       aes(x = year)) +
  geom_line(aes(y = omega_low), col = "blue") +
  geom_line(aes(y = omega_up), col = "red") +
  geom_hline(yintercept = 0.075) + ## Sigma mean = 0.075, sigma s.d. = 0.0045
  geom_hline(yintercept = 0.075 + 0.0045*2, lty = 2) + 
  geom_hline(yintercept = 0.075 - 0.0045*2, lty = 2) +
  labs(y = expression(paste("Measurement Error (", omega, ") vs Sigma (", sigma, ")")))


## Preallocate lists for loop
predictions <- list()
temp_2100 <- list()


## Loop over all four RCPs ##
for (i in 1:4) {
  
  rcp_type <- c("rcp26" , "rcp45" , "rcp60" , "rcp85")[i]
  
  ## Subset the data to the relevant RCP ##
  ## NEW: Add omega columns for measurement error
  clim_df <- 
    climate %>%
    filter(rcp == rcp_type) %>%
    mutate(had_omega = (had_975 - had_full)/2, ## Currently 95% bound, i.e. 2*sigma.
           cw_omega = cw_1sigma)
  
  ## Fill in measurement error for future values (needed for prediction) based on
  ## measurement error of last two decades
  omega_had_recent <- 
    (clim_df %>% 
       filter(!is.na(had_full)) %>% 
       tail(20))$had_omega %>% 
    mean()
  #0.045275
  #
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
      "mu_beta <- 1/3.71
      "},
    if (prior_type == "den") {
      "mu_beta <- 0
      "},
    if (convic_type == "mod") {
      "sigma_beta <- 0.25/3.71
      "},
    if (convic_type == "strong") {
      "sigma_beta <- 0.065/3.71
      "}, 
    
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
  
  bugs_file <- paste0("./Robustness/MeasError/BUGSfiles/", 
                      prior_type, "-", convic_type, "-", rcp_type, "-me.txt")
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
  mod_iters <- chain_length
  mod_samples <- parCodaSamples(cl, "jags_mod", variable.names = parameters, 
                                n.iter = mod_iters, n.chain = n_chains) ## n.chains was 10, but changed to 3
  stopCluster(cl)
  
  
  ##------------------------------------------------------------------------------
  ## SUMMARY AND DIAGNOSTICS
  
  #   head(mod_samples)
  #   summary(mod_samples)
  #   quantile(mod_samples, probs = c(5, 50, 95)/100)
  #   
  #   geweke.diag(mod_samples)
  #   heidel.diag(mod_samples)
  #   raftery.diag(mod_samples)
  
  ## Extract regression coefficients ##
  ## (NOTE: Use only those from RCP 2.6 for consistency) ##
  if (rcp_type == "rcp26") {
    
    ## Get parameters MCMC list into separate matrix for later. Combines all chains into one matrix.##
    coefs_mat <- as.matrix(mod_samples[, c(1:6)], iters = F)
    
    # sigma_post <- as.vector(as.matrix(mod_samples[, "sigma"], iters = F))
    # data.frame(mean = mean(sigma_post), 
    #            q025 = quantile(sigma_post, .0275), 
    #            q975 = quantile(sigma_post, .975), 
    #            row.names = NULL)
    
    
    ## Get summary statistics for tables ##
    coefs_tab <- jagsresults(mod_samples, params = c("alpha", "beta", "gamma", "delta", "eta"))
    coefs_tab <- rbind(coefs_tab[c("beta", "gamma", "delta", "eta", "alpha"), ]) 
    ### Mean and C.I. ###
    coefs_tab <- 
      cbind(decimals(coefs_tab[, 1], 3), 
            paste0("[",
                   decimals(coefs_tab[, 3], 3), 
                   ", ", 
                   decimals(coefs_tab[, 7], 3),
                   "]")
            )
    
    colnames(coefs_tab) <- c("Mean", "95% Credible Interval")
    rownames(coefs_tab) <- c("Total radiative forcing", "Stratospheric aerosols",
                             "SOI", "AMO", "Constant")
    coefs_tab
    
    
    #### Density plot ###
    coefs_df <- 
      tbl_df(as.data.frame(coefs_mat)) %>%
      select(alpha, beta, gamma, delta, eta, sigma) %>% 
      gather(coef, values) 
    
    ## Posterior TCRs, temp prediction at 2100 (and coefficient values) ##  
    tcr <- coefs_mat[, "beta"] * rf2x
    
    rm(coefs_mat)
    
    # coefs_plot <- 
    coefs_df %>%
      mutate(coef = gsub("alpha", "alpha[0]", coef),
             coef = gsub("beta", "beta[1]", coef),
             coef = gsub("gamma", "gamma[2]", coef),
             coef = gsub("delta", "delta[3]", coef),
             coef = ifelse(coef=="eta", "eta[4]", coef)) %>%
      mutate(coef = factor(coef, levels = c("alpha[0]", "beta[1]", "gamma[2]",
                                            "delta[3]", "eta[4]", "sigma"))) %>%
      ggplot(aes(x = values, group = coef)) +
      geom_line(stat = "density") +
      cowplot::theme_cowplot() +
      theme(
        text = element_text(family = font_type),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size=14),
        legend.position = "none",
        strip.text = element_text(size = 18, colour = "black"),
        strip.background = element_rect(fill = "white"), ## Facet strip
        panel.margin = unit(2, "lines") ## Increase gap between facet panels
        ) +
      facet_wrap(~coef, ncol = 2, scales = "free",
                 labeller = label_parsed) +
      ggsave(file = paste0("./Robustness/TablesFigures/coefs-",
                           prior_type, convic_type, "-me.pdf"),
             width = 8, height = 10, 
             device = cairo_pdf ## Need for Palatino fontspacing to work. See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466
             )

    ## Further, only report the noninformative (RCP 2.6) coefficient results for consistency ##
    if (prior_type == "ni") {    
      
      ### Regression table ###
      stargazer(coefs_tab, align = T, header = F, 
                title = "Posterior regression results: Noninformative priors",
                label = "tab:ni-coefs-me",
                notes.align = "l",
                notes.append = T,
                notes = c("Dependent variable: Global Mean Surface Temperature (GMST)."),
                out = "./Robustness/TablesFigures/ni-coefs-me.tex"
      )
    }    ## End of noninf. (RCP 2.6) "if" sub-clause
    
  } ## End of RCP 2.6 "if" clause
  
  ## Summarise temperature predictions over 1866-2100 ##
  
  pred <- jagsresults(mod_samples, params = "y_pred", regex = T)
  pred <- tbl_df(as.data.frame(pred[, c("mean", "2.5%", "97.5%")]))
  colnames(pred) <- c("mean", "q025", "q975")
  
  pred$series <- rcp_type
  # pred$prior <- prior_type
  # pred$conviction <- convic_type
  pred$year <- seq(from=1866, length.out=nrow(pred))
  pred <- pred %>% 
    gather(stat, temp, -c(year, series)) %>%
    select(year, everything())
  
  df_2100 <- tbl_df(as.data.frame(as.matrix(mod_samples[, "y_pred[235]"])))
  colnames(df_2100) <- "temp"
  df_2100$rcp <- rcp_type 
  
  rm(mod_samples)
  
  predictions[[i]] <- pred ## add it to the list
  # predictions$rcp_type <- pred ## add it to the list
  temp_2100[[i]] <- df_2100 # add it to the list
  # temp_2100$rcp_type <- df_2100 # add it to the list
  
  rm(list = setdiff(ls(), 
                    c("climate", "rf2x", "font_type",
                      "match_coefs", "match_priors", "decimals",
                      "chain_length", "n_chains",
                      "prior_type", "convic_type", "N",
                      "coefs_tab", "predictions", "temp_2100",
                      "ni_2100", "luke_mod_2100", "luke_strong_2100", "den_mod_2100", "den_strong_2100",
                      "tcr", "tcr_ni", "tcr_luke_mod", "tcr_luke_strong", "tcr_den_mod", "tcr_den_strong",
                      "coefs_ni", "coefs_luke_mod", "coefs_luke_strong", "coefs_den_mod", "coefs_den_strong"
                    )))
  } ## END OF RCP LOOP FOR JAGS SIMLUATIONS


## For later graphs ##
rcp_names <- c(expression("RCP 2.6 (420 ppmv CO"[2]*")"),
               expression("RCP 4.5 (540 ppmv CO"[2]*")"),
               expression("RCP 6.0 (670 ppmv CO"[2]*")"),
               expression("RCP 8.5 (940 ppmv CO"[2]*")"))  
rcp_cols <- c("limegreen", "orchid", "orange", "red2")


predictions <- 
  bind_rows(
    data.frame(year = climate$year[c(1:N)],
               series = "had_full",
               stat = "mean",
               temp = climate$had_full[c(1:N)]
               ),
    bind_rows(predictions)
  )

temp_2100 <- bind_rows(temp_2100)

## Temperatures only in the year 2100
ggplot(temp_2100 %>% group_by(rcp),
       aes(x = temp)) +
  cowplot::theme_cowplot() +
  theme(text = element_text(family = font_type),
        legend.title = element_blank()
        ) +
  geom_density(aes(fill = rcp, linetype = NA), alpha = 0.5) +
  geom_density(aes(col = rcp)) +
  # geom_line(stat = "density", aes(col = rcp)) +
  labs(x = expression(~degree*C), y = "Density", 
       title = paste("Temp in 2100:", prior_type, convic_type)) +
  scale_colour_brewer(palette = "Spectral", breaks = levels(as.factor(temp_2100$rcp)),
                      labels = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5")
                      ) +
  scale_fill_brewer(palette = "Spectral", breaks = levels(as.factor(temp_2100$rcp)),
                    labels = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5")
                    )



if (prior_type == "ni") {
  ni_2100 <- temp_2100
  tcr_ni <- tcr
  coefs_ni <- coefs_tab
}

if (prior_type == "luke" & convic_type == "mod") {
  luke_mod_2100 <- temp_2100
  tcr_luke_mod <- tcr
  coefs_luke_mod <- coefs_tab
}

if (prior_type == "luke" & convic_type == "strong") {
  luke_strong_2100 <- temp_2100
  tcr_luke_strong <- tcr
  coefs_luke_strong <- coefs_tab
}

if (prior_type == "den" & convic_type == "mod") {
  den_mod_2100 <- temp_2100
  tcr_den_mod <- tcr
  coefs_den_mod <- coefs_tab
}

if (prior_type == "den" & convic_type == "strong") {
  den_strong_2100 <- temp_2100
  tcr_den_strong <- tcr
  coefs_den_strong <- coefs_tab
}

# rm(temp_2100, tcr, coefs_tab)

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

## predictions plot
ggplot(data = predictions, 
       aes(x = year, col = series, fill = series, linetype = series)) +
  #ggthemes::theme_few() + 
  cowplot::theme_cowplot() +
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
  theme(
    axis.line.x = element_line(linetype = 1), ## Temporary bug(?) in cowplot theme: missing axis line
    axis.line.y = element_line(linetype = 1), ## Ditto
    text = element_text(family = font_type),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20, angle = 0),
    axis.text  = element_text(size=18),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey90", size = 1),
    # legend.position = c(.16, .81),
    legend.position = c(.2, .75),
    legend.title = element_blank(), # switch off the legend title
    legend.text = element_text(size=18),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.key.width = unit(3.75, "line"),
    legend.key.height = unit(2.25, "line"),
    legend.key.size = unit(2, "line")
    ) +
    ggsave(file = paste("./Robustness/TablesFigures/predictions-", 
                        prior_type, convic_type, "-me.pdf", sep = ""),
           width = 10, height = 6.75,
           device = cairo_pdf) ## Need for Palatino font spacing to work. See: https://github.com/wch/extrafont/issues/8#issuecomment-50245466


## Zoom in on historic record with comparison between model and measurement error
bind_rows(
  data.frame(year = climate$year[c(1:N)],
             series = "had_full",
             mean = climate$had_full[c(1:N)],
             q025 = climate$had_025[c(1:N)],
             q975 = climate$had_975[c(1:N)]
             ),
  predictions %>% filter(series == "fitted")
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
  theme(
    text = element_text(family = font_type),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20, angle = 0),
    axis.text  = element_text(size=18),
    legend.position = c(.3, .85),
    legend.title = element_blank(), # switch off the legend title
    legend.text = element_text(size=18),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.key.width = unit(3.75, "line"),
    legend.key.height = unit(2.25, "line"),
    legend.key.size = unit(2, "line")
    ) +
  ggsave(file = paste("./Robustness/TablesFigures/fit-measerror-compare-", 
                      prior_type, convic_type, ".pdf", sep = ""),
         width = 10, height = 6.75,
         device = cairo_pdf)
