a <- a + 1

if(recurse_type == "historic"){
  clim_df <- 
    climate %>%
    filter(year >=  n)
}

if(recurse_type == "future"){
  clim_df <- 
    climate %>%
    filter(year <=  n)
}

yr_min <- min(clim_df$year)

##------------------------------------------------------------------------------
## THE BUGS/JAGS MODEL.
## Note: Removing y_pred b/c predictions into the future not needed for recursive regs

N <- nrow(clim_df)

mod_string <- paste(
  "model{
  
  for(t in 1:N) {
  mu[t] <- alpha + beta*trf[t] + gamma*volc_sim[t] + delta*soi_sim[t] + eta*amo_sim[t]
  had_sim[t]  ~ dnorm(mu[t], tau)
  # y_pred[t] ~ dnorm(mu[t], tau) ## For predictions into the future
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

bugs_file <- paste0("./Recursive/BUGSFiles/", prior_type, "-", 
                    convic_type, ".txt")
if(prior_type == "ni"){bugs_file <- gsub("--", "-", bugs_file)}
writeLines(mod_string, con = bugs_file)

load.module("lecuyer") ## JAGS module uses lecuyer random number generator (to avoid overlap/correlation in a parallel format)

cl <- makeCluster(n_chains, type = "SOCK") # no. of clusters (i.e. MCMC chains), SOCK is simplest cluster
parLoadModule(cl, "lecuyer", quiet = T)

##------------------------------------------------------------------------------
## INTIALIZE THE CHAINS.

data_list <- list("N" = N, "had_sim" = clim_df$had_sim, "trf" = clim_df$trf, 
                  "volc_sim" = clim_df$volc_sim, 
                  "soi_sim" = clim_df$soi_sim, "amo_sim" = clim_df$amo_sim)
inits_list <- function() {
  list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1)
}

##------------------------------------------------------------------------------
## RUN THE CHAINS.

parameters <- c("alpha", "beta", "gamma", "delta", "eta", "sigma"#, "y_pred"
                )
par_inits <- parallel.inits(inits_list, n.chains = n_chains) # Initialisation

parJagsModel(cl, name = "jags_mod", file = bugs_file, 
             data = data_list, inits = par_inits, n.chains = n_chains, n.adapt = 1000)
parUpdate(cl, "jags_mod", n.iter = 1000) # burn-in
mod_iters <- chain_length
mod_samples <- parCodaSamples(cl, "jags_mod", variable.names = parameters, 
                              n.iter = mod_iters, n.chain = n_chains) 
stopCluster(cl)

##------------------------------------------------------------------------------
## Get TCR coefficient (i.e. beta) and then use to get TCR.

tcr <- as.matrix(mod_samples[, "beta"], iters = F) * rf2x

rm(bugs_file, cl, data_list, mod_samples, mod_string, par_inits, parameters)


## Set figure if haven't done so by calling the proportional noninformative prior previously
if(prior_type == "ni"){
  
  A <- A + 1 # For animation figures
  pdf(file = paste0("./Recursive/TablesFigures/Animation/", 
                    recurse_type,
                    "/rec-tcr-", 1000 + A,".pdf"), 
      width = 16, height = 10, family = "Palatino")
  par(las = 1, 
      mar = c(5.5, 5.7, 4, 2) + 0.1, 
      mgp = c(4, 1.5, 0)) # Change margin and tick labels to give space for different font
  plot(1, type = "n", 
       ylab = "Density", xlab = expression(~degree~C), 
       main = "", cex.lab = 2.5, cex.axis = 2.25, cex.main = 2.25, cex.sub = 2.25,
       xlim = c(-1, 3),
       ylim = c(0, 8),
       font.lab = 2)
  ## IPCC "likely" range
  x <- c(1.0, 2.5, 2.5, 1.0)
  y <- c(8, 8, 0, 0) 
  polygon(x, y, col = "gray90", border = NA)
  ## Noninformative
  lines(density(tcr), lwd = 4, col = prior_cols[5])
  ## Moderate lukewarmer prior
  curve(dnorm(x, mean = 1, sd = 0.25), 
        from = 0, to = 2, lwd = 4, lty = 2, col = prior_cols[4], add = T, yaxt = "n")
  ## Strong lukewarmer prior
  curve(dnorm(x, mean = 1, sd = 0.065), 
        from = 0.5, to = 1.5, lwd = 4, lty = 2, col = prior_cols[3], add = T, yaxt = "n")
  ## Moderate Denier prior
  curve(dnorm(x, mean = 0, sd = 0.25), 
        from = -1, to = 1, lwd = 4, lty = 2, col = prior_cols[2], add = T, yaxt = "n")
  ## Strong denier prior
  curve(dnorm(x, mean = 0, sd = 0.065),
        from = -0.5, to = 0.5, lwd = 4, lty = 2, col = prior_cols[1], add = T, yaxt = "n")
  ## Year tracker
  if(recurse_type == "historic"){
    text(1.75, 7, label = paste("Looking back to", yr_min), cex = 2)
    # text(1.75,6.5, label = paste("(Sample size =", yr_max - yr_min + 1, "years)"), cex = 2)
  }
  if(recurse_type == "future"){
    text(1.75, 7, label = paste("Looking forward to", max(clim_df$year)), cex = 2)
    # text(1.75,6.5, label = paste("(Sample size =", nrow(clim_df), "years)"), cex = 2)
  }
  text(1.75,6.5, label = paste("(Sample size =", nrow(clim_df), "years)"), cex = 2)
}
## Moderate lukewarmer
if(prior_type == "luke" & convic_type == "mod"){
  lines(density(tcr), lwd = 4, col = prior_cols[4])
}
## Strong lukewarmer
if(prior_type == "luke" & convic_type == "strong"){
  lines(density(tcr), lwd = 4, col = prior_cols[3])
}
## Moderate Denier
if(prior_type == "den" & convic_type == "mod"){
  lines(density(tcr), lwd = 4, col = prior_cols[2])
}
## Strong denier
if(prior_type == "den" & convic_type == "strong"){
  lines(density(tcr), lwd = 4, col = prior_cols[1])
  ## Legend
  legend("topleft", 
         prior_names, 
         col = prior_cols,
         cex = 2, lty = 1, lwd = 4, bty = "n")
  dev.off()
}

tcr_rec[[a]] <- 
  data.frame(year_to = ifelse(recurse_type == "historic", yr_min, max(clim_df$year)),
             samp_size = nrow(clim_df),
             series = paste0(prior_type, convic_type), 
             mean = mean(tcr),
             q025 = quantile(tcr, p = 0.025),
             q975 = quantile(tcr, p = 0.975),
             row.names = NULL
             )

rm(tcr, clim_df)
