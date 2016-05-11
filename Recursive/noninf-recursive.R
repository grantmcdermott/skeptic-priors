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
## GET PARAMETER SAMPLES 

theta_sample <- 
  blinreg(clim_df$had_sim, 
          cbind(1, clim_df$trf, clim_df$volc_sim, clim_df$soi_sim, clim_df$amo_sim), 
          chain_length * n_chains
          )

## Get coefficients MCMC list into separate matrix for later. Combines all chains into one matrix.##
coefs_mat <- as.matrix(theta_sample$beta)

tcr <- coefs_mat[, 2] * rf2x

rm(theta_sample, coefs_mat)


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
## Moderate lukewarmer
curve(dnorm(x, mean = 1, sd = 0.25), 
      from = 0, to = 2, lwd = 4, lty = 2, col = prior_cols[4], add = T, yaxt = "n")
## Strong lukewarmer
curve(dnorm(x, mean = 1, sd = 0.065), 
      from = 0.5, to = 1.5, lwd = 4, lty = 2, col = prior_cols[3], add = T, yaxt = "n")
## Moderate Denier
curve(dnorm(x, mean = 0, sd = 0.25), 
      from = -1, to = 1, lwd = 4, lty = 2, col = prior_cols[2], add = T, yaxt = "n")
## Strong denier
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
