rm(list=ls()) # Clear data

## Specify data source ##
#setwd("C:/Users/Grant/Documents/NHH/Papers/Sceptic/sceptic-priors")

library(readr)
library(ggplot2)
library(cowplot) ## For cowplot ggplot theme
require(RColorBrewer)

# scc <- read.table("./PAGE09/SCC.txt"), header=T)
scc <- read_csv("./Data/PAGE09/scc.csv")
# attach(scc)

## Set colours and names for consistency among graphs ##
prior_cols <- c(brewer.pal(12, "Paired")[c(2, 4, 8, 6)], "#000000")
prior_names <- c("Strong Denier", "Moderate Denier", 
                 "Strong Lukewarmer", "Moderate Lukewarmer", 
                 "Noninformative")

decimals <- function(x, k) format(round(x, k), nsmall=k) # Decimal function (to make sure, e.g. three decimals places are always printed)

scc_all <- rbind(
  cbind(decimals(mean(scc$ni), 2), 
        decimals(quantile(scc$ni, .5), 2), 
        paste0("[", decimals(quantile(scc$ni, .025), 2), ", ", 
               decimals(quantile(scc$ni, .975), 2), "]")
        ),
  cbind(decimals(mean(scc$lukemod), 2),
        decimals(quantile(scc$lukemod, .5), 2), 
        paste0("[", decimals(quantile(scc$lukemod, .025), 2), ", ", 
               decimals(quantile(scc$lukemod, .975), 2), "]")
        ),
  cbind(decimals(mean(scc$lukestrong), 2), 
        decimals(quantile(scc$lukestrong, .5), 2), 
        paste0("[", decimals(quantile(scc$lukestrong, .025), 2), ", ", 
               decimals(quantile(scc$lukestrong, .975), 2), "]")
        ),  
  cbind(decimals(mean(scc$denmod), 2), 
        decimals(quantile(scc$denmod, .5), 2), 
        paste0("[", decimals(quantile(scc$denmod, .025), 2), ", ", 
               decimals(quantile(scc$denmod, .975), 2), "]")
        ), 
  cbind(decimals(mean(scc$denstrong), 2), 
        decimals(quantile(scc$denstrong, .5), 2), 
        paste("[", decimals(quantile(scc$denstrong, .025), 2), ", ", 
              decimals(quantile(scc$denstrong, .975), 2), "]")
        )
  )

rownames(scc_all) <- rev(prior_names)
colnames(scc_all) <- c("Mean", "Median", "95% Probability interval")
scc_all

library(stargazer) ## For nice LaTeX tables
stargazer(scc_all,
          title = "Social cost of carbon (US\\$2005)",
          label = "tab:scc",
          align = TRUE, 
          notes.align = "l",
          out = "./TablesFigures/scc.tex"
          )


## Plot data
## Note: Truncate at 99.5th percentile (density curves only), to produce reasonable looking 
## curves. Exception is the Strong Denier, since no extreme values.
pdf(file = "./TablesFigures/scc.pdf", width = 16, height = 10, family = "Palatino")
par(las = 1, mar = c(5, 6, 4, 2) + 0.5, mgp = c(4.5, 1.5, 0)) # Rotate y-axis tick labels, change margin to give space for additional y-axis text
# Noninformative
plot(density(subset(scc$ni, scc$ni <= quantile(scc$ni, .995))), 
     ylab = "Density", xlab = expression(bold("Social cost of CO"[2]*" (US$2005)")),
     main = "", cex.lab = 2.75, cex.axis = 2.5, cex.main = 2.5, cex.sub = 2.5,
     xlim = c(0, 105),
     ylim = c(0, .4),
     font.lab = 2,
     lwd = 4,
     col = prior_cols[5] 
     )
abline(v=mean(scc$ni),lty = 1, lwd = 3, col = prior_cols[5])
abline(v=quantile(scc$ni, .5), lty = 2, lwd = 3, col = prior_cols[5])
# Moderate lukewarmer
lines(density(subset(scc$lukemod, scc$lukemod <= quantile(scc$lukemod, .995))), 
      lwd = 4, col = prior_cols[4])
abline(v = mean(scc$lukemod), lty = 1, lwd = 3, col = prior_cols[4])
abline(v = quantile(scc$lukemod, prob = .5), lty = 2, lwd = 3, col = prior_cols[4])
# Strong lukewarmer
lines(density(subset(scc$lukestrong, scc$lukestrong <= quantile(scc$lukestrong, .995))), 
      lwd = 4, col = prior_cols[3])
abline(v = mean(scc$lukestrong),lty = 1, lwd = 3, col = prior_cols[3])
abline(v = quantile(scc$lukestrong, .5), lty = 2, lwd = 3, col = prior_cols[3])
# Moderate denier
lines(density(subset(scc$denmod, scc$denmod <= quantile(scc$denmod, .995))), 
      lwd = 4, col = prior_cols[2])
abline(v = mean(scc$denmod), lty = 1, lwd = 3, col = prior_cols[2])
abline(v = quantile(scc$denmod, .5), lty = 2, lwd = 3, col = prior_cols[2])
# Strong denier
lines(density(scc$denstrong), lwd = 4, col = prior_cols[1])
abline(v = mean(scc$denstrong), lty = 1, lwd = 3, col = prior_cols[1])
abline(v = quantile(scc$denstrong, .5), lty = 2, lwd = 3, col = prior_cols[1])
## Repeat nonif mean to make clear behind mod denier mean
abline(v = mean(scc$ni), lty = 1, lwd = 2, col = prior_cols[5])
## Legend
legend("topright", 
       prior_names,#[c(1:4)], 
       col = prior_cols,#[c(1:4)], 
       cex = 2.1, lty = 1, lwd = 4, bty = "n")
## Turn off device
dev.off()