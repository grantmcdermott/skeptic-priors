## Load packages ##
library(LearnBayes) ## For simulating noninformative prior (using random multivarite normal command)
library(rjags) ## For running the MCMC (Gibbs) sampler
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) ## devtools::install_github("johnbaums/jagstools") For extracting summary statistics from MCMC chain
library(grid) ## To adjust legend key width and size in ggplot2 themes that don't naturally support a grid
library(gridExtra) ## Facilitates easier labelling in ggplot2 (e.g. annote with extrafont fonts)
library(tidyverse)
library(hrbrthemes)
library(cowplot) ## For cowplot ggplot theme
library(ggthemes) ## For additional ggplot2 themes (e.g. "few") 
library(RColorBrewer)
library(extrafont) ## For additional fonts in ggplot2
library(stargazer) ## For nice LaTeX tables
library(pbapply) ## Add progress bar to *apply functions


##################################
### GLOBAL ELEMENTS AND THEMES ###
##################################

## Choose non-standard font for plots. Installation: https://github.com/wch/extrafont 
## Will revert to ggplot2 default if not available.
font_type <- choose_font(c("Fira Sans", "Palatino Linotype")[1])
suff <- ifelse(font_type=="Palatino Linotype", "-ppl", "") ## For keeping track of exported files.

## Set global plot theme 
## Note: Specific themes for various plots at the bottom of this document
theme_set(
  theme_ipsum(
    base_size = 12,
    axis_title_size = 14,
    axis_title_just = "c"
    ) +
    theme(
      text = element_text(family = font_type),
      strip.text = element_text(hjust = 0.5)
      )
  )

## Assign colours and names for later graphs ##
rcp_names <- c("(a) RCP 2.6", "(b) RCP 4.5", "(c) RCP 6.0", "(d) RCP 8.5")
rcp_cols <- c("darkgreen", "darkorchid", "darkorange2", "darkred")
rcp_fills <- c("lightgreen", "orchid", "orange", "red")

prior_names <- c("Strong Denier", "Moderate Denier", 
                "Strong Lukewarmer", "Moderate Lukewarmer", "Noninformative")
# prior_cols <- c(brewer.pal(12, "Paired")[c(2, 4, 6, 8)], "#000000")
prior_cols <- c("Strong Denier"="#1F78B4", "Moderate Denier"="#33A02C", 
                 "Strong Lukewarmer"="#E31A1C", "Moderate Lukewarmer"="#FF7F00", 
                 "Noninformative"="#000000")


#################
### FUNCTIONS ###
#################

######################################
######################################
## Global (highest) common denominator
gcd <- 
  function(x,y) {
    r <- x%%y
    return(ifelse(r, gcd(y, r), y))
  }

#######################################
#######################################
## Negate version of %in% function
"%nin%" <- Negate("%in%")

#######################################
#######################################
### Decimal function 
### (To make sure, e.g. three decimals places are always printed in tables)

decimals <- function(x, k) {
  as.double(format(round(x, k), nsmall = k))
}

#######################################
#######################################
## Match short prior names to long prior names
match_priors <- function(x) {
  x <- gsub("ni", "Noninformative", x)
  x <- gsub("lukemod", "Moderate Lukewarmer", x)
  x <- gsub("lukestrong", "Strong Lukewarmer", x)
  x <- gsub("denmod", "Moderate Denier", x)
  x <- gsub("denstrong", "Strong Denier", x)
  return(x)
}

#######################################
#######################################
## Match short coeficient names to long coeficient names
match_coefs <- function(x) {
  x <- gsub("alpha", "Constant", x)
  x <- gsub("beta", "Total radiative forcing", x)
  x <- gsub("gamma", "Volcanic aerosols", x)
  x <- gsub("delta", "SOI", x)
  x <- gsub("eta", "AMO", x)
  x <- gsub("tcr", "Implied TCR", x)
  return(x)
}


#######################################
#######################################
## Match short RCP names to long RCP names
match_rcps <- function(x) {
  x <- gsub("rcp26", rcp_names[1], x)
  x <- gsub("rcp45", rcp_names[2], x)
  x <- gsub("rcp60", rcp_names[3], x)
  x <- gsub("rcp85", rcp_names[4], x)
  return(x)
}

#######################################
#######################################
## Specific plot themes

theme_coefs <-
  theme(
    text = element_text(family = font_type),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size=14),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size=14),
    legend.position = "none",
    strip.text = element_text(size = 18, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.spacing = unit(2, "lines") ## Increase gap between facet panels
    ) 

pred_plot_func <-
  function(predictions) {
    ggplot(
      data = predictions, 
      aes(x = year, col = series, fill = series, linetype = series)
      ) +
      ylab("Temperature anomaly (°C)\n") + xlab("Year") +
      geom_line(
        data = predictions %>% filter(series %in% c("rcp26","rcp45","rcp60","rcp85")),
        aes(y = mean), lwd = 0.5
        ) + 
      geom_ribbon(
        aes(ymin = q025, ymax = q975), lty = 0, alpha = 0.3
        ) +
      geom_line(
        data = predictions %>% filter(series %in% c("had_full", "fitted")),
        aes(y = mean), lwd = 0.5
        ) + 
      ## Historic vs Forecast period
      geom_vline(xintercept = 2005, colour = "black", linetype = "longdash") +
      annotate(
        "text", x = 1985, y = max(predictions$q975, na.rm=T), 
        label = "Hindcast", size = 4.5, family = font_type
        ) + 
      annotate(
        "text", x = 2025, y = max(predictions$q975, na.rm=T), 
        label = "Forecast", size = 4.5, family = font_type
        ) +
      scale_colour_manual(
        values = c("#0066ff", "black", rcp_cols),
        labels = c("HadCRUT4 ", "Model fit"),
        breaks = c("had_full", "fitted"),
        limits = levels(predictions$series)
        ) +
      scale_fill_manual(
        values = c(NA, "black", rcp_fills),
        labels = c("HadCRUT4 ", "Model fit"),
        breaks = c("had_full", "fitted"),
        limits = levels(predictions$series)
        ) +
      scale_linetype_manual(
        values = c(1, 2, 2, 2, 2, 2),
        labels = c("HadCRUT4 ", "Model fit"),
        breaks = c("had_full", "fitted"),
        limits = levels(predictions$series)
        ) +
      ## "Fake" secondary y-axis for line labels
      scale_y_continuous(
        sec.axis = dup_axis(
          breaks = predictions %>% filter(year==2100, grepl("rcp", series)) %>% pull(mean),
          labels = gsub(" \\(forecast\\)","",series_labs[3:6]),
          name = NULL)
        ) +
      theme(
        axis.title.x = element_blank(),
        axis.text  = element_text(size=18), 
        axis.text.y.right = element_text(color=rcp_cols, margin=margin(t=0, r=0, b=0, l=-20)),
        legend.title = element_blank(), 
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.key.width = unit(2, "line")
      )  
    }

theme_2100 <-
  theme(
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank(),
    # strip.text = element_text(size = 18, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.spacing = unit(2, "lines") ## Increase gap between facet panels
  ) 

theme_tcr <-
  theme(
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank()
  )

theme_recursive <-
  theme(
    text = element_text(family = font_type),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18, angle = 0),
    axis.text  = element_text(size=17),
    legend.position = "none",
    strip.text = element_text(size = 17, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.spacing = unit(2, "lines") ## Increase gap between facet panels
  )
