## Load packages ##
library(readr) ## For reading in data files
library(LearnBayes) ## For simulating noninformative prior (using random multivarite normal command)
library(rjags) ## For running the MCMC (Gibbs) sampler
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) ## devtools::install_github("johnbaums/jagstools") For extracting summary statistics from MCMC chain
library(ggplot2)
library(cowplot) ## For cowplot ggplot theme
library(ggthemes) ## For additional ggplot2 themes (e.g. "few") 
library(RColorBrewer)
library(grid) ## To adjust legend key width and size in ggplot2 themes that don't naturally support a grid
library(gridExtra) ## Facilitates easier labelling in ggplot2 (e.g. annote with extrafont fonts)
library(extrafont) ## For additional fonts in ggplot2
library(stargazer) ## For nice LaTeX tables
library(dplyr) ## For manipulating and munging data frames
library(tidyr) ## For tidying data frames
library(purrr) ## For manipulating vectors and functions (complements dplyr)
library(pbapply) ## Add progress bar to apply functions


#######################################
#######################################
## Choose font type for graphs (note extrafont package installation instructions)
font_type <- choose_font("Palatino Linotype") ## Will use ggplot2 default font if not available

## Assign colours and names for later graphs ##
rcp_names <- c(expression("RCP 2.6 (420 ppmv CO"[2]*")"),
               expression("RCP 4.5 (540 ppmv CO"[2]*")"),
               expression("RCP 6.0 (670 ppmv CO"[2]*")"),
               expression("RCP 8.5 (940 ppmv CO"[2]*")"))  
# rcp_cols <- c("limegreen", "orchid", "orange", "red2")
rcp_cols <- c("darkgreen", "darkorchid", "darkorange2", "darkred")
rcp_fills <- c("lightgreen", "orchid", "orange", "red")

# prior_cols <- c("dodgerblue2", "limegreen", "orange", "red2", "gray20")
prior_cols <- c(brewer.pal(12, "Paired")[c(2, 4, 8, 6)], "#000000")
prior_names <- c("Strong Denier", "Moderate Denier", 
                 "Strong Lukewarmer", "Moderate Lukewarmer", 
                 "Noninformative")

#######################################
#######################################

## Negate version of %in% function
"%nin%" <- Negate("%in%")

#######################################
#######################################

### Decimal function (to make sure, e.g. three decimals places are 
### always printed in tables) ###

decimals <- function(x, k) {
  as.double(format(round(x, k), nsmall = k))
}

#######################################
#######################################

### Match short prior names to long prior names
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

### Match short coeficient names to long coeficient names
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

### Match short RCP names to long RCP names
match_rcps <- function(x) {
  x <- gsub("rcp26", rcp_names[1], x)
  x <- gsub("rcp45", rcp_names[2], x)
  x <- gsub("rcp60", rcp_names[3], x)
  x <- gsub("rcp85", rcp_names[4], x)
  return(x)
}

#######################################
#######################################

### GGPLOT2 themes

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
    panel.margin = unit(2, "lines") ## Increase gap between facet panels
    ) 

theme_pred <-
  theme(
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
  ) 

theme_2100 <-
  theme(
    text = element_text(family = font_type),
    legend.position = "bottom",
    legend.title = element_blank(),
    # strip.text = element_text(size = 18, colour = "black"),
    strip.background = element_rect(fill = "white"), ## Facet strip
    panel.margin = unit(2, "lines") ## Increase gap between facet panels
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
    panel.margin = unit(2, "lines") ## Increase gap between facet panels
  )