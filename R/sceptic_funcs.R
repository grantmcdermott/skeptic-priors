## Load packages ##
library(LearnBayes) ## For simulating noninformative prior (using random multivarite normal command)
library(rjags) ## For running the MCMC (Gibbs) sampler
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) ## devtools::install_github("johnbaums/jagstools") Extract summary statistics from MCMC objects
library(grid) ## To adjust legend key width and size in ggplot2 themes that don't naturally support a grid
library(gridExtra) ## Easier labelling in ggplot2 (e.g. annote with extrafont fonts)
library(tidyverse)
library(hrbrthemes) ## For the theme_ipsum() ggplot2 theme
# library(cowplot) ## Added ggplot2 functionality and theme_cowplot() theme
library(ggridges) ## For ridge plots
library(RColorBrewer) ## Nice colour palettes
library(extrafont) ## For additional fonts in ggplot2
library(stargazer) ## For nice LaTeX tables
library(xtable) ## Another LaTeX table option
library(pbapply) ## Add progress bar to *apply functions


##################################
### GLOBAL ELEMENTS AND THEMES ###
##################################

## Choose non-standard font for plots. Installation: https://github.com/wch/extrafont 
## Will revert to ggplot2 default if not available.
font_type <- choose_font(c("Fira Sans", "Open Sans"))

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
rcp_cols <- scales::viridis_pal(option="plasma")(9)[c(1,3,5,7)] 

prior_names <- c("Strong Denier", "Moderate Denier", 
                "Strong Lukewarmer", "Moderate Lukewarmer", "Noninformative")
# c(brewer.pal(12, "Paired")[c(2, 1, 6, 5)], "#000000") ## Want slightly darker for light pairs
prior_cols <- c("Strong Denier"="#1F78B4", "Moderate Denier"="#8BBDDA",
                "Strong Lukewarmer"="#E31A1C", "Moderate Lukewarmer"="#F68080",
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


######################################
######################################
####    * PLOTTING FUNCTIONS *    ####
######################################
######################################

######################################
######################################
## Coeficients densities plot function

coefs_plot <-
  function(coefs_df) {
    coefs_df %>%
      mutate(coef = factor(coef, levels=c("alpha","beta","gamma", "delta","eta","sigma"))) %>%
      ggplot(aes(x = values, group = coef)) +
      geom_density(alpha=0.2, fill="black") +
      facet_wrap(~coef, ncol = 2, scales = "free", labeller = label_parsed) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )
  }


################################################
################################################
## Temperature prediction plot function (by RCP)

pred_plot <-
  function(predictions) {
    
    # series_labs <- c("HadCRUT4", "Model fit", 
    #                  "RCP 2.6 (forecast)", "RCP 4.5 (forecast)", 
    #                  "RCP 6.0 (forecast)", "RCP 8.5 (forecast)")
    series_labs <- trimws(gsub(".*\\)","",rcp_names))
    
    ggplot(
      data = predictions, 
      aes(x = year, col = series, fill = series, linetype = series)
      ) +
      ylab("Temperature anomaly (°C)\n") + xlab("Year") +
      geom_line(
        data = predictions %>% filter(series %in% c("fitted","rcp26","rcp45","rcp60","rcp85")),
        aes(y = mean), lwd = 0.5
        ) + 
      geom_ribbon(
        aes(ymin = q025, ymax = q975), lty = 0, alpha = 0.3
        ) +
      geom_line(
        data = predictions %>% filter(series %in% c("had_full")),
        aes(y = mean), lwd = 0.5
        ) + 
      ## Historic vs Forecast period
      geom_vline(xintercept = 2005, colour = "gray35", linetype = 6) +
      annotate(
        "text", x = 1985, y = max(predictions$q975, na.rm=T), 
        label = "Hindcast", size = 4.5, family = font_type, colour = "gray35"
        ) + 
      annotate(
        "text", x = 2025, y = max(predictions$q975, na.rm=T), 
        label = "Forecast", size = 4.5, family = font_type, colour = "gray35"
        ) +
      scale_colour_manual(
        values = c("black", "#377EB8", rcp_cols),
        labels = c("HadCRUT4 ", "Model fit"),
        breaks = c("had_full", "fitted"),
        limits = levels(predictions$series)
        ) +
      scale_fill_manual(
        values = c(NA, "#377EB8", rcp_cols),
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
          labels = gsub(" \\(forecast\\)","",series_labs),
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


############################
############################
## TCR density plot function

tcr_plot <-
  function(tcr) {
    
    ## Can't use stat_function() that maps to facets, ridges or other aesthetic elements.
    ## So have to create the data manually instead.
    ## See: https://github.com/tidyverse/ggplot2/issues/2357
    p_df <- 
      priors_df %>%
      mutate(prior = paste0(prior_type, convic_type)) %>%
      mutate(prior = factor(match_priors(prior), levels=prior_names)) %>% 
      filter(prior_type!="ni")
    prior_dens <- 
      lapply(prior_names[1:4], function(x){
        df <- p_df %>% filter(prior==x)
        tcr_grid <- seq(from=qnorm(0.0001,df$mu,df$sigma), to=qnorm(0.9999,df$mu,df$sigma), length=100)
        data_frame(
          tcr = tcr_grid,
          height = dnorm(tcr_grid, mean=df$mu, sd=df$sigma),
          prior = x
        )
      }) %>%
      bind_rows()
    
    tcr_density <- 
      lapply(unique(tcr$prior), function(x){
        tcr_df <- filter(tcr, prior==x)
        tcr_df <- density(tcr_df$tcr)
        out <- 
          data_frame(
            tcr=tcr_df$x, 
            height=tcr_df$y,
            prior = x
          )
      }) %>%
      bind_rows()
    
    tcr_density %>%
      mutate(prior = factor(match_priors(prior), levels=prior_names)) %>%
      ggplot(aes(x=tcr, y =fct_reorder(prior, tcr), height=height, group=prior, col=prior, fill=prior)) +
      ## Dummy data (need to plot first otherwise annotate geom doesn't work)
      geom_density_ridges(stat = "identity", scale = 1.75, alpha = 0, col=NA) +
      ## IPCC "likely" region (1.0–2.5 °C)
      annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf, alpha = .2) +
      ## Priors
      geom_density_ridges(
        stat = "identity", scale = 1.75, 
        data = prior_dens,
        lty = 2, fill = NA
      ) +
      ## Posteriors
      geom_density_ridges(stat = "identity", scale = 1.75, alpha = 0.5, lwd = 0.5) +
      ## Stylistic elements
      labs(x = expression("TCR"~"("*degree*C*")"), y = "Density") +
      xlim(-1, 3) +
      scale_colour_manual(values = prior_cols) +
      scale_fill_manual(values = prior_cols) +
      theme(
        axis.text.y = element_text(vjust = 0),
        axis.title.y = element_blank(),
        legend.position = "none"
      ) 
  }

## Priors only version of the above
tcr_plot_priors <-
  function(tcr) {
    
    ## Can't use stat_function() that maps to facets, ridges or other aesthetic elements.
    ## So have to create the data manually instead.
    ## See: https://github.com/tidyverse/ggplot2/issues/2357
    p_df <- 
      priors_df %>%
      mutate(prior = paste0(prior_type, convic_type)) %>%
      mutate(prior = factor(match_priors(prior), levels=prior_names)) %>% 
      filter(prior_type!="ni")
    prior_dens <- 
      lapply(prior_names[1:4], function(x){
        df <- p_df %>% filter(prior==x)
        tcr_grid <- seq(from=qnorm(0.0001,df$mu,df$sigma), to=qnorm(0.9999,df$mu,df$sigma), length=100)
        data_frame(
          tcr = tcr_grid,
          height = dnorm(tcr_grid, mean=df$mu, sd=df$sigma),
          prior = x
          )
      }) %>%
      bind_rows()
    
    tcr_density <- 
      lapply(unique(tcr$prior), function(x){
        tcr_df <- filter(tcr, prior==x)
        tcr_df <- density(tcr_df$tcr)
        out <- 
          data_frame(
            tcr=tcr_df$x, 
            height=tcr_df$y,
            prior = x
            )
      }) %>%
      bind_rows()
    
    tcr_density %>%
      mutate(prior = factor(match_priors(prior), levels=prior_names)) %>%
      ggplot(aes(x=tcr, y =fct_reorder(prior, tcr), height=height, group=prior, col=prior, fill=prior)) +
      ## Dummy data (need to plot first otherwise annotate geom doesn't work)
      geom_density_ridges(stat = "identity", scale = 1.75, alpha = 0, col=NA) +
      ## IPCC "likely" region (1.0–2.5 °C)
      annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = Inf, alpha = .2) +
      ## Priors
      geom_density_ridges(
        stat = "identity", scale = 1.75, 
        data = prior_dens,
        lty = 2, fill = NA
        ) +
      ## Posteriors
      # geom_density_ridges(stat = "identity", scale = 1.75, alpha = 0.5, lwd = 0.5) +
      ## Stylistic elements
      labs(x = expression("TCR"~"("*degree*C*")"), y = "Density") +
      xlim(-1, 3) +
      scale_colour_manual(values = prior_cols) +
      scale_fill_manual(values = prior_cols) +
      theme(
        axis.text.y = element_text(vjust = 0),
        axis.title.y = element_blank(),
        legend.position = "none"
        ) 
  }

#################################################
#################################################
## Temperature's in 2100 pointrange plot function
all_2100_plot_func <-
  function(all_2100) {
    all_2100 %>%
      group_by(rcp, prior) %>%
      summarise(
        temp_mean = mean(temp, na.rm=T),
        temp_q025 = quantile(temp, p=0.025, na.rm=T),
        temp_q975 = quantile(temp, p=0.975, na.rm=T)
      ) %>%
      ungroup %>%
      mutate(prior = factor(match_priors(prior))) %>%
      mutate(rcp = match_rcps(rcp)) %>%
      ggplot(aes(x=fct_reorder(prior, temp_mean), y=temp_mean, ymin=temp_q025, ymax=temp_q975, col=prior)) +
      geom_pointrange() + 
      coord_flip() +
      scale_colour_manual(values = prior_cols) +
      facet_wrap(~ rcp) +
      labs(x = "Prior", y = "Temperature anomaly by 2100 (°C)") +
      theme(
        legend.position = "none",
        axis.title.y = element_blank()
      ) 
  }


##############################
##############################
## Recursive TCR plot function
recursive_plot_func <-
  function(tcr_rec) {
    tcr_rec %>%
      mutate(prior = factor(match_priors(prior), levels=prior_names[c(5,1:4)])) %>% ## Plot NI first
      ggplot(aes(x = year_to, y = mean, col = prior, fill = prior), lwd=0.5) +
      geom_line(aes(y=q025, lty=prior)) +
      geom_line(aes(y=q975, lty=prior)) +
      geom_ribbon(
        data= tcr_rec %>%
          filter(prior!="Noninformative") %>%
          mutate(prior = factor(match_priors(prior), levels=prior_names)),
        aes(ymin = q025, ymax = q975), lty = 0, alpha = 0.5
        ) +
      geom_line(lwd = .75) +
      labs(y = "Temperature anomaly (°C)\n") + 
      ## scale_y_continuous(limits = c(-1, 3)) +
      scale_colour_manual(values = rep(prior_cols, 2)) +
      scale_fill_manual(values = rep(c(prior_cols[1:4], "Noninformative"=NA), 2)) +
      scale_linetype_manual(
        values = c(0, 0, 0, 0, 2),
        limits = prior_names) +
      facet_wrap(~ priorlab, ncol = 2) +
      theme(
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none"
        )
    }


#########################
#########################
## Evidence plot function 

## Grid
evid_plot_func <-
  function(evid) {
    evid %>% 
      mutate(yrs = ifelse(is.na(yrs), 2100-1866, yrs)) %>%
      ggplot(aes(x = mu, y = sigma)) +
      geom_raster(aes(fill = yrs + 1866 -1)) +
      scale_fill_gradient2(
        name = "Year beliefs\nconverge",
        midpoint = 2015,
        labels = c("1950"="1950", "2000"="2000", "2050"="2050", "2100"="> 2100"),
        low="#377EB8", high="#E41A1C",
        limits = c(min(evid$yrs)+1866-1, 2100)
        ) +
      guides(fill = guide_colorbar(reverse = TRUE)) +
      # labs(x = expression(mu), y = expression(sigma)) +
      labs(x = "Prior mean (μ)", y = "Prior standard deviation (σ)") +
      facet_wrap(~thresh_lab)
  }

## Using lines instead of a grid
evid_plot_lines_func <-
  function(evid) {
    evid %>%
      filter(mu %in% round(seq(0, 1, by = .2), 1)) %>%
      ggplot(aes(x = sigma, y = yrs + 1866 - 1, group = factor(mu), col = factor(mu))) +
      geom_line() + 
      geom_line(
        data = evid %>% 
          group_by(mu) %>% 
          filter(is.na(yrs) | is.na(lead(yrs))) %>% 
          filter(mu < .5) %>%
          filter(mu %in% seq(0, 1, by = .2)),
        aes(x = sigma_dash - .001, y = yrs_dash + 1866), lty = 5
        ) +
      geom_hline(yintercept = 2015, col = "red", lty = 2) +
      geom_hline(yintercept = 2100, col = "red", lty = 1) +
      geom_text(
        data = evid %>% 
          filter(mu %in% round(seq(0, 1, by = .2), 1)) %>%
          filter(sigma == min(sigma)), 
        aes(label = sprintf('mu == "%1.1f"', mu)), 
        hjust = 0, nudge_x = .001, 
        parse = T, family = font_type, size = 3.5
        ) +
      geom_text(
        data = evid %>% filter(is.na(yrs) & mu == 0),
        aes(x = sigma, y = yrs_dash + 1866 - 1, label = sprintf('mu == "%1.1f"', mu)), 
        vjust = 0, nudge_y = 2, 
        hjust = 0, nudge_x = .001,  
        parse = T, family = font_type, size = 3.5
        ) +
      labs(
        x = expression(paste("Prior convinction (", sigma, ")")), 
        y = "Year beliefs converge"
        ) +
      scale_colour_manual(values = brewer.pal(9, "Blues")[4:9]) +
      scale_x_reverse(expand = c(0.2, 0)) +
      scale_y_continuous(expand = c(0.1, 0)) +
      facet_wrap(~thresh_lab) + theme(legend.position = "none")
  }


############################
############################
## SCC density plot function

scc_plot <-
  function(scc) {
    
    scc %>%
      gather(prior, scc) %>%
      mutate(prior = factor(match_priors(prior), levels=prior_names)) %>%
      ggplot(aes(x = scc, col = prior)) +
      geom_line(stat = "density") +
      xlim(-10, 100) + ## NB: x-axis is truncated to aid visual inspection!
      labs(
        # x = expression("Social cost of CO"[2]*" (US$2005)"), 
        x = "US$ (2005)", 
        y = "Density") + 
      scale_colour_manual(values = prior_cols) +
      guides(col = guide_legend(nrow = 2)) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )
  }
