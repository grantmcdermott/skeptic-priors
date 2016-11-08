##------------------------------------------------------------------------------
## GET PARAMETER SAMPLES 

theta_sample <- 
  blinreg(clim_df$had_sim, 
          cbind(1, clim_df$trf, clim_df$volc_sim, clim_df$soi_sim, clim_df$amo_sim), 
          chain_length
          )

## Get TCR
tcr <-
  as.matrix(theta_sample$beta)[,2] %>%
  data_frame() %>% 
  magrittr::set_colnames("beta") %>%
  mutate(prior = paste0(prior_type, convic_type)) %>%
  mutate(year_to = ifelse(recurse_type == "historic", yr_min, max(clim_df$year)))

tcr$tcr <- tcr$beta * rf2x

tcr <- tcr %>% select_("-beta")

return(tcr = tcr)