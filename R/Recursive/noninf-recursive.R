##------------------------------------------------------------------------------
## GET PARAMETER SAMPLES 

theta_sample <- 
  blinreg(clim_df$had_sim, 
          cbind(1, clim_df$trf, clim_df$volc_sim, clim_df$soi_sim, clim_df$amo_sim), 
          chain_length
          )

## Get TCR
tcr <- 
  (as.matrix(theta_sample$beta)[,2] * rf2x) %>%
  as_data_frame() %>%
  magrittr::set_colnames("tcr") %>%
  mutate(
    prior = paste0(prior_type, convic_type),
    year_to = ifelse(recurse_type == "historic", yr_min, max(clim_df$year))
    )

return(tcr = tcr)