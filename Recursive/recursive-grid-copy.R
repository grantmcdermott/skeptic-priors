rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")
library(tibble)
library(viridis)

####
## NOTE: SKIP TO LINE 445 ONCE THE BELOW CODE HAS ALREADY BEEN RUN ONCE. TAKES SOME TIME TO
## RUN AND NO NEED TO REPEAT ONCE THE TARGET DATA HAS ALREADY BEEN PRODUCED


## Optional for replication
set.seed(123) 

## Load climate data
climate <- read_csv("./Data/climate.csv")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
chain_length <- 3000
n_chains <- detectCores() - 1 

## Only need to compare forcings and (relevant) predicted temps, so which RCP doesn't 
## really matter.
rcp_type <- "rcp26" 

## Load climate data and deviation DF for simulating "true" future values
climate <- 
  read_csv("./Data/climate.csv") %>%
  filter(rcp == rcp_type) #%>%
  # filter(year <= yr_max) 

y_dev <- read_csv("./Recursive/Data/y-dev.csv")

climate$noise <- sample(y_dev$dev, nrow(climate), replace = T)

## Also add noise to future "mean" covariate values; otherwise collinearity problems for some
## regressions.
climate$volc_noise <- rnorm(nrow(climate))
climate$soi_noise <- rnorm(nrow(climate))
climate$amo_noise <- rnorm(nrow(climate))

climate <- 
  climate %>% 
  mutate(volc_sim = ifelse(year <= 2005, volc_mean, volc_mean + volc_noise),
         soi_sim = ifelse(year <= 2005, soi_mean, soi_mean + soi_noise),
         amo_sim = ifelse(year <= 2005, amo_mean, amo_mean + amo_noise)
         ) %>%
  mutate(had_sim = 
           ifelse(year <= 2005, 
                  had,
                  -0.110 + 0.415*trf + 0.047*volc_sim + -0.028*soi_sim + 0.481*amo_sim + noise
                  )
         )


## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length * n_chains, mean = 3.71, sd = 0.1855) 


## Choose posterior TCR threshold. 
## Mean TCR values must be at least this big to qualify as fully converged with mainstream. 
thresh <- 1.5


#########################################################
### START OF HISTORIC RECURSIVE UPDATING #####
##############################################

mu <- seq(1, 0, length = 11)
sigma <- seq(.25, .0674, length = 11)

yr_max <- 2005
yr_min <- min(climate$year)

yrs_df <- list()
a <- 0

ptm <- proc.time()
pb <- txtProgressBar(min = 0, max = length(sigma), style = 3)
## Go from the mildest sceptic to the most hardcore denier
for(i in 1:length(sigma)){
  
  if(i == 1){
    yrs <- 95#100
  }else{
    yrs <- yrs_j
  }
  
  for(j in 1:length(mu)){
    
    # if((yrs == 140) & (round(tcr_mean, digits = 1) < thresh)){
    #   yrs = 140
    #   tcr_mean = NA
    # }else{}

    yrs <- yrs - 1#5
    tcr_mean <- 0 ## Place holder value
    
    while((round(tcr_mean, digits = 1) < thresh) & (yrs < nrow(filter(climate, year <= yr_max)))) {
    # while((round(tcr_mean, digits = 1) < thresh) & (yrs < nrow(climate))) {
      
      yrs <- yrs + 1#5
      
      clim_df <-
        climate %>%
        filter(year <= yr_max) %>%
        filter(year >= max(year) - yrs)
      
      ##------------------------------------------------------------------------------
      ## THE BUGS/JAGS MODEL.
      ## Note: Removing y_pred b/c predictions into the future not needed for recursive regs
      N <- nrow(clim_df)
      
      mu_beta <- mu[j]/3.71
      sigma_beta <- sigma[i]/3.71
      
      mod_string <- paste(
        "model{
        
        for(t in 1:N) {
        mu[t] <- alpha + beta*trf[t] + gamma*volc_sim[t] + delta*soi_sim[t] + eta*amo_sim[t]
        had_sim[t]  ~ dnorm(mu[t], tau)
        # y_pred[t] ~ dnorm(mu[t], tau) ## For predictions into the future
        }
        ", 
        paste0(
        "
        mu_beta <- ", mu_beta),
        paste0(
        "
        sigma_beta <- ", sigma_beta),
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
      
      bugs_file <- "./Recursive/BUGSfiles/recursive-bugs.txt"
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
      tcr_mean <-mean(tcr)
      
      rm(bugs_file, cl, data_list, mod_samples, mod_string, par_inits, parameters, tcr)
      
      #### end of while loop
    }
    
    if(j == 1){yrs_j <- yrs}
    
    a <- a + 1
    yrs_df[[a]] <- data_frame(mu = mu[j], 
                             sigma = sigma[i],
                             tcr_mean = tcr_mean,
                             yrs = yrs
                             # yrs = ifelse(round(tcr_mean, digits = 1) < thresh, NA, yrs)
                             )
    
  }
  setTxtProgressBar(pb, i)
}
close(pb)
ptm <- proc.time() - ptm
# user  system elapsed 
# 5.66    2.21  466.14
# user  system elapsed 
# 14.29    5.17 1119.70 

yrs_df <-
  bind_rows(yrs_df)

write_csv(yrs_df, "./Recursive/Data/tcr-rec-grid-historic.csv")


###############
## Read in data

yrs_df <-
  read_csv("./Recursive/Data/tcr-rec-grid-historic.csv") #%>% filter(mu %in% c(1.0, 0.0))

yrs_df %>%
  mutate(yrs = ifelse(round(tcr_mean, digits = 1) < thresh, NA, yrs)) %>%
  ggplot(aes(x = mu, y = sigma, z = yrs)) +
  geom_raster(aes(fill = yrs)) + 
  scale_fill_viridis(name = "Years",
                     direction = -1,
                     guide = guide_colourbar(reverse = TRUE)) +
  # scale_fill_gradient(name = "Years",
  #                     low = "lightblue", high = "blue",
  #                     guide = guide_colourbar(reverse = TRUE)) +
  labs(x = expression(mu), y = expression(sigma)) +
  theme(text = element_text(family = font_type))



#########################################################
### START OF FUTURE RECURSIVE UPDATING #####
############################################
fyrs <-
  yrs_df %>% 
  filter(round(tcr_mean, digits = 1) < thresh)

sigma <- 
  rev(sort((fyrs %>%
              distinct(sigma))$sigma))
# 0.17696 0.15870 0.14044 0.12218 0.10392 0.08566 0.06740

## The catch is these sigma values are only relevant for certain mu values.
## Need to exclude these at the relevant points in the nested loop.
fyrs <-
  fyrs %>%
  arrange(desc(mu)) 

incl_list <-
  lapply(mu, function(m){
    (fyrs %>%
       arrange(desc(mu)) %>%
       filter(mu == m))$sigma
    })


yr_max <- 2100#2050 

fyrs_df <- list()
a <- 0

ptm <- proc.time()
pb <- txtProgressBar(min = 0, max = length(sigma), style = 3)
## Go from the mildest sceptic to the most hardcore denier
for(i in 1:length(sigma)){
  
  if(i == 1){
    yrs <- 140
  }else{
    yrs <- yrs_j
  }
  
  for(j in 1:length(mu)){
    
    if(sigma[i] %nin% incl_list[[j]]){
      # yrs <- NA
      tcr_mean <- NA
    }else{
    
    yrs <- yrs - 1
    tcr_mean <- 0 ## Place holder value
    
    while((round(tcr_mean, digits = 1) < thresh) & (yrs < nrow(filter(climate, year <= yr_max)))) {
      # while((round(tcr_mean, digits = 1) < thresh) & (yrs < nrow(climate))) {
      
      yrs <- yrs + 1#5
      
      clim_df <-
        climate %>%
        filter(year <= min(year) + yrs)
      
      ##------------------------------------------------------------------------------
      ## THE BUGS/JAGS MODEL.
      ## Note: Removing y_pred b/c predictions into the future not needed for recursive regs
      N <- nrow(clim_df)
      
      mu_beta <- mu[j]/3.71
      sigma_beta <- sigma[i]/3.71
      
      mod_string <- paste(
        "model{
        
        for(t in 1:N) {
        mu[t] <- alpha + beta*trf[t] + gamma*volc_sim[t] + delta*soi_sim[t] + eta*amo_sim[t]
        had_sim[t]  ~ dnorm(mu[t], tau)
        # y_pred[t] ~ dnorm(mu[t], tau) ## For predictions into the future
        }
        ", 
        paste0(
          "
          mu_beta <- ", mu_beta),
        paste0(
          "
          sigma_beta <- ", sigma_beta),
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
      
      bugs_file <- "./Recursive/BUGSfiles/recursive-bugs-future.txt"
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
      tcr_mean <-mean(tcr)
      
      rm(bugs_file, cl, data_list, mod_samples, mod_string, par_inits, parameters, tcr)
      
      #### end of while loop
    }
    
    if(j == 1){yrs_j <- yrs}
    
    } ## End of else clause on inc_list 
    
    a <- a + 1
    fyrs_df[[a]] <- data_frame(mu = mu[j], 
                              sigma = sigma[i],
                              tcr_mean = tcr_mean,
                              yrs = yrs
                              # yrs = ifelse(round(tcr_mean, digits = 1) < thresh, NA, yrs)
    )
    
    }
  setTxtProgressBar(pb, i)
  }
close(pb)
ptm <- proc.time() - ptm
# user  system elapsed 
# 6.42    2.38  548.53
### That's for 2050. Going to 2100 yields an additional:
# user  system elapsed 
# 2.27    0.69  185.39 

fyrs_df <-
  bind_rows(fyrs_df)


write_csv(fyrs_df, "./Recursive/Data/tcr-rec-grid-future.csv")

fyrs_df <- read_csv("./Recursive/Data/tcr-rec-grid-future.csv")


fyrs_df %>%
  mutate(yrs = ifelse(round(tcr_mean, digits = 1) < thresh, NA, yrs)) %>%
  ggplot(aes(x = mu, y = sigma, z = yrs)) +
  geom_raster(aes(fill = yrs)) + 
  scale_fill_viridis(name = "Years",
                     direction = -1,
                     guide = guide_colourbar(reverse = TRUE)) +
  # scale_fill_gradient(name = "Years",
  #                     low = "lightblue", high = "blue",
  #                     guide = guide_colourbar(reverse = TRUE)) +
  labs(x = expression(mu), y = expression(sigma)) +
  theme(text = element_text(family = font_type))


yrs_comb <- 
  bind_rows(yrs_df %>% filter(round(tcr_mean, digits = 1) >= thresh),
          fyrs_df %>% filter(!is.na(tcr_mean))
          ) 

write_csv(yrs_comb, "./Recursive/Data/tcr-rec-grid.csv")


#####################################################
#####################################################
## Read in final data

yrs_comb <- read_csv("./Recursive/Data/tcr-rec-grid.csv")


p <- 
  yrs_comb %>%
  mutate(yrs = ifelse(round(tcr_mean, digits = 1) < thresh, NA, yrs)) %>%
  ggplot(aes(x = mu, y = sigma)) +
  geom_raster(aes(fill = yrs)) + 
  scale_fill_viridis(name = "Years",
                     direction = -1,
                     guide = guide_colourbar(reverse = TRUE)) +
  labs(x = expression(mu), y = expression(sigma)) +
  theme(text = element_text(family = font_type),
        axis.title.y = element_text(angle=0))


yrs_df <- read_csv("./Recursive/Data/tcr-rec-grid-historic.csv")
bound <-
  yrs_df %>%
  filter(yrs == 140) %>%
  distinct(sigma) %>%
  # group_by(mu) %>% filter(sigma == max(sigma)) %>% ungroup() %>%
  mutate(mu_1 = mu - .05, sigma_1 = sigma + 0.01826/2,
         mu_2 = mu + .05, sigma_2 = sigma + 0.01826/2,
         mu_3 = mu + .05, sigma_3 = sigma - 0.01826/2) %>%
  select(mu_1:sigma_3) %>%
  mutate(grp = row_number()) %>%
  gather(key, value, -grp) %>%
  separate(key, into = c("key", "suff")) %>%
  mutate(grp = paste0(grp, "_", suff)) %>%
  select(-suff) %>%
  spread(key, value) %>% 
  select(-grp) %>% 
  mutate(mu = round(mu, 2)) %>%
  distinct(mu, sigma) %>%
  filter(sigma > 0.07653) %>%
  arrange(mu, desc(sigma))

p +
  geom_line(data = bound, col = "red", lwd = 1.25) +
  ggsave(file = "./Recursive/TablesFigures/recursive_grid.pdf",
         width = 6, height = 4,
         device = cairo_pdf)

