rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for 
## plotting and tables
source("sceptic_funcs.R")

## Optional for replication
set.seed(123) 


## Load climate data
climate <- read_csv("Data/climate.csv")

## Decide on length of MCMC chains (including no. of chains in parallel JAGS model)
chain_length <- 9000
n_chains <- max(sapply(1:detectCores(), function(x) gcd(x, chain_length)))

## Only need to compare forcings and (relevant) predicted temps, so which RCP  
## doesn't really matter.
rcp_type <- "rcp26" 

## Load climate data and deviation DF for simulating "true" future values
climate <- 
  read_csv("Data/climate.csv") %>%
  filter(rcp == rcp_type) 

y_dev <- read_csv("Evidence/Data/y-dev.csv")

climate$noise <- sample(y_dev$dev, nrow(climate), replace = T)

## Also add noise to future "mean" covariate values; otherwise collinearity 
## problems for future regressions b/c values stay constant.
climate$volc_noise <- rnorm(nrow(climate), sd = .1)
climate$soi_noise <- rnorm(nrow(climate), sd = .1)
climate$amo_noise <- rnorm(nrow(climate), sd = .1)

## Simulate climate data based on noninformative parameters
climate <- 
  climate %>% 
  mutate(volc_sim = ifelse(year <= 2005, volc_mean, volc_mean + volc_noise),
         soi_sim = ifelse(year <= 2005, soi_mean, soi_mean + soi_noise),
         amo_sim = ifelse(year <= 2005, amo_mean, amo_mean + amo_noise)
         ) %>%
  mutate(had_sim = -0.102 + 0.418*trf + 0.051*volc_sim + -0.028*soi_sim + 0.473*amo_sim + noise)


## Set radiative forcing distribution used for calulating TCRs later in code.
## Centered around 3.71 °C +/- 10% (within 95% CI). 
## Length of disbn equals length of MCMC chain for consistency
rf2x <- rnorm(chain_length, mean = 3.71, sd = 0.1855) 

## EVIDENCE FUNCTION
## DESCRIPTION:
## The function loops along a sequence of (decreasing) sigma values for a given 
## mu. Once it reaches the end of the sequence for that mu, it moves on to the  
## next, lower mu value and repeats the procedure. It does this for two TCR
## threshold levels: 1.3 °C and 1.5°C. (Note that the function will need to be 
## adjusted  if different thresholds are chosen.) The mean posterior TCR values  
## must be at least as big as the relevant threshold to qualify as fully converged  
## with mainstream. 
evid_func <- 
  function(d){
    pblapply(1:nrow(d), function(x){
      
      m <- d[x, ]$mu
      s <- d[x, ]$sigma
      max_m <- d[x, ]$max_mu
      thresh <- d[x, ]$threshold
      
      yr_max <- 2100

      yrs <<- ifelse(thresh == 1.5 & round(m, 1) == 1.0 & round(s, 3) == 0.250, #x == 122, ## i.e. first row where thresh == 1.5
                     100, 
                     ifelse(thresh == 1.3 & x == 1, 
                            65, 
                            ifelse(max_m != 1, 
                                   yrs - 1, 
                                   yrs_j))) 
      
      tcr_mean <- 0 ## Place holder value
      
      while((round(tcr_mean, digits = 1) < thresh) & 
            (yrs < nrow(filter(climate, year <= yr_max)))) {
        
        ## NB: Use "<<-" to assign value to global workspace for next iteration
        yrs <<- yrs + 1 
        
        clim_df <-
          climate %>%
          filter(year <= min(year) + yrs)
        
        ##----------------------------------------------------------------------
        ## THE BUGS/JAGS MODEL.
        ## Note: Removing y_pred b/c predictions into the future not needed for
        ## recursive regs.
        N <- nrow(clim_df)
        
        mu_beta <- m/3.71
        sigma_beta <- s/3.71
        
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
        
        bugs_file <- paste("Evidence/BUGSFiles/evidence-bugs.txt")
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
        mod_iters <- chain_length/n_chains
        mod_samples <- parCodaSamples(cl, "jags_mod", variable.names = parameters, 
                                      n.iter = mod_iters, n.chain = n_chains) 
        stopCluster(cl)
        
        ##------------------------------------------------------------------------------
        ## Get TCR coefficient (i.e. beta) and then use to get TCR.
        
        tcr <- as.matrix(mod_samples[, "beta"], iters = F) * rf2x
        tcr_mean <- mean(tcr)
        
        rm(bugs_file, cl, data_list, mod_samples, mod_string, par_inits, parameters, tcr)
        
      } ### end of while loop
      
      ## "Restart" the iteration each time we get back to the "max mu" value for that
      ## sigma loop. (Max mu will always be at mu = 1 for historic data frame, but not 
      ## necessarily for future data frame.) Will begin at convergence year for the 
      ## previous sigma value.
      ## Again, use "<<-" to assign value to global workspace for next iteration
      if(max_m == 1){yrs_j <<- yrs - 1}
      
      df <- 
        data_frame(mu = m,
                   sigma = s,
                   tcr_mean = tcr_mean,
                   yrs = yrs,
                   thresh = thresh)
      
      return(df)
      
    } ## end of pblapply function
    ) %>% ## end of pbapply
      bind_rows()
  } ## end of entire evid_func function


###################################
## Sceptic prior range (data frame)

mu <- round(seq(1, 0, length = 11), 1)
sigma <- round( seq(.25, .065, length = 11), 4)
threshold <- c(1.3, 1.5)

df <- 
  expand.grid(mu = mu, sigma = sigma, threshold = threshold) %>%
  group_by(threshold, sigma) %>% 
  mutate(max_mu = ifelse(mu == max(mu), 1, 0)) %>%
  ungroup() 

####
## Run the function over the full range of sceptic priors
## NOTE: Takes about an hour to run on my laptop (16 GB RAM), so skip to line 224
## once it has already been run and the data exported to file.
evid <- evid_func(df)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% Elapsed time: 55m 09s
rm(yrs, yrs_j)

## Write data
write_csv(evid, "Evidence/Data/tcr-evidence.csv")

## Read data
evid <- read_csv("Evidence/Data/tcr-evidence.csv")

## Include labels for facetting and manually adjust some parts of the data for
## plotting (only applies to obs where no convergence by 2100)
evid <-
  evid %>%
  mutate(yrs_dash = ifelse(yrs == 235, yrs - 1, yrs)) %>% 
  mutate(yrs = ifelse(round(tcr_mean, 1) < thresh, NA, yrs)) %>%
  mutate(sigma_dash = ifelse(round(tcr_mean, 1) < thresh, (0.0674+(.2-mu)/30)^1.2+0.03, sigma)) %>%
  mutate(thresh_lab = paste0(thresh, " °C")) %>%
  mutate(thresh_lab = ifelse(thresh == 1.3, 
                             paste0("(a) ", thresh_lab), 
                             paste0("(b) ", thresh_lab)))

## Plot the data
## Years with red-white-blue colour scheme
evid_plot <- evid_plot_func(evid)
evid +
  ggsave(
    file = "Evidence/TablesFigures/PNGs/evidence-grid.png",
    width = 8, height = 4
    )
evid +
  ggsave(
    file = "Evidence/TablesFigures/evidence-grid.pdf",
    width = 8, height = 4,
    device = cairo_pdf
    )
rm(evid_plot)

## Lines instead of grid
evid_plot_lines <- evid_plot_lines_func(evid) 
evid_plot_lines +
  ggsave(
    file = "Evidence/TablesFigures/PNGs/evidence-lines.png",
    width = 8, height = 4
    )
evid_plot_lines +
  ggsave(
    file = "Evidence/TablesFigures/evidence-lines.pdf",
    width = 8, height = 4,
    device = cairo_pdf
    )
rm(evid_plot_lines)
