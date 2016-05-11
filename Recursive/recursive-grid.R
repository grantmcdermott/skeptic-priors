rm(list = ls()) # Clear data

## Load all packages, as well as some helper functions that will be used for plotting and tables
source("sceptic_funcs.R")
library(tibble)
library(pbapply) ## Add progress bar to *apply functions
library(viridis)

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


######################
## RECURSIVE FUNCTION

recurs_func <- 
  function(d){
    pblapply(1:nrow(d), function(x){
      
      m <- d[x, ]$mu
      s <- d[x, ]$sigma
      
      if(d[1, ]$rec == "historic"){yr_max <- 2005}
      if(d[1, ]$rec == "future"){yr_max <- 2100}
      
      if(d[1, ]$rec == "historic"){
        yrs <<- ifelse(thresh == 1.5 & x == 1,
                       95,
                       ifelse(thresh == 1.3 & x == 1,
                              45, 
                              ifelse(m != 1, 
                                     yrs - 1, 
                                     yrs_j
                                     )
                              )
                       )
      }
      
      if(d[1, ]$rec == "future"){
        yrs <<- ifelse(x == 1,
                       140,
                       ifelse(m != 1, 
                              yrs - 1, 
                              yrs_j
                              )
                       )
        }
      
      tcr_mean <- 0 ## Place holder value
      
      while((round(tcr_mean, digits = 1) < thresh) & (yrs < nrow(filter(climate, year <= yr_max)))) {
        
        ## NB: Use "<<-" to assign value to global workspace for next iteration
        yrs <<- yrs + 1 
        
        ## For historic recursion, we go backwards one year at a time
        if(d[1, ]$rec == "historic"){
          clim_df <-
            climate %>%
            filter(year <= yr_max) %>%
            filter(year >= max(year) - yrs)
        }
        
        ## For future recursion, we go forwards one year at a time
        if(d[1, ]$rec == "future"){
          clim_df <-
            climate %>%
            filter(year <= min(year) + yrs)
        }
      
        ##------------------------------------------------------------------------------
        ## THE BUGS/JAGS MODEL.
        ## Note: Removing y_pred b/c predictions into the future not needed for recursive regs
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
        
        bugs_file <- paste0("./Recursive/BUGSfiles/recursive-bugs-tcr", 
                            gsub("\\.", "_", thresh),
                            ".txt")
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
        tcr_mean <- mean(tcr)
        
        rm(bugs_file, cl, data_list, mod_samples, mod_string, par_inits, parameters, tcr)
        
      } ### end of while loop
      
      ## "Restart" the iteration each time we get back to a mu of 1
      ## Again, use "<<-" to assign value to global workspace for next iteration
      if(m == 1){yrs_j <<- yrs - 1}
      
      df <- 
        data_frame(mu = m,
                   sigma = s,
                   tcr_mean = tcr_mean,
                   yrs = yrs,
                   thresh = thresh
                   )
      
      return(df)
      
    } ## end of pblapply function
  ) %>% ## end of pbapply
  bind_rows()
    
} ## end of entire recurse_func function


#########################################################


## Choose posterior TCR threshold. 
## Mean TCR values must be at least this big to qualify as fully converged with mainstream. 
thresh <- c(1.3, 1.5)[1]


####
## NOTE: SKIP TO LINE 290 ONCE THE BELOW CODE HAS ALREADY BEEN RUN ONCE. TAKES QUITE A WHILE
## TO EXECUTE AND NO NEED TO REPEAT ONCE THE TARGET DATA HAS ALREADY BEEN PRODUCED

################################
## Historic recursive estimates

mu <- seq(1, 0, length = 11)
sigma <- seq(.25, .0674, length = 11)

hist <- 
  expand.grid(mu = mu, sigma = sigma) %>%
  mutate(rec = "historic")

hrec <- recurs_func(hist)
## For TCR = 1.3
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% Elapsed time: 29m 44s
rm(yrs, yrs_j)

## Write data
write_csv(hrec, paste0("./Recursive/Data/tcr",
                         gsub("\\.", "_", thresh),
                         "-rec-grid-historic.csv"))


################################
## Future recursive estimates

futr <-
  yrs_df %>% 
  filter(round(tcr_mean, digits = 1) < thresh) %>%
  select(mu, sigma) %>%
  mutate(rec = "future")

frec <- recurs_func(futr)
## For TCR = 1.3
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% Elapsed time: 02m 41s
rm(yrs)

## Write data
write_csv(frec, paste0("./Recursive/Data/tcr",
                         gsub("\\.", "_", thresh),
                         "-rec-grid-future.csv"))


################################
## Combined recursive estimates

rec <-bind_rows(hrec, frec) 

## Write data
write_csv(rec, paste0("./Recursive/Data/tcr",
                       gsub("\\.", "_", thresh),
                       "-rec-grid-comb.csv"))



###########################################################################
## Now compare combined recursive estimates from both 1.3 and 1.5 mean TCR

## Read data
rec <-
  bind_rows(
    lapply(c("1_3", "1_5"), function(x){
      read_csv(paste0("./Recursive/Data/tcr", x, "-rec-grid-comb.csv"))
    })
  )
  

p <-
  rec %>%
  mutate(yrs = ifelse(round(tcr_mean, digits = 1) < thresh, NA, yrs)) %>%
  mutate(thresh_lab = paste0(thresh, " °C")) %>%
  mutate(thresh_lab = ifelse(thresh == 1.3, 
                             paste0("(a) ", thresh_lab), 
                             paste0("(b) ", thresh_lab))) %>%
  ggplot(aes(x = mu, y = sigma)) +
  geom_raster(aes(fill = yrs)) +
  scale_fill_viridis(name = "Years",
                     direction = -1,
                     guide = guide_colourbar(reverse = TRUE)) +
  labs(x = expression(mu), y = expression(sigma)) +
  facet_wrap(~thresh_lab#, labeller = label_parsed
             ) +
  theme(text = element_text(family = font_type),
        axis.title.y = element_text(angle=0),
        strip.background = element_rect(fill = "white"), ## Facet strip
        panel.margin = unit(2, "lines"))


## Place bounds around future and historic values
bounds <-
  bind_rows(
    lapply(c("1_3", "1_5"), function(x){
      read_csv(paste0("./Recursive/Data/tcr", x, "-rec-grid-historic.csv")) %>%
        # arrange(desc(sigma), mu) %>%
        filter(yrs == 140) %>%
        distinct(sigma) %>%
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
        filter(sigma > ifelse(gsub("_", ".", x) == 1.5, 0.07653, 0)) %>%
        arrange(mu, desc(sigma)) %>%
        mutate(thresh = gsub("_", ".", x))
    })
  ) %>%
  mutate(thresh_lab = paste0(thresh, " °C")) %>%
  mutate(thresh_lab = ifelse(thresh == 1.3, 
                             paste0("(a) ", thresh_lab), 
                             paste0("(b) ", thresh_lab))) 
## Manual adjustment to make plot look nicer
bounds <-
  bounds %>% group_by(thresh) %>% #filter(mu == min(mu)) %>%
  mutate(mu = ifelse(mu == min(mu) & mu != -.05, -.05, mu))

p +
  geom_line(data = bounds, col = "red", lwd = 1.25) +
  # geom_text(aes(label = yrs), family = font_type) +
  ggsave(file = "./Recursive/TablesFigures/recursive-grid.pdf",
         width = 8, height = 4,
         device = cairo_pdf)

