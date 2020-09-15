data {
  int<lower=1> N;            // number of obs (historic only)
  real beta_mu;              // prior on climate resistance (via TCR)        
  real<lower=0> beta_sigma;  // prior on climate resistance (via TCR)
  real<lower=0> gamma_sigma; // weakly informative prior on volc
  real<lower=0> delta_sigma; // weakly informative prior on soi
  real<lower=0> eta_sigma;   // weakly informative prior on amo
  vector[N] gmst;
  vector[N] gmst_omega;      // measurement error (std deviation) on gmst series 
  vector[N] trf;    
  vector[N] volc;      
  vector[N] amo;      
  vector[N] soi;      
}
parameters {
  real<lower=0> sigma;
  real alpha;
  real beta;
  real gamma;
  real delta;
  real eta;
  real<lower=-1,upper=1> phi; // AR(1) param
}
transformed parameters {
  // Only estimating historic data (up to N periods)
  vector[N] mu;
  vector[N] sigmasq_tot;
  vector[N] epsilon;
  vector[N] sigma_cor;
  // First period
  mu[1] = alpha + beta*trf[1] + gamma*volc[1] + delta*soi[1] + eta*amo[1];
  epsilon[1] = gmst[1] - mu[1]; // observed minus predicted
  sigmasq_tot[1] = sigma*sigma + gmst_omega[1]*gmst_omega[1]; // Add ME component
  sigma_cor[1] = sqrt(sigmasq_tot[1] * (1-phi*phi));
  // Subsequent periods
  for(t in 2:N) {
    mu[t] = alpha + beta*trf[t] + gamma*volc[t] + delta*soi[t] + eta*amo[t] + phi*epsilon[t-1];
    sigmasq_tot[t] = sigma*sigma + gmst_omega[t]*gmst_omega[t]; // Add ME component
    epsilon[t] = gmst[t] - mu[t]; // observed minus predicted
    // Full AR(1) error
    sigma_cor[t] = sqrt(sigmasq_tot[t] * (1-phi*phi)); // Var = sigma2 * (1-rho^2)
  }
}
model {
  phi ~ normal(0,1);
  gmst ~ normal(mu, sigma_cor);
  sigma ~ cauchy(0,5);
  alpha ~ normal(0, 1);
  beta ~ normal(beta_mu, beta_sigma);
  gamma ~ normal(0, 1.92);
  delta ~ normal(0, 0.54);
  eta ~ normal(0, 3.37);
}