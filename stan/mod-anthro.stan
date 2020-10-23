data {
  int<lower=1> N;                 // number of obs (historic only)
  real beta_mu;                   // prior on climate resistance (via TCR)        
  real<lower=0> beta_sigma;       // prior on climate resistance (via TCR)
  real<lower=0> beta_other_sigma; // weakly informative prior on trf_other
  real<lower=0> gamma_sigma;      // weakly informative prior on volc
  real<lower=0> delta_sigma;      // weakly informative prior on soi
  real<lower=0> eta_sigma;        // weakly informative prior on amo
  vector[N] gmst;
  vector[N] trf_anthro;
  vector[N] trf_other;    
  vector[N] volc;     
  vector[N] soi;  
  vector[N] amo;      
}
parameters {
  real<lower=0> sigma;
  real alpha;
  real beta;
  real beta_other;
  real gamma;
  real delta;
  real eta;
  real<lower=-1,upper=1> phi; // AR(1) param
}
transformed parameters {
  // Only estimating historic data (up to N periods)
  vector[N] mu;
  vector[N] epsilon;
  real<lower=0> sigma_cor;
  // First period
  mu[1] = alpha + beta*trf_anthro[1] + beta_other*trf_other[1] + gamma*volc[1] + delta*soi[1] + eta*amo[1];
  epsilon[1] = gmst[1] - mu[1]; // observed minus predicted
  // Subsequent periods
  for(t in 2:N) {
    mu[t] = alpha + beta*trf_anthro[t] + beta_other*trf_other[t] + gamma*volc[t] + delta*soi[t] + eta*amo[t] + phi*epsilon[t-1];
    epsilon[t] = gmst[t] - mu[t]; // observed minus predicted
  }
  // Full AR(1) error
  sigma_cor = sqrt(sigma*sigma * (1-phi*phi)); // Var = sigma2 * (1-rho^2)
}
model {
  phi ~ normal(0,1);
  gmst ~ normal(mu, sigma_cor);
  sigma ~ cauchy(0,5);
  alpha ~ normal(0, 1);
  beta ~ normal(beta_mu, beta_sigma);
  beta_other ~ normal(0, beta_other_sigma);
  gamma ~ normal(0, gamma_sigma);
  delta ~ normal(0, delta_sigma);
  eta ~ normal(0, eta_sigma);
}