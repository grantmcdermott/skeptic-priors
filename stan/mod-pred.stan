data {
  int<lower=1> N1; // number of obs (historic only)
  int<lower=1> N2; // number of obs (incl. future)
  real beta_mu; // priors on climate resistance (via TCR)        
  real<lower=0> beta_sigma; // priors on climate resistance (via TCR)
  vector[N1] gmst;
  vector[N2] trf;    
  vector[N2] volc;      
  vector[N2] amo;      
  vector[N2] soi;      
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
  // Only estimating historic data (up to N1 periods)
  vector[N1] mu;
  vector[N1] epsilon;
  real<lower=0> sigma_cor;
  // First period
  mu[1] = alpha + beta*trf[1] + gamma*volc[1] + delta*soi[1] + eta*amo[1];
  epsilon[1] = gmst[1] - mu[1]; // observed minus predicted
  // Subsequent periods
  for(t in 2:N1) {
    mu[t] = alpha + beta*trf[t] + gamma*volc[t] + delta*soi[t] + eta*amo[t] + phi*epsilon[t-1];
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
  gamma ~ normal(0, 1.92);
  delta ~ normal(0, 0.54);
  eta ~ normal(0, 3.37);
}
generated quantities {
  // Generating preditions for full sample, including future (N2 periods)
  vector[N2] y_pred;
  vector[N2] mu_pred;
  vector[N2] epsilon_pred;
  // First period
  mu_pred[1] = alpha + beta*trf[1] + gamma*volc[1] + delta*soi[1] + eta*amo[1];
  y_pred[1] = normal_rng(mu_pred[1], sigma_cor);
  epsilon_pred[1] = y_pred[1] - mu_pred[1];
  // Subsequent periods (including future)
  for(t in 2:N2) {
    mu_pred[t] = alpha + beta*trf[t] + gamma*volc[t] + delta*soi[t] + eta*amo[t] + phi*epsilon_pred[t-1];
    y_pred[t] = normal_rng(mu_pred[t], sigma_cor);
    epsilon_pred[t] = y_pred[t] - mu_pred[t];
  }
}