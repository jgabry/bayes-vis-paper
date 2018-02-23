/* Simple linear regression */
data {
  int<lower=1> N;       // number of observations
  vector[N] log_sat;    // log of satellite measurements
  vector[N] log_pm;     // log of ground PM_2.5 measurements
}
parameters {
  real beta0;           // global intercept
  real beta1;           // global slope
  real<lower=0> sigma;  // error sd for Gaussian likelihood
}
model {
  // Log-likelihood
  target += normal_lpdf(log_pm | beta0 + beta1 * log_sat, sigma);

  // Log-priors
  target += normal_lpdf(sigma | 0, 1)
          + normal_lpdf(beta0 | 0, 1)
          + normal_lpdf(beta1 | 1, 1);
}
generated quantities {
  vector[N] log_lik;    // pointwise log-likelihood for LOO
  vector[N] log_pm_rep; // replications from posterior predictive dist

  for (n in 1:N) {
    real log_pm_hat_n = beta0 + beta1 * log_sat[n];
    log_lik[n] = normal_lpdf(log_pm[n] | log_pm_hat_n, sigma);
    log_pm_rep[n] = normal_rng(log_pm_hat_n, sigma);
  }
}
