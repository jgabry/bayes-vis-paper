/* Non-centered parameterization (NCP) */

functions {
  /* Compute N-vector of linear predictors
   * 
   * @param intercept Scalar global intercept.
   * @param slope Scalar global slope.
   * @param intercept_shifts Vector of group-specific intercept shifts.
   * @param intercept_shifts Vector of group-specific slope shifts.
   * @param group_ids Integer array of group indexes.
   * @param x Vector of predictor values.
   * 
   * @return N-vector of linear predictors (N = num_elements(x)).
   * 
   * @throw Reject if intercept_shifts and slope_shifts don't have length
   *        equal to max(group_ids);
   * @throw Reject if x doesn't have same number of elements as group_ids. 
   */
  vector yhat(real intercept, 
              real slope, 
              vector intercept_shifts, 
              vector slope_shifts, 
              int[] group_ids, 
              vector x) {
    if (num_elements(x) != num_elements(group_ids))
      reject("'x' and 'group_ids' must have same length.");
    if (num_elements(intercept_shifts) != max(group_ids))
      reject("'intercept_shifts' should have max(group_ids) elements.");
    if (num_elements(slope_shifts) != max(group_ids))
      reject("'slope_shifts' should have max(group_ids) elements.");
      
    return intercept + intercept_shifts[group_ids] 
            + (slope + slope_shifts[group_ids]) .* x;
  }
}
data {
  int<lower=1> N;                   // number of observations
  int<lower=1> R;                   // number of super regions
  int<lower=1,upper=R> region[N];   // region IDs
  vector[N] log_sat;                // log of satellite measurements
  vector[N] log_pm;                 // log of ground PM_2.5 measurements
}
parameters {
  real beta0;                       // global intercept         
  real beta1;                       // global slope 
  vector[R] beta0_region_raw;       // 'raw' region intercept offsets for NCP
  vector[R] beta1_region_raw;       // 'raw' region slope offsets for NCP
  real<lower=0> tau0;               // sd of beta0_region
  real<lower=0> tau1;               // sd of beta1_region
  real<lower=0> sigma;              // error sd for Gaussian likelihood
}
transformed parameters {
  // Scale and shift raw parameters for NCP 
  // (but no shift in this case since mean 0)
  vector[R] beta0_region = tau0 * beta0_region_raw; 
  vector[R] beta1_region = tau1 * beta1_region_raw;
}
model {
  // Log-likelihood
  vector[N] log_pm_hat = yhat(beta0, beta1, 
                              beta0_region, beta1_region, 
                              region, log_sat);
  target += normal_lpdf(log_pm | log_pm_hat, sigma);
    
  // Log-priors
  target += normal_lpdf(sigma | 0, 1)
          + normal_lpdf(tau0  | 0, 1)
          + normal_lpdf(tau1  | 0, 1)
          + normal_lpdf(beta0 | 0, 1)
          + normal_lpdf(beta1 | 1, 1)
          + normal_lpdf(beta0_region_raw | 0, 1)
          + normal_lpdf(beta1_region_raw | 0, 1);
}
generated quantities {
  vector[N] log_lik;     // pointwise log-likelihood for LOO
  vector[N] log_pm_rep;  // replications from posterior predictive dist

  { // local block to not save log_pm_hat
    vector[N] log_pm_hat = yhat(beta0, beta1, 
                                beta0_region, beta1_region, 
                                region, log_sat);
    for (n in 1:N) {
      log_lik[n] = normal_lpdf(log_pm[n] | log_pm_hat[n], sigma);
      log_pm_rep[n] = normal_rng(log_pm_hat[n], sigma);
    }
  }
}
