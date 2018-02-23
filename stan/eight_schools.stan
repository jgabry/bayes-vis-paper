/* 8-schools centered parameterization */
data {
  int<lower=0> J;            // number of schools
  real y[J];                 // estimated treatment effect (school j)
  vector<lower=0>[J] sigma;  // std err of effect estimate (school j)
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] theta;
}
model {
  mu ~ normal(0, 100);
  tau ~ cauchy(0, 10);
  theta ~ normal(mu, tau); 
  y ~ normal(theta, sigma);
} 
