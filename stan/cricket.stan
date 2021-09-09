//
// This Stan program fits a simple regression model.
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of data points
  vector[N] y;    // vector of responses
  vector[N] x;    // vector of predictors
}

// The parameters accepted by the model. 
parameters {
  real b0;              // intercept
  real b1;              // slope
  real<lower=0> sigma;  // SD of residual errors
}

// The model to be estimated. 
model {
  // priors
  b0 ~ normal(0, 10);
  b1 ~ normal(0, 10);
  sigma ~ inv_gamma(0.01, 0.01);
  
  // data model
  y ~ normal(b0 + b1 * x, sigma);
}

