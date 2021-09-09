//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as Poisson distributed
// with mean 'lambda'.
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> N; // number of data points
  int<lower=0> y[N]; // data vector
}

// The parameters accepted by the model. 
parameters {
  real<lower=0> lambda;
}

// The model to be estimated. We model the output
// 'y' to be Poisson distributed with mean 'lambda'
// and a gamme prior on lambda.
model {
  // priors
  lambda ~ gamma(0.01, 0.01);
  
  // data model
  y ~ poisson(lambda);
  // alt: target += poisson_lpdf(y | lambda)

}

