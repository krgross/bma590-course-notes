# Bayesian computation

This chapter of the computing companion will focus solely on the computing aspects of Bayesian computation in R.  See the course notes or relevant sections of Bolker for the underlying theory.

The landscape of computing tools available to fit Bayesian models is fluid.  Here, we will look at three tools currently available: `R2jags`, which is based on the JAGS (Just Another Gibbs Sampler) platform, `rstan`, which is based on the computer program Stan (itself based on Hamiltonian Monte Carlo, or HMC), and the recent `rstanarm`, which seeks to put much of the computational details in the background.  (The "arm" portion of the name `rstanarm` is an acronym for applied regression modeling.)  

Throughout, we will be working with two data sets: the horse kick data (again), and a data set that details how the rate at which a cricket chirps depends on the air temperature.  The horse kick data are useful in this context because a Gamma distribution is a conjugate prior for Poisson data.  Thus, if we use a Gamma prior, then we know the posterior exactly.  Therefore, we can compare the approximations provided by stochastic sampling schemes to the known posterior.  The cricket data set will be used as an example of a simple linear regression, even though the data hint that the actual relationship between temperature and the rate of chirping is nonlinear.

## Computations with conjugate priors

Suppose that we observe an iid random sample $X_1, \ldots, X_n$ from a Poisson distribution with unknown parameter $\lambda$.  (This is the setting for the horse-kick data.)  If we place a Gamma prior with shape parameter $a$ and rate parameter $r$ on $\lambda$, then the posterior distribution is also Gamma with shape parameter $a + \sum_n X_n$ and rate parameter $r + n$.  In other words,
\begin{align*}
\lambda & \sim \mbox{Gamma}(a, r) \\
X_1, \ldots, X_n & \sim \mbox{Pois}(\lambda) \\
\lambda | X_1, \ldots, X_n & \sim \mbox{Gamma}(a + \sum_n X_n, r + n) \\
\end{align*}

In the horse kick data, $\sum_n x_n = 196$ and $n = 280$.  Suppose we start with the vague Gamma prior $a=.01$, $r = .01$ on $\lambda$.  This prior has mean $a/r = 1$ and variance $a/r^2 = 100$.  The posterior distribution for $\lambda$ is then a Gamma with shape parameter $a = 196.01$ and rate parameter $280.01$.  We can plot it:


```r
horse <- read.table("data/horse.txt", 
                    header = TRUE,
                    stringsAsFactors = TRUE)

l.vals <- seq(from = 0, to = 2, length = 200)
plot(l.vals, dgamma(l.vals, shape = 196.01, rate = 280.01), type = "l", 
     xlab = expression(lambda), ylab = "")
lines(l.vals, dgamma(l.vals, shape = .01, rate = .01), lty = "dashed")

abline(v = 0.7, col = "red")

legend("topleft", leg = c("prior", "posterior"), 
       lty = c("dashed", "solid"))
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-1-1.png" width="672" />

The red line shows the MLE, which is displaced slightly from the posterior mode.

As a point estimate, we might consider any of the following.  The posterior mean can be found exactly as $a/r$ = 0.70001.  Alternatively, we might consider the posterior median

```r
qgamma(0.5, shape = 196.01, rate = 280.01)
```

```
## [1] 0.6988206
```

Finally, we might conisder the posterior mode:

```r
optimize(f = function(x) dgamma(x, shape = 196.01, rate = 280.01), interval = c(0.5, 1), maximum = TRUE)
```

```
## $maximum
## [1] 0.6964383
## 
## $objective
## [1] 7.995941
```

To find a 95\% confidence interval, we might consider the central 95\% interval:

```r
qgamma(c(0.025, 0.975), shape = 196.01, rate = 280.01)
```

```
## [1] 0.6054387 0.8013454
```

A 95\% highest posterior density (HPD) interval takes a bit more work:

```r
diff.in.pdf <- function(x){
  
  upper <- qgamma(p = x, shape = 196.01, rate = 280.01)
  lower <- qgamma(p = x - .95, shape = 196.01, rate = 280.01)
  
  dgamma(upper, shape = 196.01, rate = 280.01) - dgamma(lower, shape = 196.01, rate = 280.01)
}

(upper.qtile <- uniroot(diff.in.pdf, interval = c(0.95, 1))$root)
```

```
## [1] 0.9722176
```

```r
(hpd.ci <- qgamma(p = c(upper.qtile - .95, upper.qtile), shape = 196.01, rate = 280.01))
```

```
## [1] 0.6031732 0.7988576
```

We might also ask questions like: What is the posterior probability that $\lambda > 2/3$?  These caluclations are straightforward in a Bayesian context, and they make full sense.

```r
pgamma(2/3, shape = 196.01, rate = 280.01, lower.tail = FALSE)
```

```
## [1] 0.7434032
```

Thus we would say that there is a 0.743 posterior probability that $\lambda > 2/3$.

As an illustration, note that if we had begun with a more informative prior --- say, a gamma distribution with shape parameter $a = 50$ and rate parameter = $100$ --- then the posterior would have been more of a compromise between the prior and the information in the data:


```r
plot(l.vals, dgamma(l.vals, shape = 196 + 50, rate = 100 + 280), type = "l", 
     xlab = expression(lambda), ylab = "")
lines(l.vals, dgamma(l.vals, shape = 50, rate = 100), lty = "dashed")

abline(v = 0.7, col = "red")

legend("topleft", leg = c("prior", "posterior"), 
       lty = c("dashed", "solid"))
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-7-1.png" width="672" />

## JAGS in R

All of the computational tools that we will examine in this section involve some form of stochastic sampling from the posterior.  This computing companion will largely use the default settings, though in real practice the analyst will often have to do considerable work adjusting the settings to obtain a satisfactory approximation.

We'll use JAGS through R, using the library `r2jags`.  Here is JAGS code to approximate the posterior to $\lambda$ for the horse kick data, using the vague prior.


```r
require(R2jags)
```

```
## Loading required package: R2jags
```

```
## Warning: package 'R2jags' was built under R version 4.1.1
```

```
## Loading required package: rjags
```

```
## Warning: package 'rjags' was built under R version 4.1.1
```

```
## Loading required package: coda
```

```
## Linked to JAGS 4.3.0
```

```
## Loaded modules: basemod,bugs
```

```
## 
## Attaching package: 'R2jags'
```

```
## The following object is masked from 'package:coda':
## 
##     traceplot
```


```r
horse.model <- function() {
  
  for (j in 1:J) {             # J = 280, number of data points
    y[j] ~ dpois (lambda)      # data model:  the likelihood      
  }
  
  lambda ~ dgamma (0.01, 0.01) # prior 
                               # note that BUGS / JAGS parameterizes 
                               # gamma by shape, rate
}

jags.data <- list(y = horse$deaths, J = length(horse$deaths))

jags.params <- c("lambda")

jags.inits <- function(){
  list("lambda" = rgamma(0.01, 0.01))
}

jagsfit <- jags(data               = jags.data, 
                inits              = jags.inits, 
                parameters.to.save = jags.params,
                model.file         = horse.model,
                n.chains           = 3,
                n.iter             = 5000)
```

```
## module glm loaded
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 280
##    Unobserved stochastic nodes: 1
##    Total graph size: 283
## 
## Initializing model
```

Let's take a look at some summary statistics of the fit

```r
print(jagsfit)
```

```
## Inference for Bugs model at "C:/Users/krgross/AppData/Local/Temp/Rtmp6j4Ult/model3e8287b3912.txt", fit using jags,
##  3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
##  n.sims = 3750 iterations saved
##          mu.vect sd.vect    2.5%     25%     50%     75%   97.5%  Rhat n.eff
## lambda     0.699   0.048   0.608   0.667   0.698   0.732   0.795 1.001  3800
## deviance 629.247   1.334 628.310 628.405 628.731 629.542 633.140 1.001  3800
## 
## For each parameter, n.eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
## 
## DIC info (using the rule, pD = var(deviance)/2)
## pD = 0.9 and DIC = 630.1
## DIC is an estimate of expected predictive error (lower deviance is better).
```

The Rhat values suggest that our chains have converged, as we might hope for such a simple model.  We can generate a trace plot using `traceplot` to inspect convergence visually, but beware that visual assessment of convergence is prone to error

We can also use the `lattice` package to construct smoothed estimates of the posterior density:

```r
require(lattice)
```

```
## Loading required package: lattice
```

```r
jagsfit.mcmc <- as.mcmc(jagsfit)
densityplot(jagsfit.mcmc)
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-11-1.png" width="672" />

For a more involved example, let's take a look at the simple regression fit to the cricket data.  First, we'll make a plot of the data and fit a SLR model by least squares.  

```r
cricket <- read.table("data/cricket.txt", header = TRUE)

cricket.slr <- lm(chirps ~ temperature, data = cricket)
summary(cricket.slr)
```

```
## 
## Call:
## lm(formula = chirps ~ temperature, data = cricket)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.56009 -0.57930  0.03129  0.59020  1.53259 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -0.30914    3.10858  -0.099 0.922299    
## temperature  0.21193    0.03871   5.475 0.000107 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.9715 on 13 degrees of freedom
## Multiple R-squared:  0.6975,	Adjusted R-squared:  0.6742 
## F-statistic: 29.97 on 1 and 13 DF,  p-value: 0.0001067
```

```r
plot(chirps ~ temperature, data = cricket)
abline(cricket.slr)
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-12-1.png" width="672" />

Now we'll fit the same model in JAGS, using vague priors for all model parameters

```r
cricket.model <- function() {
  
  for (j in 1:J) {             # J = number of data points
    
    y[j] ~ dnorm (mu[j], tau)  # data model:  the likelihood 
                               # note that BUGS / JAGS uses precision
                               # instead of variance
    
    mu[j] <- b0 + b1 * x[j]    # compute the mean for each observation
  }
  
  b0 ~ dnorm (0.0, 1E-6)       # prior for intercept
  b1 ~ dnorm (0.0, 1E-6)       # prior for slope
  tau ~ dgamma (0.01, 0.01)    # prior for tau
                               # note that BUGS / JAGS parameterizes 
                               # gamma by shape, rate
  
  sigma <- pow(tau, -1/2)      # the SD of the residaul errors
}

jags.data <- list(y = cricket$chirps, 
                  x = cricket$temperature,
                  J = nrow(cricket))

jags.params <- c("b0", "b1", "tau", "sigma")

jags.inits <- function(){
  list("b0" = rnorm(1), "b1" = rnorm(1), "tau" = runif(1))
}

jagsfit <- jags(data               = jags.data, 
                inits              = jags.inits, 
                parameters.to.save = jags.params,
                model.file         = cricket.model,
                n.chains           = 3,
                n.iter             = 5000)
```

```
## module glm loaded
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 15
##    Unobserved stochastic nodes: 3
##    Total graph size: 70
## 
## Initializing model
```

```r
print(jagsfit)
```

```
## Inference for Bugs model at "C:/Users/krgross/AppData/Local/Temp/RtmpaALYhr/modele205b6f45ab.txt", fit using jags,
##  3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
##  n.sims = 3750 iterations saved
##          mu.vect sd.vect   2.5%    25%    50%    75%  97.5%  Rhat n.eff
## b0        -0.267   3.350 -6.904 -2.385 -0.289  1.802  6.401 1.001  3800
## b1         0.211   0.042  0.128  0.185  0.212  0.238  0.293 1.001  3800
## sigma      1.037   0.229  0.704  0.877  1.000  1.155  1.596 1.001  3800
## tau        1.056   0.421  0.392  0.749  1.000  1.300  2.019 1.001  3800
## deviance  42.940   2.773 39.796 40.925 42.243 44.144 50.139 1.002  2200
## 
## For each parameter, n.eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
## 
## DIC info (using the rule, pD = var(deviance)/2)
## pD = 3.8 and DIC = 46.8
## DIC is an estimate of expected predictive error (lower deviance is better).
```

```r
traceplot(jagsfit)
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-13-1.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-13-2.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-13-3.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-13-4.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-13-5.png" width="672" />
The traces for the intercept aren't great, but we haven't centered the predictor either.  In the usual way, the slope and intercept are strongly negatively correlated in the posterior.  We can visualize this posterior correlation:


```r
library(hexbin)
```

```
## Warning: package 'hexbin' was built under R version 4.1.1
```

```r
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
with(jagsfit$BUGSoutput$sims.list, hexbinplot(b0 ~ b1, colramp = rf))
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-14-1.png" width="672" />



We could make life easier on ourselves by centering the predictor and trying again:


```r
cricket$temp.ctr <- cricket$temperature - mean(cricket$temperature)

jags.data <- list(y = cricket$chirps, 
                  x = cricket$temp.ctr,
                  J = nrow(cricket))

jags.params <- c("b0", "b1", "tau", "sigma")

jags.inits <- function(){
  list("b0" = rnorm(1), "b1" = rnorm(1), "tau" = runif(1))
}

jagsfit <- jags(data               = jags.data, 
                inits              = jags.inits, 
                parameters.to.save = jags.params,
                model.file         = cricket.model,
                n.chains           = 3,
                n.iter             = 5000)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 15
##    Unobserved stochastic nodes: 3
##    Total graph size: 70
## 
## Initializing model
```

```r
print(jagsfit)
```

```
## Inference for Bugs model at "C:/Users/krgross/AppData/Local/Temp/RtmpaALYhr/modele20767654a1.txt", fit using jags,
##  3 chains, each with 5000 iterations (first 2500 discarded), n.thin = 2
##  n.sims = 3750 iterations saved
##          mu.vect sd.vect   2.5%    25%    50%    75%  97.5%  Rhat n.eff
## b0        16.656   0.279 16.097 16.484 16.653 16.830 17.217 1.001  3800
## b1         0.212   0.042  0.130  0.185  0.211  0.238  0.297 1.002  1200
## sigma      1.038   0.224  0.712  0.879  1.002  1.158  1.570 1.001  3800
## tau        1.050   0.415  0.406  0.745  0.996  1.295  1.974 1.001  3800
## deviance  42.933   2.711 39.823 40.945 42.249 44.175 50.133 1.001  3800
## 
## For each parameter, n.eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
## 
## DIC info (using the rule, pD = var(deviance)/2)
## pD = 3.7 and DIC = 46.6
## DIC is an estimate of expected predictive error (lower deviance is better).
```

```r
traceplot(jagsfit)
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-15-1.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-15-2.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-15-3.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-15-4.png" width="672" /><img src="03-BayesianComputation_files/figure-html/unnamed-chunk-15-5.png" width="672" />

The posteriors for the intercept and slope are now uncorrelated:


```r
library(hexbin)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
with(jagsfit$BUGSoutput$sims.list, hexbinplot(b0 ~ b1, colramp = rf))
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-16-1.png" width="672" />

## Stan

Stan is a separate program based on Hamiltonian Monte Carlo.  Stan can be accesses through R using the `rstan` library.

We'll start by looking at how to use `rstan` to estimate $\lambda$ for the horse-kick data using vague priors.  `rstan` requires that a Stan program be prepared as a text file and stored locally.  To fit the horse-kick data, we'll use the Stan program listed here.  This program is saved as the file `horse.stan`.
```stan
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
  // prior
  lambda ~ gamma(0.01, 0.01);
  
  // data model
  y ~ poisson(lambda);
}
```


```r
require(rstan)
options(mc.cores = parallel::detectCores())
```


```r
horse.dat <- list(N = nrow(horse),
                  y = horse$deaths)

stan.fit <- stan(file = 'stan/horse.stan', data = horse.dat)
```

Apparently this warning is not particularly problematic.  We can have a look at the fit by asking to `print` it:

```r
print(stan.fit)
```

```
## Inference for Stan model: horse.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##           mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
## lambda    0.70    0.00 0.05    0.61    0.67    0.70    0.73    0.80  1515    1
## lp__   -266.41    0.02 0.71 -268.43 -266.57 -266.12 -265.96 -265.92  1997    1
## 
## Samples were drawn using NUTS(diag_e) at Wed Sep 08 21:52:31 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Next, we'll use Stan to fit the linear regression model.  Again, we need to write a Stan program that R will call. The Stan program for the regression model is shown below.

```stan
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
```


```r
cricket.dat <- list(N = nrow(cricket),
                    y = cricket$chirps,
                    x = cricket$temp.ctr)

stan.fit <- stan(file = 'stan/cricket.stan', data = cricket.dat)
print(stan.fit)
```

```
## Inference for Stan model: cricket.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##        mean se_mean   sd   2.5%   25%   50%   75% 97.5% n_eff Rhat
## b0    16.64    0.00 0.26  16.11 16.48 16.65 16.81 17.15  3194    1
## b1     0.21    0.00 0.04   0.13  0.19  0.21  0.24  0.29  3357    1
## sigma  1.02    0.00 0.21   0.70  0.87  0.99  1.14  1.54  2103    1
## lp__  -8.96    0.03 1.26 -12.28 -9.56 -8.68 -8.02 -7.51  1354    1
## 
## Samples were drawn using NUTS(diag_e) at Wed Sep 08 22:19:07 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

We can make a plot of the posterior samples using `pairs`:

```r
pairs(stan.fit, pars = c("b0", "b1", "sigma"))
```

<img src="03-BayesianComputation_files/figure-html/unnamed-chunk-21-1.png" width="672" />



