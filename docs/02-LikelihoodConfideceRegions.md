# Beyond the MLE: Confidence regions and hypothesis tests using the likelihood function

Likelihood can be used for more than simply isolating the MLE.  The likelihood can also be used to generate confidence intervals for single parameters, or confidence regions for several parameters.  We'll start by using the horse-kick data to see how to generate a confidence interval for a single parameter, and then move on to considering models with more than one parameter.

## Confidence intervals for single parameters

First we'll read in the data and recreate the negative log likelihood function.


```r
#################
## Preparation
################

# read in the data

horse <- read.table("data/horse.txt", header = TRUE)

horse.neg.ll <- function(my.lambda) {
  ll.vals <- dpois(x = horse$deaths, lambda = my.lambda, log = TRUE)
  -1 * sum(ll.vals) 
}

# create a vector of lambda values using the 'seq'uence command
lambda.vals <- seq(from = 0.5, to = 1.0, by = 0.01)

# create an empty vector to store the values of the log-likelihood
ll.vals <- double(length = length(lambda.vals))

# use a loop to find the log-likelihood for each value in lambda.vals
for (i.lambda in 1:length(lambda.vals)) {
  ll.vals[i.lambda] <- horse.neg.ll(lambda.vals[i.lambda])
}

plot(ll.vals ~ lambda.vals, xlab = "lambda", ylab = "negative log likelihood", type = "l")
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-1-1.png" width="672" />

Now we'll find the a (asymptotic) 95\% confidence interval for $\lambda$ directly.

```r
cutoff.ll <- horse.neg.ll(0.7) + qchisq(0.95, df = 1) / 2

# recreate the plot and add a line
plot(ll.vals ~ lambda.vals, xlab = "lambda", ylab = "negative log likelihood", type = "l")
abline(h = cutoff.ll, col = "red", lty = "dashed")
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-2-1.png" width="672" />

```r
# use uniroot to find the confidence bounds precisely

my.function <- function(my.lambda){
  
  horse.neg.ll(0.7) + qchisq(0.95, df = 1) / 2 - horse.neg.ll(my.lambda)  
}

(lower <- uniroot(f = my.function, interval = c(0.6, 0.7)))
```

```
## $root
## [1] 0.6065198
## 
## $f.root
## [1] -3.556854e-05
## 
## $iter
## [1] 4
## 
## $init.it
## [1] NA
## 
## $estim.prec
## [1] 6.103516e-05
```

```r
(upper <- uniroot(f = my.function, interval = c(0.7, 0.9)))
```

```
## $root
## [1] 0.8026265
## 
## $f.root
## [1] -0.0001007316
## 
## $iter
## [1] 6
## 
## $init.it
## [1] NA
## 
## $estim.prec
## [1] 6.103516e-05
```

As an alternative programming style, we could have defined the objective function on the fly, and not bothered to create `my.function`.


```r
(lower <- uniroot(f = function(x) horse.neg.ll(0.7) + qchisq(0.95, df = 1) / 2 - horse.neg.ll(x) ,
                 interval = c(0.6, 0.7)))
```

```
## $root
## [1] 0.6065198
## 
## $f.root
## [1] -3.556854e-05
## 
## $iter
## [1] 4
## 
## $init.it
## [1] NA
## 
## $estim.prec
## [1] 6.103516e-05
```

Let's recreate the plot and add vertical lines to indicate the confidence interval.


```r
plot(ll.vals ~ lambda.vals, xlab = "lambda", ylab = "negative log likelihood", type = "l")
abline(h = cutoff.ll, col = "red", lty = "dashed")
abline(v = c(lower$root, upper$root), col = "red")
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-4-1.png" width="672" />

```r
# clean up the workspace
rm(list = ls())
```
Thus, the 95\% CI for $\lambda$ is $(0.607, 0.803)$.

## Confidence regions, profile likelihoods, and associated univariate intervals

With a 2-parameter model, we can plot a confidence region directly.  First some housekeeping to get started:


```r
library(emdbook)
data("ReedfrogFuncresp")

# rename something shorter

frog <- ReedfrogFuncresp
rm(ReedfrogFuncresp)

frog.neg.ll <- function(params){
  
  a <- params[1]
  h <- params[2]
  
  prob.vals <- a / (1 + a * h * frog$Initial)
  
  ll.vals <- dbinom(frog$Killed, size = frog$Initial, prob = prob.vals, log = TRUE)
  -1 * sum(ll.vals)
}

(frog.mle <- optim(par = c(0.5, 1/60),
                   fn  = frog.neg.ll))
```

```
## Warning in dbinom(frog$Killed, size = frog$Initial, prob = prob.vals, log =
## TRUE): NaNs produced
```

```
## $par
## [1] 0.52585566 0.01660104
## 
## $value
## [1] 46.72136
## 
## $counts
## function gradient 
##       61       NA 
## 
## $convergence
## [1] 0
## 
## $message
## NULL
```

```r
a.mle <- frog.mle$par[1]
h.mle <- frog.mle$par[2]

# plot negative likelihood contours

a.vals <- seq(from = 0.3, to = 0.75, by = 0.01)
h.vals <- seq(from = 0.001, to = 0.03, by = 0.001)

ll.vals <- matrix(nrow = length(a.vals), ncol = length(h.vals))

for (i.a in 1:length(a.vals)) {
  for(i.h in 1:length(h.vals)) {
    ll.vals[i.a, i.h] <- frog.neg.ll(c(a.vals[i.a], h.vals[i.h]))
  }
}

contour(x = a.vals, y = h.vals, z = ll.vals, nlevels = 100,
        xlab = "a", ylab = "h")

points(x = a.mle, y = h.mle, col = "red")
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-5-1.png" width="672" />

Equipped with the contour plot, graphing the appropriate confidence region is straightforward.


```r
cut.off <- frog.neg.ll(c(a.mle, h.mle)) + (1 / 2) * qchisq(.95, df = 2)

# recreate the plot and add a line for the 95% confidence region
contour(x = a.vals, y = h.vals, z = ll.vals, nlevels = 100,
        xlab = "a", ylab = "h")

points(x = a.mle, y = h.mle, col = "red")
contour(x = a.vals, y = h.vals, z = ll.vals, 
        levels = cut.off,
        add = TRUE, col = "red", lwd = 2)
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-6-1.png" width="672" />

However, there are several drawbacks to confidence regions.  First, while a two-dimensional confidence region can be readily visualized, it is hard to summarize or describe.  Second, and more importantly, most models have more than two parameters.  In these models, a confidence region would have more than 2 dimensions, and thus would be impractical to visualize.  Thus it is helpful, or even essential, to be able to generate univariate confidence intervals for single parameters from high-dimensional likelihoods.  One approach to doing so is to calculate the so-called profile likelihood for a given parameter, and then to derive the univariate interval from this profile likelihood.  We will illustrate this approach by computing a profile-based confidence interval for the attack rate $a$ in the tadpole data.


```r
# profile log-likelihood function for the attack rate a

profile.ll <- function(my.a) {
  
  # Calculate the minimum log likelihood value for a given value of a, the attack rate
  
  my.ll <- function(h) frog.neg.ll(c(my.a, h))
  my.profile <- optimize(f = my.ll, interval = c(0, 0.03), maximum = FALSE)
  
  my.profile$objective
}

# plot the profile likelihood vs. a
# not necessary for finding the CI, but useful for understanding

a.values <- seq(from = 0.3, to = 0.8, by = 0.01)
a.profile <- double(length = length(a.values))

for (i in 1:length(a.values)) {
  
  a.profile[i] <- profile.ll(a.values[i])
}

plot(x = a.values, y = a.profile, xlab = "a", ylab = "negative log-likelihood", type = "l")
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-7-1.png" width="672" />

Now we'll follow the same steps as before to compute the profile-based 95% CI.


```r
# Now follow the same steps as before to find the profile 95% CI

cut.off <- profile.ll(a.mle) + qchisq(0.95, df = 1) / 2

(lower <- uniroot(f = function(x) cut.off - profile.ll(x) ,
                  interval = c(0.3, a.mle)))
```

```
## $root
## [1] 0.4024268
## 
## $f.root
## [1] -0.0001303772
## 
## $iter
## [1] 6
## 
## $init.it
## [1] NA
## 
## $estim.prec
## [1] 6.103516e-05
```

```r
(upper <- uniroot(f = function(x) cut.off - profile.ll(x) ,
                  interval = c(a.mle, 0.8)))
```

```
## $root
## [1] 0.6824678
## 
## $f.root
## [1] -9.763258e-06
## 
## $iter
## [1] 6
## 
## $init.it
## [1] NA
## 
## $estim.prec
## [1] 6.103516e-05
```

```r
plot(x = a.values, y = a.profile, xlab = "a", ylab = "negative log-likelihood", type = "l")
abline(v = c(lower$root, upper$root), col = "blue")
abline(h = cut.off, col = "blue", lty = "dashed")
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-8-1.png" width="672" />

So, the 95\% profile CI for $a$ is (0.402, 0.682).

## Locally quadratic approximations to confidence intervals and regions

Likelihood profiling provides a straightforward way to summarize a high-dimensional confidence region by univariate confidence intervals.  However, these profile intervals can still involve quite a bit of computation.  Further, they are not able to capture possible correlations among parameter estimates, which are revealed in (two-dimensional) confidence regions.  (Recall the shape of the joint confidence region for the parameters $a$ and $h$ in the tadpole data.)  Locally quadratic approximations provide a way to approximate the (already approximate) univariate confidence intervals and bivariate confidence regions using only information about the curvature of the likelihood surface at the MLE.

We'll first start by revisiting the horse-kick data again.  Of course, with the more precise $\chi^2$ based confidence interval in hand, there is no reason to seek an approximation.  But doing so allows us to illustrate the calculations involved, and to see how well the approximation fares in this case.

First some housekeepign to read the data into memory, etc.


```r
# clean up
rm(list = ls())

# read in the data

horse <- read.table("data/horse.txt", header = TRUE)

horse.neg.ll <- function(my.lambda) {
  ll.vals <- dpois(x = horse$deaths, lambda = my.lambda, log = TRUE)
  -1 * sum(ll.vals) 
}

# use uniroot to find the confidence bounds precisely

my.function <- function(my.lambda){
  
  horse.neg.ll(0.7) + qchisq(0.95, df = 1) / 2 - horse.neg.ll(my.lambda)  
}

lower <- uniroot(f = my.function, interval = c(0.6, 0.7))
upper <- uniroot(f = my.function, interval = c(0.7, 0.9))
```

Now we will proceed to use a locally quadratic approximation to the negative log likelihood.

```r
## this function finds the second derivative at the MLE by finite differences

second.deriv <- function(delta.l) {
  
  (horse.neg.ll(0.7 + delta.l) - 2 * horse.neg.ll(0.7) + horse.neg.ll(0.7 - delta.l)) / delta.l ^ 2
}

(horse.D2 <- second.deriv(1e-04))
```

```
## [1] 400
```

```r
# see how the answer changes if we change delta
second.deriv(1e-05)
```

```
## [1] 400.0003
```

Let's compare this answer to the answer obtained by `numDeriv::hessian`.

```r
numDeriv::hessian(func = horse.neg.ll, x = 0.7)
```

```
##      [,1]
## [1,]  400
```

The approximate standard error of $\hat{\lambda}$ is the square root of the inverse of the second derivative of the likelihood function.


```r
(lambda.se <- sqrt(1 / horse.D2))
```

```
## [1] 0.05
```

Now we can approximate the 95\% confidence interval by using critical values from a standard normal distribution.

```r
(lower.approx <- 0.7 - qnorm(.975) * lambda.se)
```

```
## [1] 0.6020018
```

```r
(upper.approx <- 0.7 + qnorm(.975) * lambda.se)
```

```
## [1] 0.7979982
```

Compare the approximation to the "exact" values

```r
lower$root
```

```
## [1] 0.6065198
```

```r
upper$root
```

```
## [1] 0.8026265
```

Make a plot

```r
# create a vector of lambda values using the 'seq'uence command
lambda.vals <- seq(from = 0.5, to = 1.0, by = 0.01)

# create an empty vector to store the values of the log-likelihood
ll.vals <- double(length = length(lambda.vals))

# use a loop to find the log-likelihood for each value in lambda.vals
for (i.lambda in 1:length(lambda.vals)) {
  ll.vals[i.lambda] <- horse.neg.ll(lambda.vals[i.lambda])
}

plot(ll.vals ~ lambda.vals, xlab = "lambda", ylab = "negative log likelihood", type = "l")

###################################
## Now find the confidence interval, and plot it
####################################

cutoff.ll <- horse.neg.ll(0.7) + qchisq(0.95, df = 1) / 2

# add a line to the plot
abline(h = cutoff.ll, col = "red", lty = "dashed")
abline(v = c(lower$root, upper$root), col = "red")
abline(v = c(lower.approx, upper.approx), col = "blue")

legend(x = 0.65, y = 326, 
       leg = c("exact", "approximate"), 
       pch = 16, 
       col = c("red", "blue"),
       bty = "n")
```

<img src="02-LikelihoodConfideceRegions_files/figure-html/unnamed-chunk-15-1.png" width="672" />

```r
# clean up
rm(list = ls())
```

Notice that the full $\chi^2$-based confidence intervals capture the asymmetry in the information about $\lambda$.  The intervals based on the quadratic approximation are symmetric.

Now, use the quadratic approximation to find standard errors for $\hat{a}$ and $\hat{h}$ in the tadpole predation data.  

The first part is preparatory work from old classes.

```r
library(emdbook)
data("ReedfrogFuncresp")

# rename something shorter

frog <- ReedfrogFuncresp
rm(ReedfrogFuncresp)

frog.neg.ll <- function(params){
  
  a <- params[1]
  h <- params[2]
  
  prob.vals <- a / (1 + a * h * frog$Initial)
  
  ll.vals <- dbinom(frog$Killed, size = frog$Initial, prob = prob.vals, log = TRUE)
  -1 * sum(ll.vals)
}

frog.mle <- optim(par = c(0.5, 1/60),
                  fn  = frog.neg.ll)
```

```
## Warning in dbinom(frog$Killed, size = frog$Initial, prob = prob.vals, log =
## TRUE): NaNs produced
```

```r
(a.mle <- frog.mle$par[1])
```

```
## [1] 0.5258557
```

```r
(h.mle <- frog.mle$par[2])
```

```
## [1] 0.01660104
```

Now find the hessian:

```r
(D2 <- numDeriv::hessian(func = frog.neg.ll, x = c(a.mle, h.mle)))
```

```
##            [,1]       [,2]
## [1,]   616.5606  -7394.263
## [2,] -7394.2628 130640.685
```

The matrix inverse of the hessian is the variance-covariance matrix of the parameters.  Note that R uses the function `solve` to find the inverse of a matrix.

```r
# invert to get var-cov matrix
(var.matrix <- solve(D2))
```

```
##              [,1]         [,2]
## [1,] 0.0050493492 2.857932e-04
## [2,] 0.0002857932 2.383048e-05
```

We can use the handy `cov2cor` function to convert the variance matrix into a correlation matrix:

```r
cov2cor(var.matrix)
```

```
##           [,1]      [,2]
## [1,] 1.0000000 0.8238872
## [2,] 0.8238872 1.0000000
```

Note the large correlation between $\hat{a}$ and $\hat{h}$.  Compare with Figure 6.13 of Bolker.

The standard errors of $\hat{a}$ and $\hat{h}$ are the square roots of the diagaonal elements of the variance-covariance matrix.

```r
(a.se <- sqrt(var.matrix[1, 1]))
```

```
## [1] 0.07105877
```

```r
(h.se <- sqrt(var.matrix[2, 2]))
```

```
## [1] 0.004881647
```
Note the large correlation between $\hat{a}$ and $\hat{h}$.

Let's use the (approximate) standard error of $\hat{a}$ to calculate an (approximate) 95\% confidence interval:

```r
(ci.approx <- a.mle + qnorm(c(0.025, .975)) * a.se)
```

```
## [1] 0.3865830 0.6651283
```
Recall that the 95\% confidence interval we calculated by the profile likelihood was $(0.402, 0.682)$.  So the quadratic approximation has gotten the width of the interval more or less correct, but it has fared less at capturing the asymmetry of the interval.
