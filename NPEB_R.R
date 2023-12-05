#install.packages("sfsmisc")
library(sfsmisc) # for marginal density estimation
library(ks) # for kde and derivative of kde

# Simple normal case

X <- seq(-5,5, length = 1000)
set.seed(1234)
X_norm <- rnorm(1000)

# Density estimate
tkdensity(X_norm) # Gives you optimal kernel and bandwidth

fx <- density(X_norm, bw=0.217, kernel = "gaussian")
plot(fx)
fx_d <- approxfun(fx)
ss10 <- smooth.spline(fx$y~fx$x , df = 10)

X_norm_ordered <- sort(X_norm)
fx <- kde(X_norm_ordered, h=0.217, eval.points = X_norm_ordered)
fx_d <- kdde(X_norm_ordered, h=0.217, deriv.order = 1, eval.points = X_norm_ordered)

par(mfrow=c(2,1))
plot(fx)
plot(fx_d)

score <- fx_d$estimate/fx$estimate

plot(score)

## Tweedie's formula (E(\theta|X) = x+sigma^2*score(x))
theta_hat <- mean(fx$eval.points+score)

theta_hat_raw <- fx$eval.points+score
plot(x=fx$eval.points, y=theta_hat_raw, ylim = c(-10,10))

## MSE
MSE_NPEB_normal <- mean((theta_hat_raw)^2)


## Simulations - 1) Different sample sizes from N=10 to N=10000
MSE_Normal_N <- c()
N = 10000 
for (i in 10:N) {
  set.seed(1234)
  X_norm <- rnorm(i) # N(0,1)
  X_norm_ordered <- sort(X_norm)
  fx <- kde(X_norm_ordered, eval.points = X_norm_ordered)
  fx_d <- kdde(X_norm_ordered, deriv.order = 1, eval.points = X_norm_ordered)
  score <- fx_d$estimate/fx$estimate
  
  theta_hat_raw <- fx$eval.points+score
  MSE_Normal_N[i-9] <- mean((theta_hat_raw)^2)
}

plot(MSE_Normal_N, type="l", xlab="N", ylab = "MSE", main="MSE of Nonparametric Empirical Bayes - N(0,1)")

## Simulations - 2) Different mu values (-10 to 10)
MSE_Normal_mu <- c()
N = 5000 
mu <- seq(-10, 10, length=10000)
for (i in 1:10000) {
  set.seed(1234)
  X_norm <- rnorm(N, mean=mu[i]) # N(i,1)
  X_norm_ordered <- sort(X_norm)
  fx <- kde(X_norm_ordered, eval.points = X_norm_ordered)
  fx_d <- kdde(X_norm_ordered, deriv.order = 1, eval.points = X_norm_ordered)
  score <- fx_d$estimate/fx$estimate
  
  theta_hat_raw <- fx$eval.points+score
  MSE_Normal_mu[i] <- mean((theta_hat_raw-mu[i])^2)
}

plot(MSE_Normal_mu~mu, type="l", xlab=expression(mu), ylab = "MSE", main="MSE of Nonparametric Empirical Bayes - N(mu,1)", ylim=c(0,0.01))


## Simulations - 2) Different sigma values (1 to 10)
MSE_Normal_sigma <- c()
N = 5000 
sig <- seq(1, 10, length=10000)
for (i in 1:10000) {
  set.seed(1234)
  X_norm <- rnorm(N, sd=sig[i]) # N(0,i)
  X_norm_ordered <- sort(X_norm)
  fx <- kde(X_norm_ordered, eval.points = X_norm_ordered)
  fx_d <- kdde(X_norm_ordered, deriv.order = 1, eval.points = X_norm_ordered)
  score <- fx_d$estimate/fx$estimate
  
  theta_hat_raw <- fx$eval.points+score
  MSE_Normal_sigma[i] <- mean((theta_hat_raw)^2)
}

plot(MSE_Normal_sigma~sig, type="l", xlab=expression(sigma), ylab = "MSE", main="MSE of Nonparametric Empirical Bayes - N(mu,1)", ylim=c(0,0.01))


## Poisson case (x+1) * f(x+1)/f(x)

X_pos <- seq(0,10)
set.seed(1234)
X_poi <- rpois(1000, lambda = 1)

# f(x) can be estimated using y_x/N and f(x+1) can be estimated using y_{x+1}/N

poi <- as.data.frame(table(X_poi))
poi$X_poi <- as.numeric(poi$X_poi) - 1

lambda <- c()
for (i in 1:nrow(poi)-1){
  lambda[i] <- (poi$X_poi[i]+1) * (poi$Freq[i+1]) / (poi$Freq[i])
}

lambda
plot(lambda)

## MSE
MSE_NPEB_poi <- mean((lambda-1)^2)

## Simulations - 1) Different sample sizes from N=10 to N=10000
MSE_Poisson_N <- c()
N = 10000 
for (i in 10:N) {
  set.seed(1234)
  X_poi <- rpois(i, lambda = 1)
  poi <- as.data.frame(table(X_poi))
  poi$X_poi <- as.numeric(poi$X_poi) - 1
  
  lambda <- c()
  for (j in 1:nrow(poi)-1){
    lambda[j] <- (poi$X_poi[j]+1) * (poi$Freq[j+1]) / (poi$Freq[j])
  }
  MSE_Poisson_N[i-9] <- mean((lambda-1)^2)
}

plot(MSE_Poisson_N, type="l", xlab="N", ylab = "MSE", main="MSE of Nonparametric Empirical Bayes - Poisson(1)")

# Simulations - 2) Different lambda values (from 1 to 100), N=5000
MSE_Poisson_lambda <- c()
N=5000
for (i in 1:100) {
  set.seed(1234)
  X_poi <- rpois(N, lambda = i)
  poi <- as.data.frame(table(X_poi))
  poi$X_poi <- as.numeric(poi$X_poi) - 1
  lambda <- c()
  for (j in 1:nrow(poi)-1){
    lambda[j] <- (poi$X_poi[j]+1) * (poi$Freq[j+1]) / (poi$Freq[j])
  }
  MSE_Poisson_lambda[i] <- mean((lambda-i)^2)
}

plot(MSE_Poisson_lambda, type="l", xlab=expression(lambda), ylab = "MSE", main="MSE of Nonparametric Empirical Bayes - Poisson")