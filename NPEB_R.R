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