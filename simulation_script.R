### Simulation 1: 
### Gaussian prior, into gaussian distribution 

reps <- 1e4
n <- 1e4

true_mu <- 0
sigma <- 1

### Allocate memory for error vectors 
squared_error <- numeric(reps)
absolute_error <- numeric(reps)

for (i in 1:reps){
  theta <- rnorm(n, true_mu, sigma)
  norm_vec <- rnorm(n, mean = theta)
  mu_hat <- mean(norm_vec)
  sig_hat_m <- (length(norm_vec)^(-1))*(sum((norm_vec - mu_hat)^2))
  
  ### Wiki source
  sig_hat_pi <- sig_hat_m - 1
  mu_hat_pi <- mu_hat
  
  sq_err <- (true_mu - mu_hat_pi)^2
  squared_error[i] <- sq_err
  
  abs_err <- (true_mu - mu_hat_pi)
  absolute_error[i] <- abs_err
}
mean(absolute_error)
mean(squared_error)

###############################################################################
###############################################################################

### Simulation 2
### Beta Prior, Bernoulli distribution (beta posterior)
reps <- 1e4
n <- 1e4

### Allocate memory for error vectors 
squared_error <- numeric(reps)
absolute_error <- numeric(reps)

alpha <- 2
beta <- 2 
true_mean <- alpha/(alpha + beta)

for (i in 1:reps){
  theta <- rbeta(n, alpha, beta)
  x_real <- rbinom(n = n, prob = theta, size = 1)
  post_mean <- (sum(x_real) + 2)/(length(x_real) + 2 + 2)
  
  ### Calc Squared Error and append 
  sq_err <- (true_mean - post_mean)^2
  squared_error[i] <- sq_err
  
  
  ## Calculate absolute error and append 
  abs_err <- abs(true_mean - post_mean)
  absolute_error <- abs_err
}
mean(squared_error)
mean(absolute_error)

