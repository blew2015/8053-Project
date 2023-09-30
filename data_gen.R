### geom_gen function generates realizations from a Bernoulli distribution
### using probabilities from an input vector 
### 
### Parameters:
###     p: vector of probabilities 
###     n: number of realizations to generate from the Bernoulli distribution 
###        for each probability p
###
###
### Output: 
###     bern_mat: Matrix with n+2 columns and length(p) rows. One row for each
###         value in the p input vector. The first n columns contain realizations
###         from the Bernoulli distribution. The second to last and last columns 
###         contain the true value of p and the MLE of p (row mean of realizations).
###         respectively 

bern_gen <- function(p, n){
  # Allocate memory for matrix
  bern_mat <- matrix(NA, nrow = length(p), ncol = n + 2)
  
  # Loop of values of p
  for (i in 1:length(p)){
    # Set first n columns of ith rows to the Bernoulli realizations
    bern_mat[i,1:n] <- rbinom(n = n, prob = p[i], size = 1)
  }
  
  # Calculate MLE estimator of p
  bern_mat[,n+2] <- rowMeans(bern_mat, na.rm = TRUE) 
  
  # save true mean
  bern_mat[,n+1] <- p
  return(bern_mat)
}


### norm_gen function generates realizations from a normal distribution
### using mean/variance combinations from two input vectors
### 
### Parameters:
###     mu: vector of means 
###     sigma2: vector of variances
###     n: number of realizations to generate from the normal distribution 
###        for each mu/sigma2 pair
###
###
### Output: 
###     norm_mat: Matrix with n+4 columns and length(mu) rows. One row for each
###         value in the mu input vector. The first n columns contain realizations
###         from the N(mu, sigma2) distribution. The last four columns 
###         contain the true value of both mu and sigma and the MLE of mu and 
###         sigma respectively d
###         (n+1 col = true Mu)
###         (n+2 col = sample mean)
###         (n+3 col = true value of sigma2)
###         (n+4 col = S^2)

norm_gen <- function(mu, sigma2, n){
  mu <- rnorm(15)
  sigma2 <- rep(1, 15)
  # Allocate memory for matrix
  norm_mat <- matrix(NA, nrow = length(mu), ncol = n + 4)
  
  # Loop of values of p
  for (i in 1:length(mu)){
    # Set first n columns of ith rows to the Bernoulli realizations
    norm_mat[i,1:n] <- rnorm(n = n, mean = mu[i], sd = sqrt(sigma2[i]))
    
    # Calculate MLE estimator of p
    norm_mat[i,n+2] <- mean(norm_mat[i,], na.rm = TRUE) 
    
    # save true mean
    norm_mat[i,n+1] <- mu[i]
    
    # MLE estimator of sigma^2 (double check)
    norm_mat[i,n+4] <-  (1/(n+1))*sum((norm_mat[i,1:n] - norm_mat[i,n+2])^2, na.rm = TRUE)
    
    # True value of sigma^2
    norm_mat[i, n+3] <- sigma2[i]
    
  }
  
  return(norm_mat)
  
}



set.seed(9377)

### Bernoulli test (beta prior)
n = 1e4
reps = 30
p = rbeta(n, 3, 1)
bern_gen(n = reps, p)


### Normal test (normal prior (mean))
mu = rnorm(n, 0, 1)
sigma2 <- rep(10, n)
norm_gen(mu, sigma2, n = reps)



