set.seed(8053)
library(zoo)

# ### Simulation 1: 
# ### Gaussian prior, into gaussian distribution 
# 
# reps <- 1e4
# n <- 1e4
# 
# true_mu <- 0
# sigma <- 1
# 
# ### Allocate memory for error vectors 
# squared_error <- numeric(reps)
# absolute_error <- numeric(reps)
# 
# for (i in 1:reps){
#   theta <- rnorm(n, true_mu, sigma)
#   norm_vec <- rnorm(n, mean = theta)
#   mu_hat <- mean(norm_vec)
#   sig_hat_m <- (length(norm_vec)^(-1))*(sum((norm_vec - mu_hat)^2))
#   
#   ### Wiki source
#   sig_hat_pi <- sig_hat_m - 1
#   mu_hat_pi <- mu_hat
#   
#   sq_err <- (true_mu - mu_hat_pi)^2
#   squared_error[i] <- sq_err
#   
#   abs_err <- (true_mu - mu_hat_pi)
#   absolute_error[i] <- abs_err
# }
# mean(absolute_error)
# mean(squared_error)
# 
# ###############################################################################
# ###############################################################################
# 
# ### Simulation 2
# ### Beta Prior, Bernoulli distribution (beta posterior)
# reps <- 1e4
# n <- 100
# 
# ### Allocate memory for error vectors 
# squared_error <- numeric(reps)
# 
# alpha <- 2
# beta <- 2 
# true_mean <- alpha/(alpha + beta)
# 
# ### Estimation for parames of Beta Binomial Distribution Source: 
# ### https://www-jstor-org.ezp1.lib.umn.edu/stable/2529131?searchText=Maximum+Likelihood+Estimation+for+the+Beta-Binomial+Distribution+and+an+Application+to+the+Household+Distribution+of+the+Total+Number+of+Cases+of+a+Disease&searchUri=%2Faction%2FdoBasicSearch%3FQuery%3DMaximum%2BLikelihood%2BEstimation%2Bfor%2Bthe%2BBeta-Binomial%2BDistribution%2Band%2Ban%2BApplication%2Bto%2Bthe%2BHousehold%2BDistribution%2Bof%2Bthe%2BTotal%2BNumber%2Bof%2BCases%2Bof%2Ba%2BDisease&ab_segments=0%2Fbasic_expensive_solr_cloud%2Fcontrol&refreqid=fastly-default%3A45b6d263c009119688049c66cd7a1e2c&seq=2
# 
# for (i in 1:reps){
#   theta <- rbeta(n, alpha, beta)
#   x_real <- rbinom(n = n, prob = theta, size = 1)
#   
#   indexes <- sample(seq(1, length(theta), by = 1), size = 30)
#   
#   train_theta <- theta[indexes]
#   train_x_real <- x_real[indexes]
#   
#   test_x_real <- x_real[-indexes]
#   test_theta <- theta[-indexes]
#   
#   post_mean_tru_bayes <- (sum(train) + 2)/(length(x_real) + 2 + 2)
#   
#   post_mean_bayes
#   
#   ### Calc Squared Error and append 
#   sq_err <- (true_mean - post_mean)^2
#   squared_error[i] <- sq_err
#   
#   ## Calculate absolute error and append 
#   abs_err <- abs(true_mean - post_mean)
#   absolute_error <- abs_err
# }
# mean(squared_error)
# mean(absolute_error)
# 


### Poisson distribution with Gamma Prior:

reps <- 50

pois_gam_sim<- function(n, alpha, beta, fak_alpha, fak_beta){
  n <- n
  
  squared_error <- matrix(NA, ncol = 3, nrow = n)
  
  tru_alpha <- alpha
  tru_beta <- beta
  
  true_mean <- tru_alpha*tru_beta
  
  fak_alpha <- fak_alpha
  fak_beta <- fak_beta
  
  ## Empirical Bayes
  mse_eb <- numeric(n)
  ## True Bayes
  mse_fb <- numeric(n)
  ## misspecified Bayes
  mse_msb <- numeric(n)
  for (i in 5:n){
    lambda <- rgamma(i, shape = tru_alpha, rate = tru_beta^(-1))
    f_x <- rpois(i, lambda)
    mean(lambda)
    lam_f <- cbind(lambda, f_x)
    ### Empirical Bayes Posterior Parameter estimation
    emp_b_est_beta <- var(lam_f[,2])/mean(lam_f[,2])
    emp_b_est_alpha <- (mean(lam_f[,2])^2)/(var(lam_f[,2]))
    
    emp_b_est_beta;emp_b_est_alpha
    ### Empirical Bayes Lambda estimation
    emb_b_est_lambda <- emp_b_est_beta*emp_b_est_alpha
    
    ### Full Bayes Parameter estimation
    f_b_est_beta <- ((1/tru_beta) + nrow(lam_f))^(-1)
    f_b_est_alpha <- sum(lam_f[,2]) + tru_alpha
    
    ### Full Bayes Lambda estimation
    f_b_est_lambda <- f_b_est_beta*f_b_est_alpha
    
    ### Gamma prior, with scale = 1/2, rate = 1 
    mp_b_est_beta <- ((1/(fak_beta)) + nrow(lam_f)) ^(-1)
    mp_b_est_alpha <- (sum(lam_f[,2]) + fak_alpha - 1)
    
    mp_b_est_lambda <- mp_b_est_beta*mp_b_est_alpha
    
    ### Discrete Uniform prior 
    # (e^(-nrow(lam_f)) + e^(-2*nrow(lam_f)) * 2^(sum(lam_f[,2])) +)
    
    ### Calculate squared errors for each method and append to matrix
    squared_error[i, 2] <- (emb_b_est_lambda - lam_f[i ,1])^2
    squared_error[i, 1] <- (f_b_est_lambda - lam_f[i, 1])^2
    squared_error[i, 3] <- (mp_b_est_lambda - lam_f[i, 1])^2
    
    # ### Add test point to the observed data
    # lam_f_train <- rbind(lam_f_train, lam_f_test[j,])
  }
  return(squared_error)
}

# se <- pois_gam_sim(1000, alpha = 1, beta = 1, fak_alpha = 1/2, fak_beta = 1)
# plot(se[,1])
# plot(x <- 1:1000,y=cumsum(se[,1])/seq(1, 1000), type="l")
n = 100
reps <- 1e4

emp_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
h_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
m_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)

for (i in 1:reps){
  sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 1/5)
  # sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 1, fak_beta = 7)
  
  emp_bayes_mat[,i] <- sim[,2]
  h_bayes_mat[,i] <- sim[,1]
  m_bayes_mat[,i] <- sim[,3]
}

plot(rowMeans(emp_bayes_mat, na.rm = TRUE) - rowMeans(m_bayes_mat, na.rm = TRUE))

library(ggplot2)
library(dplyr)
meansdf <- as.data.frame(cbind(rowMeans(emp_bayes_mat, na.rm = TRUE), 
                 rowMeans(m_bayes_mat, na.rm = TRUE), 
                 rowMeans(h_bayes_mat, na.rm = TRUE)))

colnames(meansdf) <- c("Emp_Bayes", "M_bayes", "H_bayes")

meansdf <- meansdf %>% 
  dplyr::filter(!row_number() %in% c(1, 2, 3, 4)) %>%
  mutate(rname = row_number())

meansdf %>% tidyr::gather("id", "value", 1:3) %>%
  ggplot(., aes(rname, value, col = id)) +
  geom_point() + 
  geom_line()



