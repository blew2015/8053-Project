set.seed(8053)
library(zoo)
library(ggplot2)
library(dplyr)

### Poisson distribution with Gamma Prior:

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

pois_gam_sim_lambda_const<- function(n, alpha, beta, fak_alpha, fak_beta){
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
    lambda <- rgamma(1, shape = tru_alpha, rate = tru_beta^(-1))
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

#### BEGIN SIMULATION
set.seed(7)
n = 10000
reps <- 100

emp_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
h_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
m_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)

for (i in 1:reps){
  # sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 1/5)
  # sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 1, fak_beta = 7)
  
  sim <- pois_gam_sim_lambda_const(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 1/5)
  # sim <- pois_gam_sim_lambda_const(n, alpha = 2, beta = 1/2, fak_alpha = 1, fak_beta = 15)

  emp_bayes_mat[,i] <- sim[,2]
  h_bayes_mat[,i] <- sim[,1]
  m_bayes_mat[,i] <- sim[,3]
}

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

ggplot(meansdf, aes(rname, Emp_Bayes)) +
  geom_line()


