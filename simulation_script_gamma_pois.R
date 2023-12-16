set.seed(8053)
library(zoo)
library(ggplot2)
library(dplyr)

### Simulates data by drawing a lambda from a gamma distribution with parameters 
### alpha, beta (scale), then draws upup to N  x_i's from poisson distribution 
### with mean lambda. We then test the various estimates and return a matrix of
### the squared errors. 
pois_gam_sim_lambda_const<- function(n, alpha, beta, fak_alpha, fak_beta){
  
  # squared_error <- matrix(NA, ncol = 4, nrow = n)
  # squared_error <- matrix(NA, ncol = 3, nrow = n)
  squared_error <- matrix(NA, ncol = 6, nrow = n)
  
  
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
  
  for (i in 10:n){
    # i = 10
    lambda <- rgamma(1, shape = tru_alpha, rate = tru_beta^(-1))
    f_x <- rpois(i, lambda)
    lam_f <- cbind(lambda, f_x)
    
    ### Non Parametric Empirical Bayes
    poi <- as.data.frame(table(f_x), stringsAsFactors = FALSE)
    
    floor(lambda - 1.5*sqrt(lambda))
    bins <- seq(floor(lambda - sqrt(lambda)), floor(lambda+1.5*sqrt(lambda)), by=8)
    if(max(f_x) > bins[length(bins)]) {
      bins <- c(-0.1, bins, max(f_x))
    }
    
    # poi <- as.data.frame(table(f_x), stringsAsFactors = FALSE)
    # poi$f_x <- as.numeric(poi$f_x)
    # lambda_npb <- c()
    
    
    f_x <- as.data.frame(f_x)
    poi <- f_x %>%
      mutate(points_bin =
               cut(f_x,
                   breaks = bins)) %>%
      group_by(points_bin) %>% summarise(f_x = mean(f_x), Freq = n())
    
    
    poi$f_x <- as.numeric(poi$f_x)
    lambda_npb <- c()
    
    for(j in 1:nrow(poi) - 1){
      lambda_npb[j] <- (poi$f_x[j] + 1) * (poi$Freq[j+1])/poi$Freq[j]
    }
    # mean((lambda_npb - lambda)^2)
    squared_error[i, 6] <- mean((lambda_npb - lambda)^2)
    
    ### Full Bayes Parameter estimation
    f_b_est_beta <- ((1/tru_beta) + nrow(lam_f))^(-1)
    f_b_est_alpha <- sum(lam_f[,2]) + tru_alpha
    
    ### Full Bayes Lambda estimation
    f_b_est_lambda <- f_b_est_beta*f_b_est_alpha
    
    ### Gamma prior, with scale = 1/2, rate = 1 
    # mp_b_est_beta <- ((1/(fak_beta)) + nrow(lam_f)) ^(-1)
    # mp_b_est_alpha <- (sum(lam_f[,2]) + fak_alpha)
    # 
    # mp_b_est_lambda <- mp_b_est_beta*mp_b_est_alpha
    
    # emp_b_est_theta <- ##gamma_MLE(f_x)
    # emp_b_est_alpha <- ##emp_b_est_theta$alpha
    # emp_b_est_beta <- ##emp_b_est_theta$beta
    # 
    ### Empirical Bayes Lambda estimation
    emb_b_est_lambda <- ( mm_est(lam_f[,2])[1] + sum(lam_f[,2]))/
      ((1/( mm_est(lam_f[,2])[2])) + nrow(lam_f))
    
    # mm_est(f_x)[1]*mm_est(f_x)[2]
    ### Improper prior of 1
    # imp_p_est_lambda_1 <- (sum(lam_f[,2]) + 1)/nrow(lam_f)
    
    ### Improper prior of \lambda ^(-1)
    imp_p_est_lambda_2 <- sum(lam_f[,2] + 1/2)/nrow(lam_f)
    
    ### Discrete Uniform prior 
    # (e^(-nrow(lam_f)) + e^(-2*nrow(lam_f)) * 2^(sum(lam_f[,2])) +)
    
    ### Calculate squared errors for each method and append to matrix
    if(mm_est(lam_f[,2])[1] != 10000){
      squared_error[i, 2] <- (emb_b_est_lambda - lam_f[i ,1])^2
    }
    else{
      squared_error[i, 2] <- NA
    }
    squared_error[i, 1] <- (f_b_est_lambda - lam_f[i, 1])^2
    # squared_error[i, 3] <- (mp_b_est_lambda - lam_f[i, 1])^2
    # squared_error[i, 4] <- (imp_p_est_lambda_1 - lam_f[i, 1])^2
    squared_error[i, 5] <- (imp_p_est_lambda_2 - lam_f[i, 1])^2
    
    
    # ### Add test point to the observed data
    # lam_f_train <- rbind(lam_f_train, lam_f_test[j,])
  }
  squared_error[,2]
  return(squared_error)
}

#### BEGIN SIMULATION
set.seed(7)
n = 1000
reps <- 10000

### The following lines initialize matrices to store simulation results
emp_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
h_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
m_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
n_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
# i_bayes_mat_1 <- matrix(data = NA, nrow = n, ncol = reps)
i_bayes_mat_2 <- matrix(data = NA, nrow = n, ncol = reps)

# emp_bayes_mat_2 <- matrix(data = NA, nrow = n, ncol = reps)
# h_bayes_mat_2 <- matrix(data = NA, nrow = n, ncol = reps)
# m_bayes_mat_2 <- matrix(data = NA, nrow = n, ncol = reps)
# n_bayes_mat_2 <- matrix(data = NA, nrow = n, ncol = reps)
# i_bayes_mat_1_2 <- matrix(data = NA, nrow = n, ncol = reps)
# i_bayes_mat_2_2 <- matrix(data = NA, nrow = n, ncol = reps)
# 
# emp_bayes_mat_3 <- matrix(data = NA, nrow = n, ncol = reps)
# h_bayes_mat_3 <- matrix(data = NA, nrow = n, ncol = reps)
# m_bayes_mat_3 <- matrix(data = NA, nrow = n, ncol = reps)
# n_bayes_mat_3 <- matrix(data = NA, nrow = n, ncol = reps)
# i_bayes_mat_1_3 <- matrix(data = NA, nrow = n, ncol = reps)
# i_bayes_mat_2_3 <- matrix(data = NA, nrow = n, ncol = reps)
# 
# emp_bayes_mat_4 <- matrix(data = NA, nrow = n, ncol = reps)
# h_bayes_mat_4 <- matrix(data = NA, nrow = n, ncol = reps)
# m_bayes_mat_4 <- matrix(data = NA, nrow = n, ncol = reps)
# n_bayes_mat_4 <- matrix(data = NA, nrow = n, ncol = reps)
# i_bayes_mat_1_4 <- matrix(data = NA, nrow = n, ncol = reps)
# i_bayes_mat_2_4 <- matrix(data = NA, nrow = n, ncol = reps)

for (i in 1:reps){
  # sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 1/5)
  # sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 1, fak_beta = 7)
  
  #### MSE better for emp_bayes
  sim_1 <- pois_gam_sim_lambda_const(n, alpha = 20, beta = 7, fak_alpha = 2, fak_beta = 1/10)
  ### Same values Maritz looked at 
  # sim_1 <- pois_gam_sim_lambda_const(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 3)
  # sim_2 <- pois_gam_sim_lambda_const(n, alpha = 10, beta = 1/2, fak_alpha = 2, fak_beta = 3)
  # sim_3 <- pois_gam_sim_lambda_const(n, alpha = 25, beta = 1/5, fak_alpha = 2, fak_beta = 3)
  # sim_4 <- pois_gam_sim_lambda_const(n, alpha = 5, beta = 1/5, fak_alpha = 2, fak_beta = 3)
  # 
  
  ### MSE better for m_bayes
  # sim <- pois_gam_sim_lambda_const(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 1) #### 15
  
  emp_bayes_mat[,i] <- sim_1[,2]
  h_bayes_mat[,i] <- sim_1[,1]
  m_bayes_mat[,i] <- sim_1[,3]
  # i_bayes_mat_1[,i] <- sim_1[,4]
  i_bayes_mat_2[,i] <- sim_1[,5]
  n_bayes_mat[,i] <- sim_1[,6]
  
  # emp_bayes_mat_2[,i] <- sim_2[,2]
  # h_bayes_mat_2[,i] <- sim_2[,1]
  # m_bayes_mat_2[,i] <- sim_2[,3]
  # i_bayes_mat_1_2[,i] <- sim_2[,4]
  # i_bayes_mat_2_2[,i] <- sim_2[,5]
  # n_bayes_mat_2[,i] <- sim_2[,6]
  # 
  # emp_bayes_mat_3[,i] <- sim_3[,2]
  # h_bayes_mat_3[,i] <- sim_3[,1]
  # m_bayes_mat_3[,i] <- sim_3[,3]
  # i_bayes_mat_1_3[,i] <- sim_3[,4]
  # i_bayes_mat_2_3[,i] <- sim_3[,5]
  # n_bayes_mat_3[,i] <- sim_3[,6]
  # 
  # emp_bayes_mat_4[,i] <- sim_4[,2]
  # h_bayes_mat_4[,i] <- sim_4[,1]
  # m_bayes_mat_4[,i] <- sim_4[,3]
  # i_bayes_mat_1_4[,i] <- sim_4[,4]
  # i_bayes_mat_2_4[,i] <- sim_4[,5]
  # n_bayes_mat_4[,i] <- sim_4[,6]
  
}



### MM Helper function 
mm_est <- function(x){
  emp_b_est_alpha <- ifelse((mean(x^2) - mean(x) - mean(x)^2)<=0, 
                            10000, 
                            mean(x)^2/(mean(x^2) - mean(x) - mean(x)^2))
  emp_b_est_beta <- (emp_b_est_alpha/mean(x))^(-1)
  return(c(emp_b_est_alpha, emp_b_est_beta))
}

### Function to complute MLE of Gamma - Poisson distribution
gamma_MLE<-function(x, tol = .1, max.itter = 10){
  tryCatch({
    itter = 0
    n<-length(x)
    mean_x<-mean(x)

    # initiate the convergence and alpha value
    converg<-10
    alpha_last <- sqrt(mean(x)) ##mm_est(x)[1]
    beta_last <- sqrt(mean(x)) ##mm_est(x)[2]
    
    # initiate two vectors to store alpha and beta in each step
    alpha_vec <- alpha_last
    beta_vec <- beta_last
    
    # Newton-Raphson
    while(converg>tol & itter < max.itter){
      
      ## Calc first derivative
      der1_a<- log(beta_last) + sum(digamma(alpha_last + x))/n - digamma(alpha_last) - log(beta_last + 1)
      
      der1_b <- (alpha_last / beta_last ) - (mean(x) + alpha_last)/(beta_last + 1)
      
      #elements of hessian
      der2_aa<- mean(trigamma(alpha_last + x)) - trigamma(alpha_last)
      
      der2_ab <- (beta_last^2 + beta_last)^(-1)
      
      der2_bb <- -(alpha_last)/(beta_last^2) + (alpha + mean(x))/((beta_last+1)^2)
      
      ### Hessian
      jac <- cbind(c(der2_aa, der2_ab), c(der2_ab, der2_bb))
      # solve(jac)%*%as.matrix(c(der1_a, der1_b))
      
      lam_new <- c(alpha_last, beta_last ) - solve(jac)%*%as.matrix(c(der1_a, der1_b))
      
      alpha_est <- lam_new[1]
      beta_est <- lam_new[2]
      
      ### Check if converged 
      converg<-norm(c(alpha_last, beta_last) - lam_new)
      
      # store estimators in each step
      alpha_vec<-c(alpha_vec, alpha_est)
      beta_vec<-c(beta_vec, beta_est)
      alpha_last<-alpha_est
      beta_last <- beta_est
      itter = itter + 1
    }
    if (itter < max.itter){
      alpha<-alpha_est
      beta<- beta_est
      estimate_MLE<-data.frame(alpha,beta)
      return(estimate_MLE)
    }
    else{
      alpha_last <- mean(x)^2/(mean(x^2) - mean(x)^2)
      beta_last <- ((mean(x^2) - mean(x)^2)/(mean(x)^2))
      return(data.frame(alpha_last,beta_last))
    }
    
  }, error = function(e){
    alpha_last <- mean(x)^2/(mean(x^2) -  mean(x)^2)
    beta_last <- ((mean(x^2) -  mean(x)^2)/(mean(x)^2))
    return(data.frame(alpha_last,beta_last))
  })
}

### Another function to simulate data and test estimates.
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
  for (i in 10:n){
    lambda <- rgamma(1, shape = tru_alpha, rate = tru_beta^(-1))
    f_x <- rpois(i, lambda)
    lam_f <- cbind(lambda, f_x)
    ### Empirical Bayes Posterior Parameter estimation
    emp_b_est_beta <- n_raph_gamma(f_x[,2])[2]##var(lam_f[,2])/mean(lam_f[,2])
    emp_b_est_alpha <- n_raph_gamma(f_x[,2])[1]##(mean(lam_f[,2])^2)/(var(lam_f[,2]))

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


### Write files with simulation results
write.csv(emp_bayes_mat, "emp_bayes_sim_se_10000_mbet_conf.csv")
write.csv(m_bayes_mat, "ms_bayes_sim_se_10000_mbet_conf.csv")
write.csv(h_bayes_mat, "fb_bayes_sim_se_10000_mbet_conf.csv")
write.csv(i_bayes_mat_2, "ib_bayes_sim_se_10000_mbet_conf.csv")
write.csv(n_bayes_mat, "np_bayes_sim_se_10000_mbet_conf.csv")
# 
# write.csv(emp_bayes_mat_2, "emp_bayes_sim_se_10000_mbet_10_2.csv")
# write.csv(m_bayes_mat_2, "ms_bayes_sim_se_1000_mbet_10_2.csv")
# write.csv(h_bayes_mat_2, "fb_bayes_sim_se_10000_mbet_10_2.csv")
# write.csv(i_bayes_mat_2_2, "ib_bayes_sim_se_10000_mbet_10_2.csv")
# write.csv(n_bayes_mat_2, "np_bayes_sim_se_1000_mbet_10_2.csv")
# 
# write.csv(emp_bayes_mat_3, "emp_bayes_sim_se_10000_mbet_25_5.csv")
# write.csv(m_bayes_mat_3, "ms_bayes_sim_se_1000_mbet_25_5.csv")
# write.csv(h_bayes_mat_3, "fb_bayes_sim_se_10000_mbet_25_5.csv")
# write.csv(i_bayes_mat_2_3, "ib_bayes_sim_se_10000_mbet_25_5.csv")
# write.csv(n_bayes_mat_3, "np_bayes_sim_se_1000_mbet_25_5.csv")
# 
# write.csv(emp_bayes_mat_4, "emp_bayes_sim_se_10000_mbet_5_5.csv")
# write.csv(m_bayes_mat_4, "ms_bayes_sim_se_1000_mbet_5_5.csv")
# write.csv(h_bayes_mat_4, "fb_bayes_sim_se_10000_mbet_5_5.csv")
# write.csv(i_bayes_mat_2_4, "ib_bayes_sim_se_10000_mbet_5_5.csv")
# write.csv(n_bayes_mat_4, "np_bayes_sim_se_1000_mbet_5_5.csv")

### Helper function to calculate simulation based confidence intervals for MSE's 
confint <- function(n, mat){
  return(list(rowMeans(mat, na.rm = TRUE) -
           1.96 *
           apply(mat, 1, sd, na.rm=TRUE)/sqrt(n),
           rowMeans(mat, na.rm = TRUE) +
             1.96 *
             apply(mat, 1, sd, na.rm=TRUE)/sqrt(n)))
}
