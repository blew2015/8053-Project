set.seed(8053)
library(zoo)
library(ggplot2)
library(dplyr)

### Poisson distribution with Gamma Prior:

### Newton Raphson estimation
n_raph_gamma <- function(x){
  # x <- c(0, 3, 0, 1 ,0)
  converge <- 10
  
  alpha_est <- mean(x) + (1/length(x))
  beta_est<-mean(x)/alpha_est
  alpha_next_itter <- 0 
  
  alpha_vec <- numeric(2)
  while (converge < 0.01){
    d_1 <- log(alpha_est/mean(x)) - 
      digamma(alpha_est) +  
      mean(digamma(x + alpha_est)) - 
      log(alpha_est/mean(x) + 1)
    
    d_2 <- length(x)/alpha_est - length(x)*trigamma(alpha_est)
    
    alpha_vec[2] <- alpha_vec[1]
    alpha_vec[1] <- alpha_est - d_1/d_2
    
    converge<<-abs(alpha_next_itter-alpha_est)
    
    alpha_est <<- alpha_next_itter
  }
  alpha<-alpha_next_itter
  beta<-mean(x)/alpha_next_itter
  return(c(alpha, beta ))
}

d_1 <- function(alpha_est){
  x =  lam_f[,2]
  log(alpha_est/mean(x)) - 
    digamma(alpha_est) +  
    mean(digamma(x + alpha_est)) - 
    log(alpha_est/mean(x) + 1)
}

d_1 <-function(alpha_est){
  x =  lam_f[,2]
  log(alpha_est/mean(x)) - 
    digamma(alpha_est) +  
    mean(digamma(x + alpha_est)) - 
    log(alpha_est/mean(x) + 1)
}


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
    lambda <- rgamma(i, shape = tru_alpha, rate = tru_beta^(-1))
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

pois_gam_sim_lambda_const<- function(n, alpha, beta, fak_alpha, fak_beta){

  # squared_error <- matrix(NA, ncol = 4, nrow = n)
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
    # i = 10
    lambda <- rgamma(1, shape = tru_alpha, rate = tru_beta^(-1))
    f_x <- rpois(i, lambda)
    lam_f <- cbind(lambda, f_x)
   
    # ### Non Parametric Empirical Bayes
    # poi <- as.data.frame(table(f_x), stringsAsFactors = FALSE)
    # 
    # floor(lambda - 1.5*sqrt(lambda)) 
    # bins <- seq(floor(lambda - sqrt(lambda)), floor(lambda+1.5*sqrt(lambda)), by=1)
    # if(max(f_x) > bins[length(bins)]) {
    #   bins <- c(-0.1, bins, max(f_x))
    # }
    # 
    # # poi <- as.data.frame(table(f_x), stringsAsFactors = FALSE)
    # # poi$f_x <- as.numeric(poi$f_x)
    # # lambda_npb <- c()
    # 
    # 
    # f_x <- as.data.frame(f_x)
    # poi <- f_x %>% 
    #   mutate(points_bin = 
    #            cut(f_x, 
    #                breaks = bins)) %>%
    #   group_by(points_bin) %>% summarise(f_x = mean(f_x), Freq = n())
    # 
    # 
    # poi$f_x <- as.numeric(poi$f_x)
    # lambda_npb <- c()
    # 
    # for(j in 1:nrow(poi) - 1){
    #   lambda_npb[j] <- (poi$f_x[j] + 1) * (poi$Freq[j+1])/poi$Freq[j]
    # }
    # mean((lambda_npb - lambda)^2)
    # squared_error[i, 4] <- mean((lambda_npb - lambda)^2)
# 
    ### Empirical Bayes Posterior Parameter estimation
    # e_b_params <- n_raph_gamma(lam_f[,2])
    emp_b_est_alpha <- (mean(lam_f[,2])^2)/(var(lam_f[,2]))
    emp_b_est_beta <- var(lam_f[,2])/mean(lam_f[,2])
    
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
n = 1000
reps <- 10000

emp_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
h_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
m_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)
n_bayes_mat <- matrix(data = NA, nrow = n, ncol = reps)

for (i in 1:reps){
  # sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 1/5)
  # sim <- pois_gam_sim(n, alpha = 2, beta = 1/2, fak_alpha = 1, fak_beta = 7)
  
  #### MSE better for emp_bayes
  # sim <- pois_gam_sim_lambda_const(n, alpha = 5, beta = 1, fak_alpha = 2, fak_beta = 1/5)
  
  ### MSE better for m_bayes
  sim <- pois_gam_sim_lambda_const(n, alpha = 2, beta = 1/2, fak_alpha = 2, fak_beta = 1) #### 15

  emp_bayes_mat[,i] <- sim[,2]
  h_bayes_mat[,i] <- sim[,1]
  m_bayes_mat[,i] <- sim[,3]
  # n_bayes_mat[,i] <- sim[,4]
}

write.csv(emp_bayes_mat, "emp_bayes_sim_se_1000_mbet.csv")
write.csv(m_bayes_mat, "ms_bayes_sim_se_1000_mbet.csv")
write.csv(h_bayes_mat, "fb_bayes_sim_se_1000_mbet.csv")



emp_bayes_mat[-c(1:9),]
apply(emp_bayes_mat, 1, sd, na.rm=TRUE)
rowMeans(emp_bayes_mat, na.rm = TRUE) - 1.96 * apply(emp_bayes_mat, 1, sd, na.rm=TRUE)/sqrt(100)

confint <- function(n, mat){
  return(list(rowMeans(mat, na.rm = TRUE) - 
           1.96 * 
           apply(mat, 1, sd, na.rm=TRUE)/sqrt(n), 
           rowMeans(mat, na.rm = TRUE) +
             1.96 * 
             apply(mat, 1, sd, na.rm=TRUE)/sqrt(n)))
}



# p<-ggplot(data=data, aes(x=interval, y=OR, colour=Drug)) + 
#   geom_point() + geom_line()
# p+geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)

emp_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/emp_bayes_sim_se_1000.csv"))[,2:1001]
h_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/fb_bayes_sim_se_1000.csv"))[,2:1001]
n_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/n_bayes_sim_se_1000.csv"))[,2:1001]
m_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/ms_bayes_sim_se_1000.csv"))[,2:1001]

e_b_int <- confint(1000, emp_bayes_mat)
h_b_int <- confint(1000, h_bayes_mat)
m_b_int <- confint(1000, h_bayes_mat)

meansdf <- as.data.frame(cbind(rowMeans(emp_bayes_mat, na.rm = TRUE), 
                 rowMeans(m_bayes_mat, na.rm = TRUE), 
                 rowMeans(h_bayes_mat, na.rm = TRUE), 
                 rowMeans(n_bayes_mat, na.rm = TRUE)))

colnames(meansdf) <- c("Emp_Bayes", "M_Bayes", "H_Bayes", "NP_Bayes")

meansdf$e_b_conf_u <- e_b_int[[2]]
meansdf$e_b_conf_l <- e_b_int[[1]]
meansdf$h_b_conf_u <- h_b_int[[2]]
meansdf$h_b_conf_l <- h_b_int[[1]]
meansdf$m_b_conf_u <- m_b_int[[2]]
meansdf$m_b_conf_l <- m_b_int[[1]]

meansdf <- meansdf %>% 
  dplyr::filter(!row_number() %in% c(1, 2, 3, 4, 5, 6,7,8,9,10)) %>%
  mutate(rname = row_number())

### Full Plot with everything
all_10_1000 <- meansdf %>% tidyr::gather("id", "value", 1:4) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() + 
  ggtitle("Bayes Methods Simulation Comparison")+
  labs( x = "Training Set Sample Size", y = "MSE", color = "Method") + 
  xlim(10, 1000) +
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5)) ### hjust ~ 1.5
ggsave("~/Documents/GitHub/8053-Project/all_10_1000.png", all_10_1000)

### Grpah and limit x to 0 to 100, include confidence intervals, just ebayes + f_bayes
eb_v_hb <- meansdf %>% tidyr::gather("id", "value", c(1,3)) %>%
  mutate(upper = ifelse(
    id =="Emp_Bayes", e_b_conf_u, ifelse( 
      id == "H_Bayes", h_b_conf_u, ifelse(
        id == "M_Bayes", m_b_conf_u, 0
      ))),
    lower = ifelse(
      id =="Emp_Bayes", e_b_conf_l, ifelse( 
        id == "H_Bayes", h_b_conf_l, ifelse(
          id == "M_Bayes", m_b_conf_l, 0
        )))) %>% 
  select(c("rname", "value", "id", "upper", "lower")) %>% 
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id)), size = 1)+ 
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1) +
  xlim(10, 100) + 
  theme_light() +
  labs(title = "Empirical Bayes vs Hierarchical Bayes (10 < n < 100)", 
       x = "Training Set Sample Size", y = "MSE", color = "Method", 
       linetype = "Method") + 
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  scale_colour_brewer(palette = "Dark2", name="Method") 
ggsave("~/Documents/GitHub/8053-Project/eb_v_hb_10_100.png", eb_v_hb)

### Graph and lit x to 900 to 1000
eb_v_hb_100_1000 <- meansdf %>% tidyr::gather("id", "value", c(1,3)) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id))) + 
  xlim(100, 1000) + 
  ylim(0, 0.075) +
  theme_light() +   
  labs(title = "Empirical Bayes vs Hierarchical Bayes (100 < n < 1000)", 
       x = "Training Set Sample Size", y = "MSE", color = "Method", 
       linetype = "Method") + 
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  scale_colour_brewer(palette = "Dark2", name="Method") ### hjust ~ 1.5
ggsave("~/Documents/GitHub/8053-Project/eb_v_hb_100_1000.png", eb_v_hb_100_1000)

### FUll Plot with everything
all_0_100 <- meansdf %>% tidyr::gather("id", "value", 1:4) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() + 
  labs(title = "Bayes Methods Simulation Comparison",
        x = "Training Set Sample Size", y = "MSE", color = "Method") + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5)) +  ### hjust ~ 1.5
  xlim(10, 100)
ggsave("~/Documents/GitHub/8053-Project/all_10_100.png", all_0_100)


### Full Plot without nonparametric 
all_no_np_10_100 <- meansdf %>% tidyr::gather("id", "value", 1:3) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id)), size = 1) + 
  ggtitle("Bayes Methods Simulation Comparrison (10 < n < 100)")+
  labs( x = "Training Set Sample Size", y = "MSE", color = "Method", 
        linetype = "Method") + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5)) +  ### hjust ~ 1.5
  xlim(10, 100) +
  scale_colour_brewer(palette = "Dark2", name="Method") ### hjust ~ 1.5
ggsave("~/Documents/GitHub/8053-Project/all_no_np_10_100.png", all_no_np_10_100)

  ### Full Plot without nonparametric 
all_no_np_100_1000<- meansdf %>% tidyr::gather("id", "value", 1:3) %>%
    ggplot(., aes(rname + 10, value, col = id)) +
    geom_line(aes(linetype = as.factor(id)), size = 1) + 
    ggtitle("Bayes Methods Simulation Comparrison (100 < n < 1000)")+
    labs( x = "Training Set Sample Size", y = "MSE", color = "Method", 
          linetype = "Method") + 
    xlim(100, 1000) + 
    ylim(0, 0.5) +
    theme_light() + 
    theme(plot.title = element_text(hjust = 0))   ### hjust ~ 1.5
ggsave("~/Documents/GitHub/8053-Project/all_no_np_100_1000.png", all_no_np_100_1000)



