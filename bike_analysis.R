library(dplyr)



bikes <- read.csv("~/Documents/GitHub/8053-Project/bike+sharing+dataset/day.csv")

alpha_guess <- 70
beta_guess <- 60
bad_alpha <- 70
bad_beta <- 1/60

mm_est_bike <- function(x){
  emp_b_est_alpha <- ifelse(mean(x)^2/(mean(x^2) - mean(x) - mean(x)^2)<=0, 
                            1e10, 
                            mean(x)^2/(mean(x^2) - mean(x) - mean(x)^2))
  emp_b_est_beta <- (mean(x^2) - mean(x) - mean(x)^2) / mean(x)
  return(c(emp_b_est_alpha, emp_b_est_beta))
}

### Poisson distribution with Gamma Prior:
gamma_MLE<-function(x, tol = 0.02, max.itter = 1000){
  tryCatch({
    itter = 0
    n<-length(x)
    mean_x<-mean(x)
    
    # initiate the convergence and alpha value
    converg<-10
    alpha_last <- alpha_guess
    beta_last <- beta_guess
    
    
    # initiate two vectors to store alpha and beta in each step
    alpha_vec<-alpha_last
    beta_vec <- beta_last
    
    # Newton-Raphson
    while(converg>tol & itter < max.itter){
      
      ## Calc first derivative
      der1_a<- log(beta_last) + mean(digamma(alpha_last + x)) - digamma(alpha_last) - log(beta_last + 1)
      
      der1_b <- (alpha_last / beta_last ) - (mean(x) + alpha_last)/(beta_last + 1)
      
      #elements of hessian
      der2_aa<- mean(trigamma(alpha_last + x)) - trigamma(alpha_last)
      
      der2_ab <- (beta_last^2 + beta_last)
      
      der2_bb <- -(alpha_last)/(beta_last^2) + (alpha + mean(x))/(beta_last+1)^2
      
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
      print(itter)
    }
    if (itter < max.itter){
      alpha<-alpha_est
      beta<- beta_est
      estimate_MLE<-data.frame(alpha,beta)
      return(estimate_MLE)
    }
    else{
      alpha <- mean(x)^2/(mean(x^2) - mean(x)^2)
      beta <- ((mean(x^2) - mean(x)^2)/(mean(x)^2))
      return(data.frame(alpha,beta))
    }
    
  }, error = function(e){
    alpha <- mean(x)^2/(mean(x^2) -  mean(x)^2)
    beta <- ((mean(x^2) -  mean(x)^2)/(mean(x)^2))
    return(data.frame(alpha_last,beta_last))
  })
  
}

##### Seed
set.seed(1234)

### Initialize B and matrices to hold standard error estimates
B = 1000 
mat_f <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
mat_e <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
mat_npb <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
mat_msp <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
mat_ip <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)

### Calculate parameter estimates and their associated squared errors
for (b in 1:B){
  for (i in 10:(nrow(bikes)-1)){
    ### Get bootstrapped sample 
    # 
    # b = 1
    # i = 10
    bstrap_ind <- sample(1:nrow(bikes), size = i, replace = TRUE)
    bstrap_set <- bikes[bstrap_ind, ]
    test_set <- bikes[-bstrap_ind,]
    
    # #NPEB results
    # poi <- bstrap_set %>% 
    #   mutate(points_bin = 
    #            cut(bstrap_set$cnt, 
    #                breaks=c(0, 1000, 2000, 4000, 5000, 6000, 7000, 8000, max(bikes$cnt)))) %>% 
    #   group_by(points_bin) %>% summarise(mean = mean(cnt), freq = n())
    # lambda <- c()
    # for(j in 1:nrow(poi) - 1){
    #   lambda[j] <- (poi$mean[j] + 1) * (poi$freq[j+1])/poi$freq[j]
    # }
    # 
    # mat_npb[i-10 + 1, b] <- mean((lambda - mean(test_set$cnt))^2)
    
    ### Calculate the empirical bayes estimates

    # mle_obj <- gamma_MLE(bstrap_set$cnt)
    # emp_b_est_beta <-  #mle_obj$beta   ###((mean(bstrap_set$cnt^2) - mean(bstrap_set$cnt) - mean(bstrap_set$cnt)^2)/(mean(bstrap_set$cnt)^2))
    # emp_b_est_alpha <- #mle_obj$alpha  ###mean(bstrap_set$cnt)^2/(mean(bstrap_set$cnt^2) - mean(bstrap_set$cnt) - mean(bstrap_set$cnt)^2)
    mm_est <- mm_est_bike(bstrap_set$cnt)
    mm_beta <- mm_est[2]
    mm_alpha <- mm_est[1]
    
    emb_b_est_beta <- (1/mm_beta + nrow(bstrap_set))^-1
    emb_b_est_alpha <- sum(bstrap_set$cnt) + mm_alpha
    emb_b_est_lam <-  emb_b_est_beta * emb_b_est_alpha
    ### Calculate the hierarchical bayes estimates
    f_b_est_beta <- ((1/beta_guess) + nrow(bstrap_set))^(-1)
    f_b_est_alpha <- sum(bstrap_set$cnt) + alpha_guess
    f_b_est_lam <-f_b_est_beta*f_b_est_alpha
    
    ### Calculate the hierarchical bayes estimates
    m_b_est_beta <- ((1/bad_beta) + nrow(bstrap_set))^(-1)
    m_b_est_alpha <- sum(bstrap_set$cnt) + bad_alpha
    m_b_est_lam <-m_b_est_beta*m_b_est_alpha
    
    ## Jeffrey's prior
    i_b_est_lam <- sum(bstrap_set$cnt + 1/2)/nrow(bstrap_set)
    mat_ip[i-10 + 1, b] <- mean((test_set$cnt - i_b_est_lam)^2)
    
    
    mat_msp[i-10 + 1, b] <- mean((test_set$cnt - m_b_est_lam)^2)
    
    ### Calculate and store emp bayes MSE
    if(mm_alpha != 1e10){
      mat_e[i-10 + 1, b] <- mean((test_set$cnt - emb_b_est_lam)^2)
    }
    else {
      mat_e[i-10 + 1, b] <- NA
    }
    ### Calculate and store hierarchical bayes MSE 
    mat_f[i-10 + 1, b] <- mean((test_set$cnt - f_b_est_lam)^2)
  }
}
# write.csv(mat_f, "f_b_bike_1000.csv")
# write.csv(mat_e, "e_b_bike_1000.csv")
# write.csv(mat_npb, "npb_bike_1000.csv")
write.csv(mat_msp, "m_b_bike_1000.csv")
# write.csv(mat_ip, "i_p_bike_1000.csv")
# 
# 
# # bikes_eb <- read.csv("~/Documents/GitHub/8053-Project/e_b_bike_100.csv")
bikes_eb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/e_b_bike_1000.csv"))[,2:1001]
bikes_fb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/f_b_bike_1000.csv"))[,2:1001]
bikes_npb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/npb_bike_1000 copy.csv"))[,2:1001]
bikes_mb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/m_b_bike_1000.csv"))[,2:1001]
bikes_ib <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/i_p_bike_1000.csv"))[,2:1001]
meansdf <- as.data.frame(cbind(rowMeans(bikes_eb, na.rm = TRUE),
                               rowMeans(bikes_fb, na.rm = TRUE),
                               rowMeans(bikes_npb, na.rm = TRUE),
                               rowMeans(bikes_mb, na.rm = TRUE),
                               rowMeans(bikes_ib, na.rm = TRUE)))

diff <- rowMeans(bikes_eb, na.rm = TRUE) - rowMeans(bikes_fb, na.rm = TRUE)
diff <- rowMeans(bikes_eb, na.rm = TRUE) - rowMeans(bikes_ib, na.rm = TRUE)
diff <- rowMeans(bikes_ib, na.rm = TRUE) - rowMeans(bikes_fb, na.rm = TRUE)

### Plot of difference between empirical bayes and full bayes
diff_plt <- ggplot(as.data.frame(diff) %>% mutate(index = row_number()),
       aes(x = index + 10, y = diff)) +
  geom_line(aes(y = diff), color = "blue") +
  xlim(10,735) +
  theme_light() +
  labs(x = "Training Set Sample Size",
       y = "MSE(E Bayes) - MSE(Strong Bayes)",
       title = "MSE(E Bayes) - MSE(Strong Bayes)") +
  theme(plot.title = element_text( size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))


ggsave("~/Documents/GitHub/8053-Project/diff_bike_ebfb.png", diff_plt)


colnames(meansdf) <- c("Param EBayes", "Strong Prior", "Nonparam EB", "Poor Prior", "Improper Prior")


meansdf <- meansdf %>%
  mutate(rname = row_number())

### Grpah and limit x to 0 to 100
all_bike_0_731 <- meansdf %>% tidyr::gather("id", "value", 1:5) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() +
  labs(x = "Training Set Sample Size" ,
       y = "MSE",
       title = "Real Data Example: Bike Rentals in DC, (10 < n < 731)",
       color = "Method") +
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  theme(plot.title = element_text( size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values = c("red", "purple", "orange", "blue", "black"),
                     name = "Method") 
ggsave("~/Documents/GitHub/8053-Project/all_bike_0_731.png", all_bike_0_731)

### Graph and limit x to 0 to 100
all_bike_0_100 <- meansdf %>% tidyr::gather("id", "value", 1:) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() +
  xlim(10, 100) +
  labs(x = "Training Set Sample Size" ,
       y = "MSE",
       title = "Real Data Example: Bike Rentals in DC, (10 < n < 100)",
       color = "Method") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  scale_color_manual(values = c("red", "purple", "orange", "blue", "black"),
                     name = "Method")
ggsave("~/Documents/GitHub/8053-Project/all_bike_0_100.png", all_bike_0_100)


### Grpah and limit x to 0 to 100
emp_fb_bike_0_100 <- meansdf %>% tidyr::gather("id", "value", c(1:2, 5)) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id)), size = 1)+
  xlim(10, 100) +
  labs(x = "Training Set Sample Size" ,
       y = "MSE", linetype="Method", color= "Method",
       title = "Real Data Example: Bike Rentals in DC") +
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  theme(plot.title = element_text( size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values = c("blue", "purple", "orange", "blue", "black"),
                     name = "Method")
ggsave("~/Documents/GitHub/8053-Project/emp_fb_bike_0_100.png", emp_fb_bike_0_100)


### Emp_Bayes and F_Bayes method at high n
emp_fb_bike_650_731 <- meansdf %>% tidyr::gather("id", "value", c(1:2, 5)) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id)), size = 1)+
  xlim(200, 731) +
  ylim(3725000, 3800000) +
  labs(x = "Training Set Sample Size" ,
       y = "MSE", linetype="Method", color= "Method",
       title = "Real Data Example: Bike Rentals in DC") +
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  theme(plot.title = element_text( size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values = c("red", "purple", "orange", "blue", "black"),
                     name = "Method")
ggsave("~/Documents/GitHub/8053-Project/emp_fb_bike_650_731.png", emp_fb_bike_650_731)


# meansdf %>% tidyr::gather("id", "value", 1) %>%
#   ggplot(., aes(rname + 10, value, col = id)) +
#   geom_line() +
#   xlim(100, 200) +
#   ylim(3750000, 3790000) +
#   labs(x = "Training Set Sample Size" ,
#        y = "MSE",
#        title = "Real Data Example: Rental Bikes in DC") +
#   theme_light()



