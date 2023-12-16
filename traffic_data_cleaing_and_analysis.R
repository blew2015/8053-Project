library(dplyr)
library(rootSolve)
accidents <- read.csv("~/Documents/GitHub/8053-Project/Accidents0515.csv")

summary(accidents)

max(accidents$Speed_limit)
count_by_month_2 <- accidents %>% mutate(
  Year = substr(Date, nchar(Date)-4+1, nchar(Date)),
  Month = substr(Date, nchar(Date)-7+1, nchar(Date) - 5),
  Day = substr(Date, 1, 2)) %>% 
  select(-Date) %>%
  filter(Year == 2014, Speed_limit < 60, Speed_limit > 45) %>% 
  group_by(Year, X1st_Road_Number, X1st_Road_Class,  X2nd_Road_Number, X2nd_Road_Class) %>%
  summarise(count = n()) %>% 
  filter(Year == "2014") 
hist(count_by_month_2$count)
write.csv(count_by_month_2, "~/Documents/GitHub/8053-Project/Accidents_filtered_01_2014_2.csv")

set.seed(1234)

tru_beta <- 1
tru_alpha <- 1

# for (i in 10:nrow(count_by_month_2)-1) {
#   real <- perm[1:i,] 
#   emp_b_est_b <- var(real$count)/mean(real$count)
#   emp_b_est_a <- mean(real$count)^2/var(real$count)
#   emp_b_est_lam <- emp_b_est_a*emp_b_est_b
#   
#   f_b_est_beta <- ((1/tru_beta) + nrow(real))^(-1)
#   f_b_est_alpha <- sum(real$count) + tru_alpha
#   f_b_est_lam <-f_b_est_beta*f_b_est_alpha
#   
#   eb_se_vec[i-9] <- (perm[i+1,]$count - emp_b_est_lam)^2
#   f_se_vec[i-9] <- (perm[i+1,]$count - f_b_est_lam)^2
# }
n_raph_gamma <- function(x){
  x <- bstrap_set$count
  converge <- 10
  alpha_est <- mean(x) + (1/length(x))
  beta_est <- mean(x)/alpha_est
  i <- 0
  while (converge < 0.01){
    d_1 <- log(alpha_est/mean(x)) - 
      digamma(alpha_est) +  
      mean(digamma(x + alpha_est)) - 
      log(alpha_est/mean(x) + 1)
    
    d_2 <- length(x)/alpha_est - length(x)*trigamma(alpha_est)
    
    alpha_next_itter <- alpha_est - d_1/d_2
    
    converge<-abs(alpha_next_itter-alpha_est)
    
    alpha_est <- alpha_next_itter
  }
  # alpha<-alpha_next_itter
  beta<-mean(x)/alpha_est
  return(c(alpha, beta ))
}
B = 1000
mat_f <- matrix(data = NA, nrow = nrow(count_by_month_2) - 10, ncol = B)
mat_e <- matrix(data = NA, nrow = nrow(count_by_month_2) - 10, ncol = B)
### B = 100 
for (b in 1:B){
  for (i in 10:nrow(count_by_month_2)){
    ### Get bootstrapped sample 
    # i = 10
    bstrap_ind <- sample(1:nrow(count_by_month_2), size = i, replace = TRUE)
    bstrap_set <- count_by_month_2[bstrap_ind, ]
    # bstrap_set$count
    test_set <- count_by_month_2[-bstrap_ind,]
    
    ### Calculate the empirical bayes estimates
    emp_b_est_b <- (var(bstrap_set$count) + 1/i)/mean(bstrap_set$count)
    emp_b_est_a <- mean(bstrap_set$count)^2/(var(bstrap_set$count) + 1/i)
    emp_b_est_lam <- emp_b_est_a*emp_b_est_b
    emp_b_est_lam
   
    
    ### Calculate the hierarchical bayes estimates
    f_b_est_beta <- ((1/tru_beta) + nrow(bstrap_set))^(-1)
    f_b_est_alpha <- sum(bstrap_set$count) + tru_alpha
    f_b_est_lam <-f_b_est_beta*f_b_est_alpha
    
    ### Calculate and store emp bayes MSE
    mat_e[i-10, b] <- mean((test_set$count - emp_b_est_lam)^2)
    
    ### Calculate and store hierarchical bayes MSE 
    mat_f[i-10, b] <- mean((test_set$count - f_b_est_lam)^2)
  }
}
write.csv(mat_f, "f_b_boot_1000.csv")
write.csv(mat_e, "e_b_boot_1000.csv")



# c = c(1, 1, 1, 1)
# ### https://canvas.harvard.edu/files/350816/download?download_frd=1&verifier=Jlwvto2cvU9DHcN6DKYKWJ1rXeCVZrQNTdhGJMKP
# ss <- multiroot(f = model, start = c(1, 1))
# 
# 
# par(mfrow = c(2, 1))
# plot(rowMeans(mat_e, na.rm = TRUE))
# plot(rowMeans(mat_f, na.rm = TRUE))
# 
# plot(rowMeans(mat_f, na.rm = TRUE) - rowMeans(mat_e, na.rm = TRUE))
# hist(count_by_month_2$count)
