library(dplyr)
accidents <- read.csv("~/Documents/GitHub/8053-Project/Accidents0515.csv")

summary(accidents)

max(accidents$Speed_limit)
count_by_month_2 <- accidents %>% mutate(
  Year = substr(Date, nchar(Date)-4+1, nchar(Date)),
  Month = substr(Date, nchar(Date)-7+1, nchar(Date) - 5),
  Day = substr(Date, 1, 2)) %>% 
  select(-Date) %>%
  filter(Year == 2014, Speed_limit < 60, Speed_limit > 45) %>% 
  group_by(Year, Month, X1st_Road_Number, X1st_Road_Class,  X2nd_Road_Number, X2nd_Road_Class) %>%
  summarise(count = n()) %>% 
  filter(Year == "2014", Month == "01") 

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
B = 1000
mat_f <- matrix(data = NA, nrow = nrow(count_by_month_2) - 10, ncol = B)
mat_e <- matrix(data = NA, nrow = nrow(count_by_month_2) - 10, ncol = B)
### B = 100 
for (b in 1:B){
  for (i in 10:nrow(count_by_month_2)){
    ### Get bootstrapped sample 
    bstrap_ind <- sample(1:nrow(count_by_month_2), size = i, replace = TRUE)
    bstrap_set <- count_by_month_2[bstrap_ind, ]
    test_set <- count_by_month_2[-bstrap_ind,]
    
    ### Calculate the empirical bayes estimates
    emp_b_est_b <- var(bstrap_set$count)/mean(bstrap_set$count)
    emp_b_est_a <- mean(bstrap_set$count)^2/var(bstrap_set$count)
    emp_b_est_lam <- emp_b_est_a*emp_b_est_b
    
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

plot(rowMeans(mat_e, na.rm = TRUE))
plot(rowMeans(mat_f, na.rm = TRUE))

plot(rowMeans(mat_f, na.rm = TRUE) - rowMeans(mat_e, na.rm = TRUE))
hist(count_by_month_2$count)
