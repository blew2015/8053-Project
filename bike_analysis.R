library(dplyr)
# library(rootSolve)
bikes <- read.csv("~/Documents/GitHub/8053-Project/bike+sharing+dataset/day.csv")

alpha_guess <- sqrt(mean(bikes$cnt))
beta_guess <- sqrt(mean(bikes$cnt))

bad_alpha <- sqrt(mean(bikes$cnt))
bad_beta <- 1/sqrt(mean(bikes$cnt))

set.seed(1234)
B = 100 
mat_f <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
mat_e <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
mat_npb <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
mat_msp <- matrix(data = NA, nrow = nrow(bikes) - 10, ncol = B)
for (b in 1:B){
  for (i in 10:(nrow(bikes)-1)){
    ### Get bootstrapped sample 
    # 
    # b = 1
    bstrap_ind <- sample(1:nrow(bikes), size = i, replace = TRUE)
    bstrap_set <- bikes[bstrap_ind, ]
    test_set <- bikes[-bstrap_ind,]
    
    #NPEB results
    poi <- bstrap_set %>% 
      mutate(points_bin = 
               cut(bstrap_set$cnt, 
                   breaks=c(0, 1000, 2000, 4000, 5000, 6000, 7000, 8000, max(bikes$cnt)))) %>% 
      group_by(points_bin) %>% summarise(mean = mean(cnt), freq = n())
    lambda <- c()
    for(j in 1:nrow(poi) - 1){
      lambda[j] <- (poi$mean[j] + 1) * (poi$freq[j+1])/poi$freq[j]
    }
    
    mat_npb[i-10 + 1, b] <- mean((lambda - mean(test_set$cnt))^2)
    
    ### Calculate the empirical bayes estimates
    emp_b_est_b <- (var(bstrap_set$cnt))/mean(bstrap_set$cnt)
    emp_b_est_a <- mean(bstrap_set$cnt)^2/(var(bstrap_set$cnt))
    emp_b_est_lam <- emp_b_est_a*emp_b_est_b

    ### Calculate the hierarchical bayes estimates
    f_b_est_beta <- ((1/beta_guess) + nrow(bstrap_set))^(-1)
    f_b_est_alpha <- sum(bstrap_set$cnt) + alpha_guess
    f_b_est_lam <-f_b_est_beta*f_b_est_alpha
    
    ### Calculate the hierarchical bayes estimates
    m_b_est_beta <- ((1/bad_beta) + nrow(bstrap_set))^(-1)
    m_b_est_alpha <- sum(bstrap_set$cnt) + bad_alpha
    m_b_est_lam <-m_b_est_beta*m_b_est_alpha
    
    mat_msp[i-10 + 1, b] <- mean((test_set$cnt - m_b_est_lam)^2)
    
    ### Calculate and store emp bayes MSE
    mat_e[i-10 + 1, b] <- mean((test_set$cnt - emp_b_est_lam)^2)
    
    ### Calculate and store hierarchical bayes MSE 
    mat_f[i-10 + 1, b] <- mean((test_set$cnt - f_b_est_lam)^2)
  }
}
# write.csv(mat_f, "f_b_bike_100.csv")
# write.csv(mat_e, "e_b_bike_100.csv")
# write.csv(mat_npb, "npb_bike_100.csv")
# write.csv(mat_msp, "m_b_bike_100.csv")

# bikes_eb <- read.csv("~/Documents/GitHub/8053-Project/e_b_bike_100.csv")
bikes_eb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/e_b_bike_1000.csv"))[,2:1001]
bikes_fb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/f_b_bike_1000.csv"))[,2:1001]
bikes_npb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/npb_bike_1000.csv"))[,2:1001]
bikes_mb <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/m_b_bike_1000.csv"))[,2:1001]

meansdf <- as.data.frame(cbind(rowMeans(bikes_eb, na.rm = TRUE),
                               rowMeans(bikes_fb, na.rm = TRUE),
                               rowMeans(bikes_npb, na.rm = TRUE),
                               rowMeans(bikes_mb, na.rm = TRUE)))

diff <- rowMeans(bikes_eb, na.rm = TRUE) - rowMeans(bikes_fb, na.rm = TRUE)

### Plot of difference between empirical bayes and full bayes
diff_plt <- ggplot(as.data.frame(diff) %>% mutate(index = row_number()), 
       aes(x = index + 10, y = diff)) + 
  geom_line(aes(y = diff), color = "blue") + 
  xlim(10,100) + 
  theme_light() + 
  labs(x = "Training Set Sample Size", 
       y = "MSE(E Bayes) - MSE(H Bayes)", 
       title = "MSE(Empirical Bayes) - MSE(Hierarchical Bayes)") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) ### hjust ~ 1.5

ggsave("~/Documents/GitHub/8053-Project/diff_bike.png", diff_plt)


colnames(meansdf) <- c("Emp_Bayes", "F_bayes", "NP_Bayes", "M_Bayes")


meansdf <- meansdf %>%
  mutate(rname = row_number())

### Grpah and limit x to 0 to 100
all_bike_0_731 <- meansdf %>% tidyr::gather("id", "value", 1:4) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() +
  labs(x = "Training Set Sample Size" , 
       y = "MSE", 
       title = "Real Data Example: Bike Rentals in DC, (10 < n < 731)",
       color = "Method") + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  scale_colour_brewer(palette = "Dark2", name="Method")  ### hjust ~ 1.5
ggsave("~/Documents/GitHub/8053-Project/all_bike_0_731.png", all_bike_0_731)

### Graph and limit x to 0 to 100
all_bike_0_100 <- meansdf %>% tidyr::gather("id", "value", 1:4) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() +
  xlim(10, 100) + 
  labs(x = "Training Set Sample Size" , 
       y = "MSE", 
       title = "Real Data Example: Bike Rentals in DC, (10 < n < 100)",
       color = "Method") + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  scale_colour_brewer(palette = "Dark2", name="Method") 
ggsave("~/Documents/GitHub/8053-Project/all_bike_0_100.png", all_bike_0_100)


### Grpah and limit x to 0 to 100
emp_fb_bike_0_100 <- meansdf %>% tidyr::gather("id", "value", 1:2) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id)), size = 1)+ 
  xlim(10, 100) + 
  labs(x = "Training Set Sample Size" , 
       y = "MSE", linetype="Method", color= "Method",
       title = "Real Data Example: Bike Rentals in DC, (10 < n < 100)") + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  scale_colour_brewer(palette = "Dark2", name="Method") 
ggsave("~/Documents/GitHub/8053-Project/emp_fb_bike_0_100.png", emp_fb_bike_0_100)


### Emp_Bayes and F_Bayes method at high n 
emp_fb_bike_650_731 <- meansdf %>% tidyr::gather("id", "value", 1:2) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id)), size = 1)+ 
  xlim(650, 731) + 
  ylim(3725000, 3800000) +
  labs(x = "Training Set Sample Size" , 
       y = "MSE", linetype="Method", color= "Method",
       title = "Real Data Example: Bike Rentals in DC, (650 < n < 731)") + 
  theme_light()+ 
  theme(plot.title = element_text(hjust = 0.5)) +### hjust ~ 1.5
  scale_colour_brewer(palette = "Dark2", name="Method") 
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



