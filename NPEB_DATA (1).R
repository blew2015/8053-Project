## Function to obtain MSE estimates - Parametric EB, Hierarchical bayes
case_1_ext <- function(n,p,sigma,mu,tau) {
  X <- c()
  bayes_a_post_mean <- c()
  bayes_b_post_mean <- c()
  EB_post_mean <- c()
  theta <- c()
  for (k in 1:p) {
    norm0 <- norm_gen_list(n, sigma, mu, tau) # simulate n Xis
    X <- cbind(X,norm0[[2]]) # store in X
    theta[k] <- norm0[[1]] # store theta
    bayes_a_post_mean[k] <- (n*tau^2/(n*tau^2+sigma^2))*mean(norm0[[2]]) + (sigma^2/(n*tau^2+sigma^2))*mu
    bayes_b_post_mean[k] <- mean(norm0[[2]])
  }
  s2 <- sum((colMeans(X)-mean(as.matrix(X)))^2)
  num <- (p-3)*(sigma^2)
  ratio <- num/s2 
  EB_post_mean <- (1-ratio)*colMeans(X)+(ratio)*mean(as.matrix(X))
  MSE_bayes_a <- mean((theta-bayes_a_post_mean)^2)
  MSE_bayes_b <- mean((theta-bayes_b_post_mean)^2)
  MSE_EB <- mean((theta-EB_post_mean)^2)
  MSE <- cbind(MSE_bayes_a, MSE_bayes_b, MSE_EB)
  return(MSE)
}

library(ks)

## DATA EXAMPLE - Student grades (NORMAL)
data_grades <- read.csv("/Users/hjkoo/Desktop/8053 Final Project/student_marks.csv")

data_grades <- data_grades[,-c(1)]
B = 1000
N_sample=nrow(data_grades)
mat_bayes <- matrix(data = NA, nrow = N_sample - 10, ncol = B)
mat_imp <- matrix(data = NA, nrow = N_sample - 10, ncol = B)
mat_eb <- matrix(data = NA, nrow = N_sample - 10, ncol = B)
mat_np <- matrix(data = NA, nrow = N_sample - 10, ncol = B)
mat_bayes_inc <- matrix(data = NA, nrow = N_sample - 10, ncol = B)
set.seed(1234)
### B = 100 
for (b in 1:B){
  for (i in 10:nrow(data_grades)){
    ### Get bootstrapped sample 
    bstrap_ind <- sample(1:nrow(data_grades), size = i, replace = TRUE)
    bstrap_set <- data_grades[bstrap_ind, ]
    # bstrap_set$count
    test_set <- data_grades[-bstrap_ind,]
    
    np_est <- c()
    n = nrow(bstrap_set)
    p = ncol(bstrap_set)
    sigma = sd(as.matrix(bstrap_set))
    ## Nonparametric EB
    #for (j in 1:ncol(bstrap_set)) {
     # scores_ordered <- sort(bstrap_set[,j])
      #fx <- kde(scores_ordered,h=2, eval.points = scores_ordered)
      #fx_d <- kdde(scores_ordered,h=2, deriv.order = 1, eval.points = scores_ordered)
      
      #score <- fx_d$estimate/fx$estimate
      
      ## Tweedie's formula (E(\theta|X) = x+sigma^2*score(x))
      #np_est[j] <- mean(fx$eval.points+sigma^2*score)
      
      #theta_hat_raw <- fx$eval.points+sigma^2*score
    #}
    
    ## Nonparametric EB
    scores_ordered <- sort(unlist(bstrap_set))
    fx <- kde(scores_ordered,h=2, eval.points = colMeans(bstrap_set))
    fx_d <- kdde(scores_ordered,h=2, deriv.order = 1, eval.points = colMeans(bstrap_set))
    score <- fx_d$estimate/fx$estimate
    theta_hat_raw <- fx$eval.points+sigma^2*score
    
    
    ## Empirical Bayes (Parametric)
    s2 <- sum((colMeans(bstrap_set)-mean(as.matrix(bstrap_set)))^2)*n
    num <- (p-3)*(sigma^2)
    ratio <- num/s2 
    EB_post_mean <- ((1-ratio)*colMeans(bstrap_set))+((ratio)*mean(as.matrix(bstrap_set)))
    
    bayes_a_post_mean <- c()
    bayes_b_post_mean <- c()
    bayes_i_post_mean <- c()
    
    ## Hierarchical Bayes (Normal Prior - N(70, tau = 3)), incorrect prior
    tau = 3
    mu = 70
    mu_inc = 30
    n = nrow(bstrap_set)
    for (k in 1:p) {
      bayes_a_post_mean[k] <- (n*tau^2/(n*tau^2+sigma^2))*mean(bstrap_set[,k]) + (sigma^2/(n*tau^2+sigma^2))*mu
      bayes_i_post_mean[k] <- (n*tau^2/(n*tau^2+sigma^2))*mean(bstrap_set[,k]) + (sigma^2/(n*tau^2+sigma^2))*mu_inc
      bayes_b_post_mean[k] <- mean(bstrap_set[,k])
    }
    
    theta <- colMeans(test_set)
    mat_bayes[i-10, b] <- mean((theta - bayes_a_post_mean)^2)
    mat_imp[i-10, b] <- mean((theta - bayes_b_post_mean)^2)
    mat_eb[i-10, b] <- mean((theta - EB_post_mean)^2)
    mat_np[i-10, b] <- mean((theta - theta_hat_raw)^2)
    mat_bayes_inc[i-10, b] <- mean((theta-bayes_i_post_mean)^2)
  }
}

spar(mfrow = c(2, 2))
plot(rowMeans(mat_bayes, na.rm = TRUE))
plot(rowMeans(mat_imp, na.rm = TRUE))
plot(rowMeans(mat_eb, na.rm = TRUE))
plot(rowMeans(mat_np, na.rm = TRUE))

mat_bayes_mean <- rowMeans(mat_bayes, na.rm = TRUE)
mat_imp_mean <- rowMeans(mat_imp, na.rm = TRUE)
mat_eb_mean <- rowMeans(mat_eb, na.rm = TRUE)
mat_np_mean <- rowMeans(mat_np, na.rm = TRUE)
mat_bayes_inc_mean <- rowMeans(mat_bayes_inc, na.rm = TRUE)
x <- seq(nrow(data_grades)-10)

plot_data <- data.frame(cbind(x, mat_bayes_mean, mat_imp_mean, mat_eb_mean, mat_np_mean, mat_bayes_inc_mean))

## We have B=1000 samples
library(ggplot2)
png(file="/Users/hjkoo/Desktop/8053 Final Project/Plots/NN_Data.png",
    width=600, height=350)
ggplot(plot_data, aes(x=x)) + geom_line(aes(y=mat_bayes_mean,color="Group 1", linetype="Group 1"),size=1.0) + geom_line(aes(y=mat_imp_mean,color="Group 2", linetype="Group 2"),size=1.0) + geom_line(aes(y=mat_eb_mean,color="Group 3", linetype="Group 3"),size=1.0) + geom_line(aes(y=mat_np_mean,color="Group 4", linetype="Group 4"),size=1.0) + geom_line(aes(y=mat_bayes_inc_mean,color="Group 5", linetype="Group 5"),size=1.0) +  labs(title = "Student Grades data",x = "Number of training samples",y = "Average MSE") +
  scale_color_manual(values = c("red", "blue", "green4", "purple","gray18"),
                     name = "Groups",
                     labels = c("Normal", "Improper", "Parametric EB", "Nonparametric EB","Incorrect mean Normal")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash","longdash"),
                        name = "Groups",
                        labels = c("Normal", "Improper", "Parametric EB", "Nonparametric EB", "Incorrect mean Normal"))
dev.off()
png(file="/Users/hjkoo/Desktop/8053 Final Project/Plots/NN_Data_wo_inc.png",
    width=600, height=350)
ggplot(plot_data, aes(x=x)) + geom_line(aes(y=mat_bayes_mean,color="Group 1", linetype="Group 1"),size=1) + geom_line(aes(y=mat_imp_mean,color="Group 2", linetype="Group 2"),size=1) + geom_line(aes(y=mat_eb_mean,color="Group 3", linetype="Group 3"),size=1) +labs(title = "Student Grades data", x = "Number of training samples", y = "Average MSE") +
  scale_color_manual(values = c("red", "blue", "green4"),
                     name = "Groups",
                     labels = c("Normal", "Improper", "Parametric EB")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                        name = "Groups",
                        labels = c("Normal", "Improper", "Parametric EB"))
dev.off()

png(file="/Users/hjkoo/Desktop/8053 Final Project/Plots/NN_Data.png",
    width=600, height=350)
ggplot(test_mu_plot, aes(x=test_mu_1)) + geom_line(aes(y=MSE_bayes_a_total_mu_mean,color="Group 1", linetype="Group 1"), alpha=0.3) + geom_line(aes(y=MSE_bayes_b_total_mu_mean,color="Group 2", linetype="Group 2")) + geom_line(aes(y=MSE_EB_total_mu_mean,color="Group 3", linetype="Group 3")) + geom_line(aes(y=MSE_NP_total_mu_mean,color="Group 4", linetype="Group 4")) +  labs(title = "Varying values of mu",
                                                                                                                                                                                                                                                                                                                                                                                        x = "Mu",
                                                                                                                                                                                                                                                                                                                                                                                        y = "Average MSE") +
  scale_color_manual(values = c("red", "blue", "green", "purple"),
                     name = "Groups",
                     labels = c("Normal", "Improper", "Parametric EB", "Nonparametric EB")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"),
                        name = "Groups",
                        labels = c("Normal", "Improper", "Parametric EB", "Nonparametric EB"))
dev.off()




## DATA EXAMPLE - Traffic data (POISSON)
data_traffic <- read.csv("/Users/hjkoo/Desktop/STAT 8053/Project/Accidents_filtered_01_2014_2.csv")

hist(data_traffic$count) # Looks one inflated
table(data_traffic$count)

## Split training and testing 

# f(x) can be estimated using y_x/N and f(x+1) can be estimated using y_{x+1}/N

poi <- as.data.frame(table(X_poi))
poi$X_poi <- as.numeric(poi$X_poi) - 1

lambda <- c()
for (i in 1:nrow(poi)-1){
  lambda[i] <- (poi$X_poi[i]+1) * (poi$Freq[i+1]) / (poi$Freq[i])
}

lambda
plot(lambda)

## MSE
MSE_NPEB_poi <- mean((lambda-1)^2)


## DATA EXAMPLE - Length of Stay at hospital (NORMAL)



## DATA EXAMPLE - Bike sharing (POISSON)
