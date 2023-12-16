### Script to calculate values in the tables where we moved the hyperparameters 
### of the gamma distribution conjugate around
emp_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/emp_bayes_sim_se_1000_mbet_10_2.csv"))[,2:1001]
m_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ms_bayes_sim_se_1000_mbet_10_2.csv"))[,2:1001]
h_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/fb_bayes_sim_se_1000_mbet_10_2.csv"))[,2:1001]
i_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ib_bayes_sim_se_1000_mbet_10_2.csv"))[,2:1001]

### Squared Error estimates: alpha 10,  beta 1/2
round(rowMeans(emp_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
round(rowMeans(i_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
round(rowMeans(h_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)


emp_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/emp_bayes_sim_se_1000_mbet_25_5.csv"))[,2:1001]
m_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ms_bayes_sim_se_1000_mbet_25_5.csv"))[,2:1001]
h_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/fb_bayes_sim_se_1000_mbet_25_5.csv"))[,2:1001]
i_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ib_bayes_sim_se_1000_mbet_25_5.csv"))[,2:1001]

### Squared Error estimates: alpha 25,  beta 1/5
round(rowMeans(emp_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
round(rowMeans(i_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
round(rowMeans(h_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)


emp_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/emp_bayes_sim_se_1000_mbet_5_5.csv"))[,2:1001]
m_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ms_bayes_sim_se_1000_mbet_5_5.csv"))[,2:1001]
h_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/fb_bayes_sim_se_1000_mbet_5_5.csv"))[,2:1001]
i_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ib_bayes_sim_se_1000_mbet_5_5.csv"))[,2:1001]

### Squared Error estimates: alpha 5,  beta 1/5
round(rowMeans(emp_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
round(rowMeans(i_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
round(rowMeans(h_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)

emp_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/emp_bayes_sim_se_1000_mbet_2_2.csv"))[,2:1001]
m_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ms_bayes_sim_se_1000_mbet_2_2.csv"))[,2:1001]
h_bayes_mat <-as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/fb_bayes_sim_se_1000_mbet_2_2.csv"))[,2:1001]
i_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/ib_bayes_sim_se_1000_mbet_2_2.csv"))[,2:1001]

### Squared Error estimates: alpha 2,  beta 1/2
rowMeans(emp_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)]
rowMeans(i_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)]
rowMeans(h_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)]

### Squared Error estimats for non parametric EB methods
n_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/np_bayes_sim_se_1000_mbet_2_2.csv"))[,2:1001]
round(rowMeans(n_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
n_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/np_bayes_sim_se_1000_mbet_5_5.csv"))[,2:1001]
round(rowMeans(n_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
n_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/np_bayes_sim_se_1000_mbet_10_2.csv"))[,2:1001]
round(rowMeans(n_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)
n_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/moving_params_res/np_bayes_sim_se_1000_mbet_25_5.csv"))[,2:1001]
round(rowMeans(n_bayes_mat[-c(1:9),], na.rm = TRUE)[c(10 - 9, 50 - 9, 100 - 9, 200 - 9, 1000 - 9)], 5)

