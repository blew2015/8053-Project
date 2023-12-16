### Read in data
emp_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/emp_bayes_sim_se_1000_mbet.csv"))[,2:1001]
h_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/fb_bayes_sim_se_1000_mbet.csv"))[,2:1001]
n_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/np_bayes_sim_se_1000_mbet.csv"))[,2:1001]
m_bayes_mat <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/ms_bayes_sim_se_1000_mbet.csv"))[,2:1001]
i_bayes_mat_2 <- as.matrix(read.csv("~/Documents/GitHub/8053-Project/Simulation Results/ib_bayes_sim_se_1000_mbet.csv"))[,2:1001]

### Generate confidence intervals
e_b_int <- confint(10000, emp_bayes_mat)
h_b_int <- confint(10000, h_bayes_mat)
m_b_int <- confint(10000, m_bayes_mat)
i_b_int <- confint(10000, i_bayes_mat_1)

### Create dataframe
meansdf <- as.data.frame(cbind(rowMeans(m_bayes_mat, na.rm = TRUE),
                               rowMeans(emp_bayes_mat, na.rm = TRUE),
                               rowMeans(h_bayes_mat, na.rm = TRUE),
                               rowMeans(i_bayes_mat_2, na.rm = TRUE),
                               rowMeans(n_bayes_mat, na.rm = TRUE)))

colnames(meansdf) <- c("Misspecified Bayes", "P EBayes", "Gamma Pois", "Imp Prior", "NP EBayes")

meansdf$e_b_conf_u <- e_b_int[[2]]
meansdf$e_b_conf_l <- e_b_int[[1]]
meansdf$h_b_conf_u <- h_b_int[[2]]
meansdf$h_b_conf_l <- h_b_int[[1]]
meansdf$m_b_conf_u <- m_b_int[[2]]
meansdf$m_b_conf_l <- m_b_int[[1]]
meansdf$i_b_conf_u <- i_b_int[[2]]
meansdf$i_b_conf_l <- i_b_int[[1]]

meansdf <- meansdf %>%
  dplyr::filter(!row_number() %in% c(1, 2, 3, 4, 5, 6,7,8,9,10)) %>%
  mutate(rname = row_number())

### Generate Plots

### Full Plot with everything
all_10_1000 <- meansdf %>% tidyr::gather("id", "value", c(1:5)) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() +
  ggtitle("Bayes Methods Simulation Comparison")+
  labs( x = "Training Set Sample Size", y = "MSE", color = "Method") +
  xlim(0, 1000) +
  # theme_light() +
  # ylim(0, 100)+
  
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values = c("red", "purple", "orange", "blue", "black"),
                     name = "Method") 
ggsave("~/Documents/GitHub/8053-Project/all_10_1000.png", all_10_1000)

no_npb_10_1000 <- meansdf %>% tidyr::gather("id", "value", c(2:4)) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line() +
  ggtitle("Bayes Methods Simulation Comparison")+
  labs( x = "Training Set Sample Size", y = "MSE", color = "Method") +
  xlim(0, 100) +
  # theme_light() +
  # ylim(0, 100)+
  
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values = c("red", "purple", "orange", "blue", "black"),
                     name = "Method") 
ggsave("~/Documents/GitHub/8053-Project/no_npb_10_1000.png", no_npb_10_1000)



### Grpah and limit x to 0 to 100, include confidence intervals, just ebayes + f_bayes
eb_v_hb <- meansdf %>% tidyr::gather("id", "value", c(2:4)) %>%
  # mutate(upper = ifelse(
  #   id =="Emp_Bayes", e_b_conf_u, ifelse(
  #     id == "H_Bayes", h_b_conf_u, ifelse(
  #       id == "IP_Bayes", i_b_conf_u, 0
  #     ))),
  #   lower = ifelse(
  #     id =="Emp_Bayes", e_b_conf_l, ifelse(
  #       id == "H_Bayes", h_b_conf_l, ifelse(
  #         id == "IP_Bayes", i_b_conf_l, 0
  #       )))) %>%
  select(c("rname", "value", "id")) %>%
  ggplot(., aes(rname + 9, value, col = id)) +
  geom_line(aes(linetype = as.factor(id)), size = 1)+
  # geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1) +
  xlim(10, 100) +
  # theme_light() +
  labs(title = "Bayes Methods Simulation Comparrison (10 < n < 100)",
       x = "Training Set Sample Size", y = "MSE", color = "Method",
       linetype = "Method") +
  theme(plot.title = element_text(hjust = 0.25, size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values = c("blue", "red", "brown", "orange", "black"),
                     name = "Method") 
ggsave("~/Documents/GitHub/8053-Project/eb_v_hb_10_100.png", eb_v_hb)

### Graph and lit x to 900 to 1000
eb_v_hb_100_1000 <- meansdf %>% tidyr::gather("id", "value", c(2:4)) %>%
  ggplot(., aes(rname + 10, value, col = id)) +
  geom_line(aes(linetype = as.factor(id))) +
  xlim(100, 1000) +
  ylim(0, 2.5) +
  # theme_light() +
  labs(title = "Bayes Methods Simulation Comparrison (100 < n < 1000)",
       x = "Training Set Sample Size", y = "MSE", color = "Method",
       linetype = "Method") +
  theme(plot.title = element_text(hjust = 0.25, size = 18, face = "bold")) +  ### hjust ~ 1.5
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  scale_color_manual(values = c("blue", "red", "brown", "orange", "black"),
                     name = "Method") 
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



