
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project")
#source("mixfunctions.R")
library("Clarity")
source("extendedfunctions.R")

#if (!"mvnfast" %in% installed.packages()) install.packages("mvnfast")
library(mvnfast)
library(foreach)
library(doParallel)
library(reshape2)
library(ggplot2)

################################################################################

#save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")


## Simulation studies ##

annealing = 0
maxiter = 250


# Simulation 1

pop5_mix0_feat10_10 = runsim(n_features1 = 10,n_features2 = 10,n_pop = 5,n_edges = 0,n_trials = 100)

pop5_mix0_feat25_25 = runsim(n_features1 = 25,n_features2 = 25,n_pop = 5,n_edges = 0,n_trials = 100)

pop5_mix0_feat50_50 = runsim(n_features1 = 50,n_features2 = 50,n_pop = 5,n_edges = 0,n_trials = 100)

pop5_mix0_feat100_100 = runsim(n_features1 = 100,n_features2 = 100,n_pop = 5,n_edges = 0,n_trials = 100)

pop5_mix0_feat250_250 = runsim(n_features1 = 250,n_features2 = 250,n_pop = 5,n_edges = 0,n_trials = 100)

pop5_mix0_feat500_500 = runsim(n_features1 = 500,n_features2 = 500,n_pop = 5,n_edges = 0,n_trials = 100)

pop5_mix0_feat1000_1000 = runsim(n_features1 = 1000,n_features2 = 1000,n_pop = 5,n_edges = 0,n_trials = 100)


# Simulation 2

# n_pop = 5
# n_edges = 1

pop5_mix1_feat10_10 = runsim(n_features1 = 10,n_features2 = 10,n_pop = 5,n_edges = 1,n_trials = 100)

pop5_mix1_feat25_25 = runsim(n_features1 = 25,n_features2 = 25,n_pop = 5,n_edges = 1,n_trials = 100)
save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")
pop5_mix1_feat50_50 = runsim(n_features1 = 50,n_features2 = 50,n_pop = 5,n_edges = 1,n_trials = 100)
save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")
pop5_mix1_feat100_100 = runsim(n_features1 = 100,n_features2 = 100,n_pop = 5,n_edges = 1,n_trials = 100)
save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")
pop5_mix1_feat250_250 = runsim(n_features1 = 250,n_features2 = 250,n_pop = 5,n_edges = 1,n_trials = 100)
save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")
pop5_mix1_feat500_500 = runsim(n_features1 = 500,n_features2 = 500,n_pop = 5,n_edges = 1,n_trials = 100)
save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")
pop5_mix1_feat1000_1000 = runsim(n_features1 = 1000,n_features2 = 1000,n_pop = 5,n_edges = 1,n_trials = 100)
save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")




# Simulation 3

# n_pop = 5
# n_edges = 2

pop5_mix2_feat10_10_add = runsim(n_features1 = 10,n_features2 = 10,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat25_25_add = runsim(n_features1 = 25,n_features2 = 25,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat50_50_add = runsim(n_features1 = 50,n_features2 = 50,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat100_100_add = runsim(n_features1 = 100,n_features2 = 100,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat250_250_add = runsim(n_features1 = 250,n_features2 = 250,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat500_500_add = runsim(n_features1 = 500,n_features2 = 500,n_pop = 5,n_edges = 2,n_trials = 75)




temp10 = c(pop5_mix2_feat10_10, pop5_mix2_feat10_10_add)

temp25 = c(pop5_mix2_feat25_25, pop5_mix2_feat25_25_add)

temp50 = c(pop5_mix2_feat50_50, pop5_mix2_feat50_50_add)

temp100 = c(pop5_mix2_feat100_100, pop5_mix2_feat100_100_add)

temp250 = c(pop5_mix2_feat250_250, pop5_mix2_feat250_250_add)

temp500 = c(pop5_mix2_feat500_500, pop5_mix2_feat500_500_add)

temp1000 = c(pop5_mix2_feat1000_1000,pop5_mix2_feat1000_1000_add1,pop5_mix2_feat1000_1000_add2,pop5_mix2_feat1000_1000_add3)



pop5_mix2_feat10_10 = temp10

pop5_mix2_feat25_25 = temp25

pop5_mix2_feat50_50 = temp50

pop5_mix2_feat100_100 = temp100

pop5_mix2_feat250_250 = temp250

pop5_mix2_feat500_500 = temp500

pop5_mix2_feat1000_1000 = temp1000


# ####
# pop5_mix2_feat1000_1000_add1 = runsim(n_features1 = 1000,n_features2 = 1000,n_pop = 5,n_edges = 2,n_trials = 25)
# save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")
# 
# pop5_mix2_feat1000_1000_add2 = runsim(n_features1 = 1000,n_features2 = 1000,n_pop = 5,n_edges = 2,n_trials = 25)
# save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")
# 
# pop5_mix2_feat1000_1000_add3 = runsim(n_features1 = 1000,n_features2 = 1000,n_pop = 5,n_edges = 2,n_trials = 25)
# save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")



#plot.cg_(pop5_mix2_feat500_500_add[[1]]$true_record$g1)





# Experiment 4

# n_pop = 6
# n_edges = 1

####


## Plots - Simulation 1 ##

oose_pop5_mix0_feat10_10 = outofsample_errors(pop5_mix0_feat10_10)
oose_pop5_mix0_feat25_25 = outofsample_errors(pop5_mix0_feat25_25)
oose_pop5_mix0_feat50_50 = outofsample_errors(pop5_mix0_feat50_50)
oose_pop5_mix0_feat100_100 = outofsample_errors(pop5_mix0_feat100_100)
oose_pop5_mix0_feat250_250 = outofsample_errors(pop5_mix0_feat250_250)
oose_pop5_mix0_feat500_500 = outofsample_errors(pop5_mix0_feat500_500)
oose_pop5_mix0_feat1000_1000 = outofsample_errors(pop5_mix0_feat1000_1000)


generate_ooseplot(oose_pop5_mix0_feat10_10, title = "5 populations, 0 admixtures, 10 features")
generate_ooseplot(oose_pop5_mix0_feat25_25, title = "5 populations, 0 admixtures, 25 features")
generate_ooseplot(oose_pop5_mix0_feat50_50, title = "5 populations, 0 admixtures, 50 features")
generate_ooseplot(oose_pop5_mix0_feat100_100, title = "5 populations, 0 admixtures, 100 features")
generate_ooseplot(oose_pop5_mix0_feat250_250, title = "5 populations, 0 admixtures, 250 features")
generate_ooseplot(oose_pop5_mix0_feat500_500, title = "5 populations, 0 admixtures, 500 features")
generate_ooseplot(oose_pop5_mix0_feat1000_1000, title = "5 populations, 0 admixtures, 1000 features")


create_pairs_plot(oose_pop5_mix0_feat10_10, title = "5 populations, 0 admixtures, 10 features")
create_pairs_plot(oose_pop5_mix0_feat25_25, title = "5 populations, 0 admixtures, 25 features")
create_pairs_plot(oose_pop5_mix0_feat50_50, title = "5 populations, 0 admixtures, 50 features")
create_pairs_plot(oose_pop5_mix0_feat100_100, title = "5 populations, 0 admixtures, 100 features")
create_pairs_plot(oose_pop5_mix0_feat250_250, title = "5 populations, 0 admixtures, 250 features")
create_pairs_plot(oose_pop5_mix0_feat500_500, title = "5 populations, 0 admixtures, 500 features")
create_pairs_plot(oose_pop5_mix0_feat1000_1000, title = "5 populations, 0 admixtures, 1000 features")


data_list_sim1 <- list(oose_pop5_mix0_feat10_10, oose_pop5_mix0_feat25_25,
                       oose_pop5_mix0_feat50_50, oose_pop5_mix0_feat100_100,
                       oose_pop5_mix0_feat250_250, oose_pop5_mix0_feat500_500,
                       oose_pop5_mix0_feat1000_1000)

n_features <- c(10, 25, 50, 100, 250, 500, 1000)

plot_mean_oose_byfeat(n_features, data_list_sim1, title = "5 populations, 0 admixtures")


###

# Simulation 2 processing

oose_pop5_mix1_feat10_10 = outofsample_errors(pop5_mix1_feat10_10)
oose_pop5_mix1_feat25_25 = outofsample_errors(pop5_mix1_feat25_25)
oose_pop5_mix1_feat50_50 = outofsample_errors(pop5_mix1_feat50_50)
oose_pop5_mix1_feat100_100 = outofsample_errors(pop5_mix1_feat100_100)
oose_pop5_mix1_feat250_250 = outofsample_errors(pop5_mix1_feat250_250)
oose_pop5_mix1_feat500_500 = outofsample_errors(pop5_mix1_feat500_500)
oose_pop5_mix1_feat1000_1000 = outofsample_errors(pop5_mix1_feat1000_1000)



generate_ooseplot(oose_pop5_mix1_feat10_10, title = "5 populations, 1 admixture, 10 features")
generate_ooseplot(oose_pop5_mix1_feat25_25, title = "5 populations, 1 admixture, 25 features")
generate_ooseplot(oose_pop5_mix1_feat50_50, title = "5 populations, 1 admixture, 50 features")
generate_ooseplot(oose_pop5_mix1_feat100_100, title = "5 populations, 1 admixture, 100 features")
generate_ooseplot(oose_pop5_mix1_feat250_250, title = "5 populations, 1 admixture, 250 features")
generate_ooseplot(oose_pop5_mix1_feat500_500, title = "5 populations, 1 admixture, 500 features")
generate_ooseplot(oose_pop5_mix1_feat1000_1000, title = "5 populations, 1 admixture, 1000 features")



create_pairs_plot(oose_pop5_mix1_feat10_10, title = "5 populations, 1 admixture, 10 features")
create_pairs_plot(oose_pop5_mix1_feat25_25, title = "5 populations, 1 admixture, 25 features")
create_pairs_plot(oose_pop5_mix1_feat50_50, title = "5 populations, 1 admixture, 50 features")
create_pairs_plot(oose_pop5_mix1_feat100_100, title = "5 populations, 1 admixture, 100 features")
create_pairs_plot(oose_pop5_mix1_feat250_250, title = "5 populations, 1 admixture, 250 features")
create_pairs_plot(oose_pop5_mix1_feat500_500, title = "5 populations, 1 admixture, 500 features")
create_pairs_plot(oose_pop5_mix1_feat1000_1000, title = "5 populations, 1 admixture, 1000 features")

data_list_sim2 <- list(oose_pop5_mix1_feat10_10, oose_pop5_mix1_feat25_25,
                       oose_pop5_mix1_feat50_50, oose_pop5_mix1_feat100_100,
                       oose_pop5_mix1_feat250_250, oose_pop5_mix1_feat500_500,
                       oose_pop5_mix1_feat1000_1000)

n_features <- c(10, 25, 50, 100, 250, 500, 1000)

plot_mean_oose_byfeat(n_features, data_list_sim2, title = "5 populations, 1 admixture")

####

# Simulation 3 processing

oose_pop5_mix2_feat10_10 = outofsample_errors(pop5_mix2_feat10_10)
oose_pop5_mix2_feat25_25 = outofsample_errors(pop5_mix2_feat25_25)
oose_pop5_mix2_feat50_50 = outofsample_errors(pop5_mix2_feat50_50)
oose_pop5_mix2_feat100_100 = outofsample_errors(pop5_mix2_feat100_100)
oose_pop5_mix2_feat250_250 = outofsample_errors(pop5_mix2_feat250_250)
oose_pop5_mix2_feat500_500 = outofsample_errors(pop5_mix2_feat500_500)
oose_pop5_mix2_feat1000_1000 = outofsample_errors(pop5_mix2_feat1000_1000)


generate_ooseplot(oose_pop5_mix2_feat10_10, title = "5 populations, 2 admixtures, 10 features")
generate_ooseplot(oose_pop5_mix2_feat25_25, title = "5 populations, 2 admixtures, 25 features")
generate_ooseplot(oose_pop5_mix2_feat50_50, title = "5 populations, 2 admixtures, 50 features")
generate_ooseplot(oose_pop5_mix2_feat100_100, title = "5 populations, 2 admixtures, 100 features")
generate_ooseplot(oose_pop5_mix2_feat250_250, title = "5 populations, 2 admixtures, 250 features")
generate_ooseplot(oose_pop5_mix2_feat500_500, title = "5 populations, 2 admixtures, 500 features")
generate_ooseplot(oose_pop5_mix2_feat1000_1000, title = "5 populations, 2 admixtures, 1000 features")


create_pairs_plot(oose_pop5_mix2_feat10_10, title = "5 populations, 2 admixtures, 10 features")
create_pairs_plot(oose_pop5_mix2_feat25_25, title = "5 populations, 2 admixtures, 25 features")
create_pairs_plot(oose_pop5_mix2_feat50_50, title = "5 populations, 2 admixtures, 50 features")
create_pairs_plot(oose_pop5_mix2_feat100_100, title = "5 populations, 2 admixtures, 100 features")
create_pairs_plot(oose_pop5_mix2_feat250_250, title = "5 populations, 2 admixtures, 250 features")
create_pairs_plot(oose_pop5_mix2_feat500_500, title = "5 populations, 2 admixtures, 500 features")
create_pairs_plot(oose_pop5_mix2_feat1000_1000, title = "5 populations, 2 admixtures, 1000 features")


data_list_sim3 <- list(oose_pop5_mix2_feat10_10, oose_pop5_mix2_feat25_25,
                       oose_pop5_mix2_feat50_50, oose_pop5_mix2_feat100_100,
                       oose_pop5_mix2_feat250_250, oose_pop5_mix2_feat500_500,
                       oose_pop5_mix2_feat1000_1000)

n_features <- c(10, 25, 50, 100, 250, 500, 1000)

plot_mean_oose_byfeat(n_features, data_list_sim3, title = "5 populations, 2 admixtures")

oose_pop5_mix2_feat10_10[[1:2]]


################################################################################
# 
# data_list_sim2 <- list(oose_pop5_mix2_feat10_10, oose_pop5_mix2_feat25_25,
#                        oose_pop5_mix2_feat50_50, oose_pop5_mix2_feat100_100,
#                        oose_pop5_mix2_feat250_250, oose_pop5_mix2_feat500_500,
#                        oose_pop5_mix2_feat1000_1000)
# 
# n_features <- c(10, 25, 50, 100, 250, 500, 1000)
# 
# plot_mean_oose_byfeat(n_features, data_list_sim2, title = "5 populations, 2 admixtures")
# 



p5m2_10_10 = pop5_mix0_feat10_10
#oose_p5m2_50_50 = outofsample_errors(pop5_mix0_feat10_10)
generate_ooseplot(oose_p5m2_10_10)
create_pairs_plot(oose_p5m2_10_10)




n_features = c(10,25,50,100,250,500,1000)

mean_error1_single_10 = mean(errors_pop5_mix2_feat10_10$res_single_error1)
mean_error1_joint_10 = mean(errors_pop5_mix2_feat10_10$res_joint_error1)

mean_error1_single_25 = mean(errors_pop5_mix2_feat25_25$res_single_error1)
mean_error1_joint_25 = mean(errors_pop5_mix2_feat25_25$res_joint_error1)

mean_error1_single_50 = mean(errors_pop5_mix2_feat50_50$res_single_error1)
mean_error1_joint_50 = mean(errors_pop5_mix2_feat50_50$res_joint_error1)

mean_error1_single_100 = mean(errors_pop5_mix2_feat100_100$res_single_error1)
mean_error1_joint_100 = mean(errors_pop5_mix2_feat100_100$res_joint_error1)

mean_error1_single_250 = mean(errors_pop5_mix2_feat250_250$res_single_error1)
mean_error1_joint_250 = mean(errors_pop5_mix2_feat250_250$res_joint_error1)

mean_error1_single_500 = mean(errors_pop5_mix2_feat500_500$res_single_error1)
mean_error1_joint_500 = mean(errors_pop5_mix2_feat500_500$res_joint_error1)

mean_error1_single_1000 = mean(errors_pop5_mix2_feat1000_1000$res_single_error1)
mean_error1_joint_1000 = mean(errors_pop5_mix2_feat1000_1000$res_joint_error1)

mean_error_single1 = c(mean_error1_single_10, mean_error1_single_25, mean_error1_single_50,
                        mean_error1_single_100, mean_error1_single_250, mean_error1_single_500, 
                        mean_error1_single_1000)

mean_error_joint1 = c(mean_error1_joint_10, mean_error1_joint_25, mean_error1_joint_50,
                       mean_error1_joint_100, mean_error1_joint_250, mean_error1_joint_500,
                       mean_error1_joint_1000)


df1 = data.frame(n_features, mean_error_single1, mean_error_joint1)

df1_long = reshape2::melt(df1, id.vars = "n_features")

ggplot(df1_long, aes(x = n_features, y = value, color = variable)) +
  geom_line() +
  labs(x = "Number of Features", y = "Mean Predictive Error", color = "Method") +
  theme_minimal()


# Standard deviation
std_error_single1 = c(sd(errors_pop5_mix2_feat10_10$res_single_error1), sd(errors_pop5_mix2_feat25_25$res_single_error1), sd(errors_pop5_mix2_feat50_50$res_single_error1),
                      sd(errors_pop5_mix2_feat100_100$res_single_error1), sd(errors_pop5_mix2_feat250_250$res_single_error1), sd(errors_pop5_mix2_feat500_500$res_single_error1), 
                      sd(errors_pop5_mix2_feat1000_1000$res_single_error1))

std_error_joint1 = c(sd(errors_pop5_mix2_feat10_10$res_joint_error1), sd(errors_pop5_mix2_feat25_25$res_joint_error1), sd(errors_pop5_mix2_feat50_50$res_joint_error1),
                     sd(errors_pop5_mix2_feat100_100$res_joint_error1), sd(errors_pop5_mix2_feat250_250$res_joint_error1), sd(errors_pop5_mix2_feat500_500$res_joint_error1),
                     sd(errors_pop5_mix2_feat1000_1000$res_joint_error1))

df1 = data.frame(n_features,std_error_single1, std_error_joint1)

df1_long = reshape2::melt(df1, id.vars = "n_features")

ggplot(df1_long, aes(x = n_features, y = value, color = variable)) +
  geom_line() +
  labs(x = "Number of Features", y = "Error", color = "Method") +
  theme_minimal()

################################################################################

mean_error2_single_10 = mean(errors_pop5_mix2_feat10_10$res_single_error2)
mean_error2_joint_10 = mean(errors_pop5_mix2_feat10_10$res_joint_error2)

mean_error2_single_25 = mean(errors_pop5_mix2_feat25_25$res_single_error2)
mean_error2_joint_25 = mean(errors_pop5_mix2_feat25_25$res_joint_error2)

mean_error2_single_50 = mean(errors_pop5_mix2_feat50_50$res_single_error2)
mean_error2_joint_50 = mean(errors_pop5_mix2_feat50_50$res_joint_error2)

mean_error2_single_100 = mean(errors_pop5_mix2_feat100_100$res_single_error2)
mean_error2_joint_100 = mean(errors_pop5_mix2_feat100_100$res_joint_error2)

mean_error2_single_250 = mean(errors_pop5_mix2_feat250_250$res_single_error2)
mean_error2_joint_250 = mean(errors_pop5_mix2_feat250_250$res_joint_error2)

mean_error2_single_500 = mean(errors_pop5_mix2_feat500_500$res_single_error2)
mean_error2_joint_500 = mean(errors_pop5_mix2_feat500_500$res_joint_error2)

mean_error2_single_1000 = mean(errors_pop5_mix2_feat1000_1000$res_single_error2)
mean_error2_joint_1000 = mean(errors_pop5_mix2_feat1000_1000$res_joint_error2)

mean_error_single2 = c(mean_error2_single_10, mean_error2_single_25, mean_error2_single_50,
                       mean_error2_single_100, mean_error2_single_250, mean_error2_single_500, 
                       mean_error2_single_1000)

mean_error_joint2 = c(mean_error2_joint_10, mean_error2_joint_25, mean_error2_joint_50,
                      mean_error2_joint_100, mean_error2_joint_250, mean_error2_joint_500,
                      mean_error2_joint_1000)


df2 = data.frame(n_features, mean_error_single2, mean_error_joint2)

df2_long = reshape2::melt(df2, id.vars = "n_features")

ggplot(df2_long, aes(x = n_features, y = value, color = variable)) +
  geom_line() +
  labs(x = "Number of Features", y = "Mean Predictive Error", color = "Method") +
  theme_minimal()


# Standard deviation
std_error_single2 = c(sd(errors_pop5_mix2_feat10_10$res_single_error2), sd(errors_pop5_mix2_feat25_25$res_single_error2), sd(errors_pop5_mix2_feat50_50$res_single_error2),
                      sd(errors_pop5_mix2_feat100_100$res_single_error2), sd(errors_pop5_mix2_feat250_250$res_single_error2), sd(errors_pop5_mix2_feat500_500$res_single_error2), 
                      sd(errors_pop5_mix2_feat1000_1000$res_single_error2))

std_error_joint2 = c(sd(errors_pop5_mix2_feat10_10$res_joint_error2), sd(errors_pop5_mix2_feat25_25$res_joint_error2), sd(errors_pop5_mix2_feat50_50$res_joint_error2),
                     sd(errors_pop5_mix2_feat100_100$res_joint_error2), sd(errors_pop5_mix2_feat250_250$res_joint_error2), sd(errors_pop5_mix2_feat500_500$res_joint_error2),
                     sd(errors_pop5_mix2_feat1000_1000$res_joint_error2))

df2 = data.frame(n_features,std_error_single2, std_error_joint2)

df2_long = reshape2::melt(df2, id.vars = "n_features")

ggplot(df2_long, aes(x = n_features, y = value, color = variable)) +
  geom_line() +
  labs(x = "Number of Features", y = "Error", color = "Method") +
  theme_minimal()


# # 10 features each
# errors_pop5_mix2_feat10_10 = prediction_errors(pop5_mix2_feat10_10)
# generate_plots(errors_pop5_mix2_feat10_10, title = "5 populations, 2 admixture, 10 features each")
# 
# # 25 features each
# errors_pop5_mix2_feat25_25 = prediction_errors(pop5_mix2_feat25_25)
# generate_plots(errors_pop5_mix2_feat25_25, title = "5 populations, 2 admixture, 25 features each")
# 
# #  50 features each
# errors_pop5_mix2_feat50_50 = prediction_errors(pop5_mix2_feat50_50)
# generate_plots(errors_pop5_mix2_feat50_50, title = "5 populations, 2 admixture, 50 features each")
# 
# # 100 features each
# errors_pop5_mix2_feat100_100 = prediction_errors(pop5_mix2_feat100_100)
# generate_plots(errors_pop5_mix2_feat100_100, title = "5 populations, 2 admixture, 100 features each")
# 
# # 250 features each
# errors_pop5_mix2_feat250_250 = prediction_errors(pop5_mix2_feat250_250)
# generate_plots(errors_pop5_mix2_feat250_250, title = "5 populations, 2 admixture, 250 features each")
# 
# # 500 features each
# errors_pop5_mix2_feat500_500 = prediction_errors(pop5_mix2_feat500_500)
# generate_plots(errors_pop5_mix2_feat500_500, title = "5 populations, 2 admixture, 500 features each")
# 
# # 1000 features each
# errors_pop5_mix2_feat1000_1000 = prediction_errors(pop5_mix2_feat1000_1000)
# generate_plots(errors_pop5_mix2_feat1000_1000, title = "5 populations, 2 admixture, 1000 features each")

#######
# 5 populations, 2 admixture, dataset1 = 10 features, dataset2 = 50 features"
#errors_pop5_mix2_feat10_50 = prediction_errors(pop5_mix2_feat10_50)
generate_plots(errors_pop5_mix2_feat10_50, title = "5 populations, 2 admixture, dataset1 = 10 features, dataset2 = 50 features")


# 5 populations, 2 admixture, dataset1 = 10 features, dataset2 = 1000 features"
#errors_pop5_mix2_feat10_1000 = prediction_errors(pop5_mix2_feat10_1000)
generate_plots(errors_pop5_mix2_feat10_1000, title = "5 populations, 2 admixture, dataset1 = 10 features, dataset2 = 1000 features")


# 5 populations, 1 admixture, 25 features each
errors_pop5_mix1_feat25_25 = prediction_errors(pop5_mix1_feat25_25)
generate_plots(errors_pop5_mix1_feat25_25, title = "5 populations, 1 admixture, 25 features each")

#######

oose_pop5_mix2_feat10_10 = outofsample_errors(pop5_mix2_feat10_10)
generate_ooseplot(oose_pop5_mix2_feat10_10)

# 5 populations, 2 admixture, 25 features each
oose_pop5_mix2_feat25_25 = outofsample_errors(pop5_mix2_feat25_25)
generate_ooseplot(oose_pop5_mix2_feat25_25)

# 5 populations, 2 admixture, 50 features each
oose_pop5_mix2_feat50_50 = outofsample_errors(pop5_mix2_feat50_50)
generate_ooseplot(oose_pop5_mix2_feat50_50)

# 5 populations, 2 admixture, 100 features each
oose_pop5_mix2_feat100_100 = outofsample_errors(pop5_mix2_feat100_100)
generate_ooseplot(oose_pop5_mix2_feat100_100)

# 5 populations, 2 admixture, 250 features each
oose_pop5_mix2_feat250_250 = outofsample_errors(pop5_mix2_feat250_250)
generate_ooseplot(oose_pop5_mix2_feat250_250)

# 5 populations, 2 admixture, 500 features each
oose_pop5_mix2_feat500_500 = outofsample_errors(pop5_mix2_feat500_500)
generate_ooseplot(oose_pop5_mix2_feat500_500)

# 5 populations, 2 admixture, 1000 features each
oose_pop5_mix2_feat1000_1000 = outofsample_errors(pop5_mix2_feat1000_1000)
generate_ooseplot(oose_pop5_mix2_feat1000_1000)

means_oose_joint = c(
  mean(oose_pop5_mix2_feat10_10$res_joint_error),
  mean(oose_pop5_mix2_feat25_25$res_joint_error),
  mean(oose_pop5_mix2_feat50_50$res_joint_error),
  mean(oose_pop5_mix2_feat100_100$res_joint_error),
  mean(oose_pop5_mix2_feat250_250$res_joint_error),
  mean(oose_pop5_mix2_feat500_500$res_joint_error),
  mean(oose_pop5_mix2_feat1000_1000$res_joint_error)
)

means_oose_data1 = c(
  mean(oose_pop5_mix2_feat10_10$res_single_error1),
  mean(oose_pop5_mix2_feat25_25$res_single_error1),
  mean(oose_pop5_mix2_feat50_50$res_single_error1),
  mean(oose_pop5_mix2_feat100_100$res_single_error1),
  mean(oose_pop5_mix2_feat250_250$res_single_error1),
  mean(oose_pop5_mix2_feat500_500$res_single_error1),
  mean(oose_pop5_mix2_feat1000_1000$res_single_error1)
)

means_oose_data2 = c(
  mean(oose_pop5_mix2_feat10_10$res_single_error2),
  mean(oose_pop5_mix2_feat25_25$res_single_error2),
  mean(oose_pop5_mix2_feat50_50$res_single_error2),
  mean(oose_pop5_mix2_feat100_100$res_single_error2),
  mean(oose_pop5_mix2_feat250_250$res_single_error2),
  mean(oose_pop5_mix2_feat500_500$res_single_error2),
  mean(oose_pop5_mix2_feat1000_1000$res_single_error2)
)


oose_df = data.frame(n_features, means_oose_joint, means_oose_data1,means_oose_data2)

oose_df_long = reshape2::melt(oose_df, id.vars = "n_features")


ggplot(oose_df_long, aes(x = n_features, y = value, color = variable)) +
  geom_line() +
  #scale_y_log10() +
  labs(x = "Number of Features", y = "Mean Out-of-Sample Predictive Error", color = "Method") +
  theme_minimal()




################################################################################

# Pairs log-log plots (Out of Sample)




#########



# 10 features
create_pairs_plot(oose_pop5_mix2_feat10_10)

# 25 features
create_log_log_plot(oose_pop5_mix2_feat25_25)

# 50 features
create_log_log_plot(oose_pop5_mix2_feat50_50)

# 100 features
create_log_log_plot(oose_pop5_mix2_feat100_100)

# 250 features
create_log_log_plot(oose_pop5_mix2_feat250_250)

# 500 features
create_log_log_plot(oose_pop5_mix2_feat500_500)

# 1000 features
create_pairs_plot(oose_pop5_mix2_feat1000_1000)







################################################################################

# pred_pop5_mix2_feat100000_10 = prediction_errors(pop5_mix2_feat100000_10,n_trials = 10)
# 
# 
# pred_pop5_mix1_feat500 = prediction_errors(sim1_feat500 ,n_trials = 50)
# 
# 
# generate_plots(pred_pop5_mix2_feat100000_10,n_trials = 10,title = "10,000 features")
# 
# 
# generate_plots(pred_pop5_mix1_feat500,n_trials = 50,title = "5 populations, 1 admixture, 500 features each")
# 
# 
# 
# 
# 
# generate_plots(results1,n_trials = 50,title = "10,000 features")
# 
# 
# generate_plots(sim1_feat500,n_trials = 50,title = "500 features")
# 
# generate_plots(results4,n_trials = 50,title = "500 features")
# 
# 
# 
# generate_plots(results5,n_trials = 50,title = "100 features")
# 
# generate_plots(results7,n_trials = 50,title = "50 features")
# 
# generate_plots(results6,n_trials = 50,title = "25 features")
# 
# 
# generate_plots(results8,n_trials = 50,title = "10 features")
# 
# 
# 
# 
# 
