
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


## Load simulation data :
#load("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Simulation 1 data.RData")


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

pop5_mix1_feat50_50 = runsim(n_features1 = 50,n_features2 = 50,n_pop = 5,n_edges = 1,n_trials = 100)

pop5_mix1_feat100_100 = runsim(n_features1 = 100,n_features2 = 100,n_pop = 5,n_edges = 1,n_trials = 100)

pop5_mix1_feat250_250 = runsim(n_features1 = 250,n_features2 = 250,n_pop = 5,n_edges = 1,n_trials = 100)

pop5_mix1_feat500_500 = runsim(n_features1 = 500,n_features2 = 500,n_pop = 5,n_edges = 1,n_trials = 100)

pop5_mix1_feat1000_1000 = runsim(n_features1 = 1000,n_features2 = 1000,n_pop = 5,n_edges = 1,n_trials = 100)




# Simulation 3

# n_pop = 5
# n_edges = 2

pop5_mix2_feat10_10_add = runsim(n_features1 = 10,n_features2 = 10,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat25_25_add = runsim(n_features1 = 25,n_features2 = 25,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat50_50_add = runsim(n_features1 = 50,n_features2 = 50,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat100_100_add = runsim(n_features1 = 100,n_features2 = 100,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat250_250_add = runsim(n_features1 = 250,n_features2 = 250,n_pop = 5,n_edges = 2,n_trials = 75)

pop5_mix2_feat500_500_add = runsim(n_features1 = 500,n_features2 = 500,n_pop = 5,n_edges = 2,n_trials = 75)



####


## Simulation 1 processing ##

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


par(mfrow=c(2,2))
generate_oosehist(oose_pop5_mix0_feat10_10, title = "5 populations, 0 admixtures, 10 features")
generate_oosehist(oose_pop5_mix0_feat25_25, title = "5 populations, 0 admixtures, 25 features")
generate_oosehist(oose_pop5_mix0_feat50_50, title = "5 populations, 0 admixtures, 50 features")
generate_oosehist(oose_pop5_mix0_feat100_100, title = "5 populations, 0 admixtures, 100 features")
generate_oosehist(oose_pop5_mix0_feat250_250, title = "5 populations, 0 admixtures, 250 features")
generate_oosehist(oose_pop5_mix0_feat500_500, title = "5 populations, 0 admixtures, 500 features")
generate_oosehist(oose_pop5_mix0_feat1000_1000, title = "5 populations, 0 admixtures, 1000 features")



create_pairs_plot(oose_pop5_mix0_feat10_10, title = "5 populations, 0 admixtures, 10 features")
create_pairs_plot(oose_pop5_mix0_feat25_25, title = "5 populations, 0 admixtures, 25 features")
create_pairs_plot(oose_pop5_mix0_feat50_50, title = "5 populations, 0 admixtures, 50 features")
create_pairs_plot(oose_pop5_mix0_feat100_100, title = "5 populations, 0 admixtures, 100 features")
create_pairs_plot(oose_pop5_mix0_feat250_250, title = "5 populations, 0 admixtures, 250 features")

par(mfrow=c(1,2))
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



generate_oosehist(oose_pop5_mix1_feat10_10, title = "5 populations, 1 admixture, 10 features")
generate_oosehist(oose_pop5_mix1_feat25_25, title = "5 populations, 1 admixture, 25 features")
generate_oosehist(oose_pop5_mix1_feat50_50, title = "5 populations, 1 admixture, 50 features")
generate_oosehist(oose_pop5_mix1_feat100_100, title = "5 populations, 1 admixture, 100 features")
generate_oosehist(oose_pop5_mix1_feat250_250, title = "5 populations, 1 admixture, 250 features")
generate_oosehist(oose_pop5_mix1_feat500_500, title = "5 populations, 1 admixture, 500 features")
generate_oosehist(oose_pop5_mix1_feat1000_1000, title = "5 populations, 1 admixture, 1000 features")


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




generate_oosehist(oose_pop5_mix2_feat10_10, title = "5 populations, 2 admixtures, 10 features")
generate_oosehist(oose_pop5_mix2_feat25_25, title = "5 populations, 2 admixtures, 25 features")
generate_oosehist(oose_pop5_mix2_feat50_50, title = "5 populations, 2 admixtures, 50 features")
generate_oosehist(oose_pop5_mix2_feat100_100, title = "5 populations, 2 admixtures, 100 features")
generate_oosehist(oose_pop5_mix2_feat250_250, title = "5 populations, 2 admixtures, 250 features")
generate_oosehist(oose_pop5_mix2_feat500_500, title = "5 populations, 2 admixtures, 500 features")
generate_oosehist(oose_pop5_mix2_feat1000_1000, title = "5 populations, 2 admixtures, 1000 features")




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

####

data_list_bymix <- list(oose_pop5_mix0_feat1000_1000, 
                       oose_pop5_mix1_feat1000_1000, oose_pop5_mix2_feat1000_1000)

calculate_se(data_list_bymix,"res_joint_error")

n_admix = c(0,1,2)


plot_mean_oose_bymix(n_admix ,data_list_bymix, title = " (1000 features)")






