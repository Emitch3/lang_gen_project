
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


## Simulation study 1 ##


runsim <- function(n_features1, n_features2, n_pop, n_edges, n_trials, maxiter=250){

  poplabels =c("A","B","C","D", "E","F","G","H","I","J","K","L","M","N")
  #maxiter = 250
  # Set number of populations
  # n_pop = 5
  labels = c(poplabels[1:n_pop-1],"O")
  # Set number of mixture edges
  # n_edges = 2
  
  # Set number of features 
  # n_features1 = 50
  # n_features2 = 50
  
  
  # Parallelisation 
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # Number of trials
  n_trials = 25
  
  annealing = 0
  
  pop5_mix2_feat50_50 = foreach(k = 1:n_trials) %dopar% {
    library(mvnfast)
    setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project")
    source("extendedfunctions.R")
    
    # Simulate a population model where there is an outgroup
    g0 <- simCoal_(n_pop, labels = labels, outgroup = "O")
    
    # Simulate a new graph with same topology but different parameters 
    g1 <- reparameterise(g0)
    
    # Add both to list
    trueglist0 <- list()
    trueglist0$g1 <- g0
    trueglist0$g2 <- g1
    class(trueglist0) <- "cglist"
    
    # Add mixture edges to both graphs
    trueglist <- add_mixedges(g = trueglist0, n_edges = n_edges)
    
    # Compute covariance matrices from both graphs
    truedatalist <- lapply(trueglist, ccov_dag_)
    
    # Record true observations
    true_record <- list(
      g1 = trueglist$g1,
      g2 = trueglist$g2,
      truedata1 = truedatalist$g1,
      truedata2 = truedatalist$g2
    )
    
    # Generate "noisy" synthetic data from the true data
    simdata1 <- cov(rmvn(n = n_features1, mu = rep(0, n_pop), sigma = truedatalist$g1))
    colnames(simdata1) <- labels
    rownames(simdata1) <- labels
    
    simdata2 <- cov(rmvn(n = n_features2, mu = rep(0, n_pop), sigma = truedatalist$g2))
    colnames(simdata2) <- labels
    rownames(simdata2) <- labels
    
    # Add both simulated data sets to list
    simdatalist <- list()
    simdatalist$g1 <- simdata1
    simdatalist$g2 <- simdata2
    
    # Simulate proposal graph for inference
    proposal0 <- simCoal_(n_pop, labels = labels, outgroup = "O")
    proposal <- add_mixedges(proposal0, n_edges = n_edges)
    proplist <- list()
    proplist$p1 <- proposal
    proplist$p2 <- proposal
    class(proplist) <- "cglist"
    
    # Make predictions with the simulated data
    
    # Joint inference with both data sets
    joint_inf <- infer_graph(g = proplist, data =  simdatalist, maxiter = maxiter, movefreqs = c(1/4, 1/4, 1/4, 1/4),
                             verbose = FALSE, plotprogress = FALSE, losstol = 1e-30, MHtol = annealing)
    
    # Inference with data set 1
    data_inf1 <- infer_graph(g = proplist$p1, data =  simdatalist$g1, maxiter = maxiter, movefreqs = c(1/4, 1/4, 1/4, 1/4),
                             verbose = FALSE, plotprogress = FALSE, losstol = 1e-30, MHtol = annealing)
    
    # Inference with data set 2
    data_inf2 <- infer_graph(g = proplist$p1, data =  simdatalist$g2, maxiter = maxiter, movefreqs = c(1/4, 1/4, 1/4, 1/4),
                             verbose = FALSE, plotprogress = FALSE, losstol = 1e-30, MHtol = annealing)
    
    # Record results
    sim_joint_results <- list(
      loss = joint_inf$loss,
      loss_list = joint_inf$loss_list,
      g1 = joint_inf$g[[1]],
      g2 = joint_inf$g[[2]]
    )
    
    sim_data1_results <- list(
      loss = data_inf1$loss,
      loss_list = data_inf1$loss_list,
      g = data_inf1$g
    )
    
    sim_data2_results <- list(
      loss = data_inf2$loss,
      loss_list = data_inf2$loss_list,
      g = data_inf2$g
    )
    
    return(list(true_record=true_record, sim_joint_results=sim_joint_results,
                sim_data1_results=sim_data1_results,sim_data2_results=sim_data2_results))
  }
  stopCluster(cl)
}

test = runsim(n_features1 = 1000,n_features2 = 1000,n_pop = 5,n_edges = 1,n_trials = 1,maxiter = 1)

## Record ##

# populations = 5, mix = 1, data1 features = 25, data2 features = 25
# populations = 5, mix = 2, data1 features = 10, data2 features = 50

# populations = 5, mix = 2, data1 features = 10, data2 features = 10
# populations = 5, mix = 2, data1 features = 25, data2 features = 25

# populations = 5, mix = 2, data1 features = 10, data2 features = 1000
# populations = 5, mix = 2, data1 features = 100, data2 features = 100
# populations = 5, mix = 2, data1 features = 250, data2 features = 250



# 5 populations, 2 admixture, 10 features each
errors_pop5_mix2_feat10_10 = prediction_errors(pop5_mix2_feat10_10)
generate_plots(errors_pop5_mix2_feat10_10, title = "5 populations, 2 admixture, 10 features each")

# 5 populations, 2 admixture, 25 features each
errors_pop5_mix2_feat25_25 = prediction_errors(pop5_mix2_feat25_25)
generate_plots(errors_pop5_mix2_feat25_25, title = "5 populations, 2 admixture, 25 features each")

# 5 populations, 2 admixture, 50 features each
errors_pop5_mix2_feat50_50 = prediction_errors(pop5_mix2_feat50_50)
generate_plots(errors_pop5_mix2_feat50_50, title = "5 populations, 2 admixture, 50 features each")

# 5 populations, 2 admixture, 100 features each
errors_pop5_mix2_feat100_100 = prediction_errors(pop5_mix2_feat100_100)
generate_plots(errors_pop5_mix2_feat100_100, title = "5 populations, 2 admixture, 100 features each")

# 5 populations, 2 admixture, 250 features each
errors_pop5_mix2_feat250_250 = prediction_errors(pop5_mix2_feat250_250)
generate_plots(errors_pop5_mix2_feat250_250, title = "5 populations, 2 admixture, 250 features each")

# 5 populations, 2 admixture, 500 features each
errors_pop5_mix2_feat500_500 = prediction_errors(pop5_mix2_feat500_500)
generate_plots(errors_pop5_mix2_feat500_500, title = "5 populations, 2 admixture, 500 features each")

# 5 populations, 2 admixture, 1000 features each
errors_pop5_mix2_feat1000_1000 = prediction_errors(pop5_mix2_feat1000_1000)
generate_plots(errors_pop5_mix2_feat1000_1000, title = "5 populations, 2 admixture, 1000 features each")

####

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


create_log_log_plot <- function(errors, title=NULL) {
  df = data.frame(
  joint= (errors$res_joint_error),
  average_data1_data2 = ((errors$res_single_error1 + errors$res_single_error2)/2)
  )
  
  ggplot(df, aes(x = joint, y = average_data1_data2)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    theme_minimal() +
    annotate("text", x = Inf, y = -Inf,
             label = paste("Gradient: ", round(coef(lm(average_data1_data2 ~ joint, data = df))[2], 2)),
             hjust = 1, vjust = 0) +
    labs(x = "Average single log error", y = "Joint log Error", title = "Log-Log plot")
}


#########



# 10 features
create_log_log_plot(oose_pop5_mix2_feat10_10)

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
create_log_log_plot(oose_pop5_mix2_feat1000_1000)







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
