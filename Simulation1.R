
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

poplabels =c("A","B","C","D", "E","F","G","H","I","J","K","L","M","N")
maxiter = 200
# Set number of populations
n_pop = 5
labels = c(poplabels[1:n_pop-1],"O")
# Set number of mixture edges
n_edges = 1 


# Set number of features 
n_features = 100

# Parallelisation 
cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

# Number of trials
n_trials = 50



results5 = foreach(k = 1:n_trials) %dopar% {
  library(mvnfast)
  
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
  
  # Compute true data losses
  pg1 = evaluate_proposal(trueglist$g1, truedatalist$g1)
  g1loss = ctree_loss_(pg1,truedatalist$g1)
  
  pg2 = evaluate_proposal(trueglist$g2, truedatalist$g2)
  g2loss = ctree_loss_(pg2,truedatalist$g2)
  
  pglist = list(pg1 = pg1,pg2 = pg2)
  class(pglist) = "cglist"
  
  jointloss = ctree_loss_(pglist,truedatalist)
  
  # Record true observations
  true_record <- list(
    g1 = pg1,
    g2 = pg2,
    g1loss = g1loss,
    g2loss = g2loss,
    jointloss = jointloss,
    truedata1 = truedatalist$g1,
    truedata2 = truedatalist$g2
  )
  
  # Generate "noisy" data from the true data
  simdata1 <- cov(rmvn(n = n_features, mu = rep(0, n_pop), sigma = truedatalist$g1))
  colnames(simdata1) <- labels
  rownames(simdata1) <- labels
  
  simdata2 <- cov(rmvn(n = n_features, mu = rep(0, n_pop), sigma = truedatalist$g2))
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
                           verbose = FALSE, plotprogress = FALSE, losstol = jointloss, MHtol = 0.5)
  
  # Inference with data set 1
  data_inf1 <- infer_graph(g = proplist$p1, data =  simdatalist$g1, maxiter = maxiter, movefreqs = c(1/4, 1/4, 1/4, 1/4),
                           verbose = FALSE, plotprogress = FALSE, losstol = g1loss, MHtol = 0.5)
  
  # Inference with data set 2
  data_inf2 <- infer_graph(g = proplist$p1, data =  simdatalist$g2, maxiter = maxiter, movefreqs = c(1/4, 1/4, 1/4, 1/4),
                           verbose = FALSE, plotprogress = FALSE, losstol = g2loss, MHtol = 0.5)
  
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

# result1 = 10000 features
# result8 = 10 features
# result7 = 50 features

# res8_joint_error = NULL
# res8_single_error1 = NULL
# res8_single_error2 = NULL
# 
# for(i in 1:n_trials){
#   res8_joint_error[i] =  (results8[[i]]$true_record$jointloss - results8[[i]]$sim_joint_results$loss)^2
#   res8_single_error1[i] = (results8[[i]]$true_record$g1loss - results8[[i]]$sim_data1_results$loss)^2
#   res8_single_error2[i] = (results8[[i]]$true_record$g2loss - results8[[i]]$sim_data2_results$loss)^2
# }


generate_plots <- function(results, n_trials, title=NULL) {
  res_joint_error = NULL
  res_single_error1 = NULL
  res_single_error2 = NULL
  
  for(i in 1:n_trials){
    res_joint_error[i] =  (0.5*results[[i]]$true_record$jointloss - 0.5*results[[i]]$sim_joint_results$loss)^2
    res_single_error1[i] = (results[[i]]$true_record$g1loss - results[[i]]$sim_data1_results$loss)^2
    res_single_error2[i] = (results[[i]]$true_record$g2loss - results[[i]]$sim_data2_results$loss)^2
  }
  
  df <- data.frame(
    Method = rep(c("Joint Error", "Single Error 1", "Single Error 2"), each = n_trials),
    Error = c(res1_joint_error, res_single_error1, res_single_error2)
  )
  
  iterations <- 1:n_trials
  
  df1 = data.frame(iterations, res_joint_error, res_single_error1, res_single_error2)
  
  df1_long = melt(df1, id.vars = "iterations", variable.name = "method", value.name = "error")
  
  p1 <- ggplot(df1_long, aes(x = iterations, y = error, color = method, shape = method, group = iterations)) +
    geom_point(size=2.5 ) +
    geom_line(aes(group = iterations), color = "gray", linetype = "dashed") +
    scale_shape_manual(values=c(16, 17, 18)) +
    scale_y_log10() +
    labs(title = paste0("Predictive Error Comparison"," - ", title), x = "Iterations", y = "Log(Error)") +
    theme_minimal()
  
  p2 <- ggplot(df1_long, aes(x = method, y = error, fill = method)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = paste0("Predictive Error Comparison"," - ", title), x = "Method", y = "Log(Error)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  list(plot1=p1, plot2=p2,  res_joint_error = res_joint_error,
       res_single_error1 = res_single_error1,
       res_single_error2 = res_single_error2)
}


plots = generate_plots(results8,n_trials = 50,title = )

plots$plot1
plots$plot2
plots$res_joint_error
plots$res_single_error1
plots$res_single_error2


plots = generate_plots(results7,n_trials = 50)

plots$plot1
plots$plot2
plots$res_joint_error
plots$res_single_error1
plots$res_single_error2


generate_plots(results1,n_trials = 50,title = "10,000 features")

generate_plots(results7,n_trials = 50,title = "50 features")
generate_plots(results8,n_trials = 50,title = "10 features")

generate_plots(results6,n_trials = 50,title = "50 features")


generate_plots(results5,n_trials = 50,title = "25 features")








################################################################################

res1_joint_error = NULL
res1_single_error1 = NULL
res1_single_error2 = NULL

for(i in 1:n_trials){
  res1_joint_error[i] =  (results1[[i]]$true_record$jointloss - results1[[i]]$sim_joint_results$loss)^2
  res1_single_error1[i] = (results1[[i]]$true_record$g1loss - results1[[i]]$sim_data1_results$loss)^2
  res1_single_error2[i] = (results1[[i]]$true_record$g2loss - results1[[i]]$sim_data2_results$loss)^2
}


# Combine the vectors into a data frame
df <- data.frame(
  Method = rep(c("Joint Error", "Single Error 1", "Single Error 2"), each = length(joint_error)),
  Error = c(res1_joint_error, res1_single_error1, res1_single_error2)
)


iterations <- 1:n_trials

df1 = data.frame(iterations, res1_joint_error, res1_single_error1, res1_single_error2)

df1_long = melt(df1, id.vars = "iterations", variable.name = "method", value.name = "error")


# Create a scatter plot using ggplot2 with a logarithmic scale on the y-axis
ggplot(df1_long, aes(x = iterations, y = error, color = method, shape = method, group = iterations)) +
  geom_point(size=2.5 ) +
  geom_line(aes(group = iterations), color = "gray", linetype = "dashed") +  # Add gray dashed lines connecting the points vertically
  scale_shape_manual(values=c(16, 17, 18)) +
  scale_y_log10() +  # Use a logarithmic scale for the y-axis
  labs(title = "Predictive Error Comparison", x = "Iterations", y = "Log(Error)") +
  theme_minimal()



# Create box plots using ggplot2 with a logarithmic scale on the y-axis
ggplot(df1_long, aes(x = method, y = error, fill = method)) +
  geom_boxplot() +
  scale_y_log10() +  # Use a logarithmic scale for the y-axis
  labs(title = "Predictive Error Comparison", x = "Method", y = "Log(Error)") +
  theme_minimal() +
  theme(legend.position = "none") 

################################################################################






################################################################################ 

# for(k in 1:n_trials) {
# 
# 
#   # Simulate a population model where there is an outgroup
#   g0=simCoal_(n_pop,labels=labels, outgroup="O")
#   
#   # Simulate a new graph with same topology but different parameters 
#   g1 = reparameterise(g0)
#   
#   # Add both to list
#   trueglist0 = list()
#   trueglist0$g1 = g0
#   trueglist0$g2 = g1
#   class(trueglist0) = "cglist"
#   
#   # Add mixture edges to both graphs
#   trueglist = add_mixedges(g= trueglist0,n_edges = n_edges)
#   
#   # Plot both graphs
#   plot.cg_(trueglist$g1,main = "True graph 1")
#   plot.cg_(trueglist$g2,main = "True graph 2")
#   
#   # Compute covariance matrices from both graphs
#   truedatalist = list()
#   truedatalist = lapply(trueglist, ccov_dag_)
#   
#   
#   # Simulate proposal graph for inference
#   proposal0 = simCoal_(n_pop,labels=labels, outgroup="O")
#   proposal = add_mixedges(proposal0,n_edges = n_edges)
#   proplist = list()
#   proplist$p1 = proposal
#   proplist$p2 = proposal
#   class(proplist) = "cglist"
#   
#   
#   # Make predictions with the true data
#   
#   # Joint inference with both data sets
#   joint_inf = infer_graph(g = proplist, data =  truedatalist, maxiter = 1, movefreqs = c(1/4, 1/4, 1/4, 1/4),
#                           verbose = F, plotprogress = F, losstol = 1e-30, MHtol = 0.5)
#   
#   # Inference with data set 1
#   data_inf1 = infer_graph(g = proplist$p1, data =  truedatalist$g1, maxiter = 1, movefreqs = c(1/4, 1/4, 1/4, 1/4),
#                           verbose = F, plotprogress = F, losstol = 1e-30, MHtol = 0.5)
#   
#   # Inference with data set 2
#   data_inf2 = infer_graph(g = proplist$p1, data =  truedatalist$g2, maxiter = 1, movefreqs = c(1/4, 1/4, 1/4, 1/4),
#                           verbose = F, plotprogress = F, losstol = 1e-30, MHtol = 0.5)
#   
#   # Record results
#   true_joint_results[[k]] = list()
#   true_joint_results[[k]]$loss = joint_inf$loss
#   true_joint_results[[k]]$loss_list = joint_inf$loss_list
#   true_joint_results[[k]]$g1 = joint_inf$g[[1]]
#   true_joint_results[[k]]$g2 = joint_inf$g[[2]]
#   
#   true_data1_results[[k]] = list()
#   true_data1_results[[k]]$loss = data_inf1$loss
#   true_data1_results[[k]]$loss_list = data_inf1$loss_list
#   true_data1_results[[k]]$g = data_inf1$g
#   
#   
#   true_data2_results[[k]] = list()
#   true_data2_results[[k]]$loss = data_inf2$loss
#   true_data2_results[[k]]$loss_list = data_inf2$loss_list
#   true_data2_results[[k]]$g = data_inf2$g
#   
#   
#   # Generate "noisy" data from the true data
#   
#   simdata1 = cov(rmvn(n = n_features, mu = rep(0, n_pop), sigma = truedatalist$g1 ))
#   colnames(simdata1) = labels
#   rownames(simdata1) = labels
#   
#   simdata2 = cov(rmvn(n = n_features, mu = rep(0, n_pop), sigma = truedatalist$g2 ))
#   colnames(simdata2) = labels
#   rownames(simdata2) = labels
#   
#   # Add both simulated data sets to list
#   simdatalist = list()
#   simdatalist$g1 = simdata1
#   simdatalist$g2 = simdata2
#   
#   
#   # Simulate proposal graph for inference
#   proposal0 = simCoal_(n_pop,labels=labels, outgroup="O")
#   proposal = add_mixedges(proposal0,n_edges = n_edges)
#   proplist = list()
#   proplist$p1 = proposal
#   proplist$p2 = proposal
#   class(proplist) = "cglist"
#   
#   # Make predictions with the simulated data
#   
#   # Joint inference with both data sets
#   joint_inf = infer_graph(g = proplist, data =  simdatalist, maxiter = 1, movefreqs = c(1/4, 1/4, 1/4, 1/4),
#                           verbose = F, plotprogress = F, losstol = 1e-30, MHtol = 0.5)
#   
#   # Inference with data set 1
#   data_inf1 = infer_graph(g = proplist$p1, data =  simdatalist$g1, maxiter = 1, movefreqs = c(1/4, 1/4, 1/4, 1/4),
#                           verbose = F, plotprogress = F, losstol = 1e-30, MHtol = 0.5)
#   
#   # Inference with data set 2
#   data_inf2 = infer_graph(g = proplist$p1, data =  simdatalist$g2, maxiter = 1, movefreqs = c(1/4, 1/4, 1/4, 1/4),
#                           verbose = F, plotprogress = F, losstol = 1e-30, MHtol = 0.5)
#   
#   
#   
#   # Record results
#   
#   sim_joint_results[[k]] = list()
#   sim_joint_results[[k]]$loss = joint_inf$loss
#   sim_joint_results[[k]]$loss_list = joint_inf$loss_list
#   sim_joint_results[[k]]$g1 = joint_inf$g[[1]]
#   sim_joint_results[[k]]$g2 = joint_inf$g[[2]]
#   
#   sim_data1_results[[k]] = list()
#   sim_data1_results[[k]]$loss = data_inf1$loss
#   sim_data1_results[[k]]$loss_list = data_inf1$loss_list
#   sim_data1_results[[k]]$g = data_inf1$g
#   
#   sim_data2_results[[k]] = list()
#   sim_data2_results[[k]]$loss = data_inf2$loss
#   sim_data2_results[[k]]$loss_list = data_inf2$loss_list
#   sim_data2_results[[k]]$g = data_inf2$g
# 
# }
# 
# 
# (true_joint_results[[1]]$loss - sim_joint_results[[1]]$loss)^2
# (true_data1_results[[1]]$loss - sim_data1_results[[1]]$loss)^2
# (true_data2_results[[1]]$loss - sim_data2_results[[1]]$loss)^2
# 


################################################################################
















































set.seed(12)
par(mfrow=c(1,2))


## Simulate a 5 population model where there is an outgroup
tg0=simCoal_(5,labels=c("A","B","C","D", "O"), outgroup="O")
plot.cg_(tg0)

tg1=simCoal_(5,labels=c("A","B","C","D", "O"), outgroup="O")
plot.cg_(tg1)

# tg1 = reparameterise(tg0)
# plot.cg_(tg1)

tg2=mixedge_(tg0,2,3,0.5,0.2) ## Firstly with weight 0.2
tg3=mixedge_(tg1,2,3,0.5,0.5) ## Secondly with weight 0.5

tg4=mixedge_(tg2,1,4,0.5,0.3) ## Firstly with weight 0.3
tg5=mixedge_(tg3,1,4,0.5,0.4) ## Secondly with weight 0.4

plot.cg_(tg4)
plot.cg_(tg5)


cglist1 = list(tg4, tg5)
class(cglist1) = "cglist"

dataref = lapply(cglist1, ccov_dag_)

## Generate proposal graph

trand0=simCoal_(5,labels=c("A","B","C","D","O"), outgroup = "O")
trand1 = add_mixedges(trand0, n_edges = 2)
plot.cg_(trand1)


trandlist=list(trand1,trand1)
class(trandlist)="cglist"


################################################################################

## Simulation study 1: Number of features same for both graphs ##


# trials for each data set
n_trials = 10 

features = c(10, 20, 50, 100, 500, 1000, 2000, 10000) 
simdatalist1 = list()
simdatalist2 = list()

for(i in 1:length(features)){
  simdata = cov(rmvn(n = features[i], mu = rep(0, 5), sigma = dataref[[1]] ))
  colnames(simdata) = colnames(dataref[[1]])
  rownames(simdata) = rownames(dataref[[1]])
  
  simdata2 = cov(rmvn(n = features[i], mu = rep(0, 5), sigma = dataref[[2]] ))
  colnames(simdata2) = colnames(dataref[[1]])
  rownames(simdata2) = rownames(dataref[[1]])
  
  simdatalist1[[i]] = simdata
  simdatalist2[[i]] = simdata2
  
}


# Generate proposal list 
proplist = list()
for(i in 1:n_trials){
  trand0 = simCoal_(5, labels = c("A", "B", "C", "D", "O"), outgroup = "O")
  trand1 = add_mixedges(trand0, n_edges = 2)
  trandlist = list(trand1,trand1)
  class(trandlist) = "cglist"
  proplist[[i]] = trandlist
}


jointlist=list()
singlelist1=list()
singlelist2=list()
cl <- makeCluster(detectCores())
registerDoParallel(cl)
  
for (k in 1:n_trials) {

  jointres = foreach(i = 1:n_trials) %dopar% {
    
    # Apply the infer_graph function
    trandlist = list(proplist[[i]][[1]],proplist[[i]][[1]])
    dataref = list(simdatalist1[[k]], simdatalist2[[k]])
    class(trandlist) = "cglist"
    
    tinf1 = infer_graph(g = trandlist, data =  dataref, maxiter = 200, movefreqs = c(1/4, 1/4, 1/4, 1/4),
                        verbose = F, plotprogress = F, losstol = 0.0001, MHtol = 0)
    tinf1
  }

  singleres1 = foreach(i = 1:n_trials) %dopar% {
    
    # Apply the infer_graph function
    tinf2 = infer_graph(g = proplist[[i]][[1]], data = simdatalist1[[k]], maxiter = 200,
                        movefreqs = c(1/4, 1/4, 1/4, 1/4),
                        verbose = F, plotprogress = F, losstol = 0.0001, MHtol = 0)
    tinf2
  }

  singleres2 = foreach(i = 1:n_trials) %dopar% {
    
    # Apply the infer_graph function
    tinf3 = infer_graph(g = proplist[[i]][[1]],data = simdatalist2[[k]], maxiter = 200,
                        movefreqs = c(1/4, 1/4, 1/4, 1/4),
                        verbose = F, plotprogress = F, losstol = 0.0001, MHtol = 0)
    tinf3
  }
  
  jointlist[[k]] = jointres
  singlelist1[[k]] = singleres1
  singlelist2[[k]] = singleres2
}

stopCluster(cl)


## Run time = 4 hours 15 mins

jointlist[[8]][[10]]$loss

single_results2[[7]]$loss

gtrue = cglist1[[1]]



results_table = data.frame(
  Noise_Level = 1:8,
  Joint = sapply(1:8, function(noise) sum(sapply(1:10, function(i) check_topology(jointlist[[noise]][[i]]$g, gtrue)))),
  Single1 = sapply(1:8, function(noise) sum(sapply(1:10, function(i) check_topology(singlelist1[[noise]][[i]]$g, gtrue)))),
  Single2 = sapply(1:8, function(noise) sum(sapply(1:10, function(i) check_topology(singlelist2[[noise]][[i]]$g, gtrue))))
)

results_table


long_data <- melt(results_table, id.vars = "Noise_Level")


ggplot(long_data, aes(x = Noise_Level, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_continuous(breaks = 1:length(features), labels = features) +  
  labs(x = "Number of features", y = "Number of correct inferrences", fill = "Method") +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.1))+
  ggtitle("")



# Perform paired t-tests
ttest_joint_single1 <- t.test(results_table$Joint, results_table$Single1, paired = TRUE)
ttest_joint_single2 <- t.test(results_table$Joint, results_table$Single2, paired = TRUE)
ttest_single1_single2 <- t.test(results_table$Single1, results_table$Single2, paired = TRUE)


print(paste("Joint vs Single1: ", ttest_joint_single1$p.value))
print(paste("Joint vs Single2: ", ttest_joint_single2$p.value))
print(paste("Single1 vs Single2: ", ttest_single1_single2$p.value))






tinf2 = infer_graph(g = proplist[[1]][[1]], data = cglist1[[1]], maxiter = 1,
                    movefreqs = c(1/4, 1/4, 1/4, 1/4),
                    verbose = F, plotprogress = F, losstol = 0.0001, MHtol = 0)
tinf2



################################################################################


