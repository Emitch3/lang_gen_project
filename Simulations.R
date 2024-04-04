
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


