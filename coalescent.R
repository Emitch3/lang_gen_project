
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project")
source("mixfunctions.R")
library("Clarity")

## TO DO
# Fix weights for removeedges_ ?


# DONE
# Forbid mixtures between siblings i.e populations with same parent node - done
# Implement posterior distribution for number of children a node has for simcoal, with alpha parameter - not needed anymore
# Update edges.cg_ to set admix lengths to 0 - done

################################################################################

## Simulate a 4 population model where there is an outgroup
tg0 = simCoal_(4,labels=c("A","B","C", "O"), alpha = 1,outgroup="O")
## Add to this graph a "mixture edge" from node 1 to node 3
tg1 <- mixedge_(tg0, 1, 3, 0.5, 0.2)

par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)

# simulate data
data=ccov_dag_(tg1)

## Make a random pair of alternative graphs
trand0=simCoal_(4,labels=c("A","B","C","O"),outgroup="O")
trand1=mixedge_(trand0,2,3,0.5,0.5)  ## Careful not to involve the outgroup

plot.cg_(trand0, main="Alternative graph",showedges = T )
plot.cg_(trand1,main="Alternative mixture graph", showedges = T)

## Perform inference
tinf1=infer_topo(trand1, data, maxiter = 50 ,verbose = T, patience = 1, losstol = 0.0005)
plot.cg_(tg1,main="Data graph", showedges = T)
plot.cg_(tinf1[[1]], main="Predicted graph")
data
ccov_dag_(tinf1[[1]])
losses = tinf1[[3]]
losses


plot.cg_(tg1)
x = dagstep_(tg1, data, freqs = c(0,0,1), verbose = T)
plot.cg_(x)

y = dagstep_(tg1, data, freqs = c(0,1,0), verbose = T)
plot.cg_(y)

dagstep_(trand1, data, freqs = c(0,0,0,1))

################################################################################

g0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")

## Add to this graph a "mixture edge" from node 2 to node 3
g1=mixedge(g0,2,3,0.5,0.2) ## Firstly with weight 0.2
g2=mixedge(g0,2,3,0.5,0.8) ## Seconldly with weight 0.8

par(mfrow=c(1,2))
plot.cg(g0, main="Original graph",showedges = T )
plot.cg(g1,main="Mixture graph", showedges = T)

## Simulate covariances from these two models
pred1=ccov_dag(g1)
pred2=ccov_dag(g2)

## Make a random pair of alternative graphs
rand0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")
rand1=mixedge(rand0,2,3,0.5,0.5)  ## Careful not to involve the outgroup

plot.cg(rand0, main="Alternative graph",showedges = T )
plot.cg(rand1,main="Alternative mixture graph", showedges = T)

## Perform inference
inf1=infer_dag(rand1,pred1)

plot.cg(g1,main="Data graph", showedges = T)
plot.cg(inf1[[1]], main="Predicted graph")



################################################################################

printType <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"type:", g$nl[[i]]$type))
  }
}


printFamily <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"children:", g$nl[[i]]$children, "parent:",g$nl[[i]]$parents))
  }
}

printDepths <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i, "parent node", g$nl[[i]]$parents, "d:", g$nl[[i]]$d))
  }
}

printTimes <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"t:", g$nl[[i]]$t))
  }
}

printPos <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"p:", g$nl[[i]]$p))
  }
}

printTmp <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"tmp:", g$nl[[i]]$tmp))
  }
}

printType <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"type:", g$nl[[i]]$type))
  }
}

printWeight <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"weight:", g$nl[[i]]$w))
  }
}


printPos <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"p:", g$nl[[i]]$p))
  }
}

printTmp <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"tmp:", g$nl[[i]]$tmp))
  }
}
