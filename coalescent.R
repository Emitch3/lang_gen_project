
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project")
source("mixfunctions.R")
library("Clarity")

## TO DO
# Implement posterior distribution for number of children a node has for simcoal, with alpha parameter
# Fix weights for removeedges_
# Update edges.cg_ to set admix lengths to 0

################################################################################

tg0 <- simCoal_(4,labels=c("A","B","C","O"), times = "test",outgroup="O")
tg1 <- mixedge_(tg0, 1, 2, 0.5, 0.2)
tg2<- mixedge_(tg1, 1, 2, 0.5, 0.3)

par(mfrow=c(1,2))
plot.cg_(tg0)
plot.cg_(tg1)
plot.cg_(tg2)

par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)
plot.cg_(tg2, main = "Mixture graph 2")

g0 <- simCoal(4,labels=c("A","B","C","O"), times = "coal",outgroup="O")
g1 <- mixedge(g0, 1, 2, 0.5, 0.2)

ccov_dag_(tg0)
ccov_dag_(tg1)

ccov_dag_(g0,old = T)

ccov_dag(g1)
ccov_dag_(g1,old = T)

plot(g1)

ccov_dag_(tg0)
ccov_dag_(tg1)
ccov_dag_(tg2)



get_weightmatrix = function(g) {
  edges = edges.cg_(g)
  p = NA
  for(i in g$tips) {p[[i]] = paths(edges, i, g$root)}
  
  # make dictionary mapping edges to column numbers
  dict = list()
  l = 1
  for(i in 1:length(g$nl)){
    node = g$nl[[i]]
    if(node$id == g$root) next
    parents = node$parents
    for(j in 1:length(parents)){
      parent = parents[j]
      edge =  paste0(node$id,"-",parent)
      dict[as.character(edge)] = l
      l = l + 1
    }
  }
  
  W = NULL  
  for(i in g$tips){
    rowvec = NULL  
    for(j in 1:(length(g$nl)-1 + length(g$mix))) rowvec[j] = 0  
    
    for(k in 1:length(p[[i]])){
      for (l in 1:length(p[[i]])){
        
        overlap = p[[i]][[k]]$edge[ p[[i]][[k]]$edge %in% p[[i]][[l]]$edge] 

        index = as.numeric(dict[overlap])
        rowvec[index] = as.numeric(rowvec[index]) + unique(p[[i]][[k]]$weight)*unique(p[[i]][[l]]$weight)
        #print(c(overlap, index, unique(p[[i]][[k]]$weight)*unique(p[[i]][[l]]$weight)))
      }
    }
    W = rbind(W,unlist(rowvec))
  }
  colnames(W) = sort(names(dict))
  return(W)
}

V = ccov_dag_(tg1)
V
W = get_weightmatrix(tg1)
W

find_lengths = function(g){
  V = ccov_dag_(g)
  W = get_weightmatrix(g)
  sol = t(W)%*%solve(W%*%t(W))%*%diag(V)
  return(sol)
}

find_lengths(tg1)

sol

edges.cg_(tg1)

U = c(1.5,1,0,1,6,2,3,1.5)

W%*%U
V

################################################################################

printType <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"type:", g$nl[[i]]$type))
  }
}

diag(V, nrow = 4, ncol = 4)
printFamily <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"children:", g$nl[[i]]$children, "parent:",g$nl[[i]]$parents))
  }
}

printDepths <- function(g){
  for(i in 1:length(g$nl)){
    print(paste("node:",i,"d:", g$nl[[i]]$d))
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
