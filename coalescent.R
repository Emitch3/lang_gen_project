
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project")
source("mixfunctions.R")
library("Clarity")

## TO DO
# Implement posterior distribution for number of children a node has for simcoal, with alpha parameter
# Fix weights for removeedges_
# Update edges.cg_ to set admix lengths to 0

################################################################################

tg0 <- simCoal_(4,labels=c("A","B","C","O"), times = "coal",outgroup="O")
tg1 <- mixedge_(tg0, 3, 2, 0.5, 0.2)
tg2<- mixedge_(tg1, 1, 2, 0.5, 0.3)

par(mfrow=c(1,2))
plot.cg_(tg1)
plot.cg_(tg2)

par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)
plot.cg_(tg2, main = "Mixture graph 2")



g0 <- simCoal(4,labels=c("A","B","C","O"), times = "coal",outgroup="O")
g1 <- mixedge(g0, 3, 2, 0.5, 0.2)

ccov_dag_(tg0)
ccov_dag_(tg1)

ccov_dag_(g0,old = T)

system.time(ccov_dag_(g1, old = T))
system.time(ccov_dag(g1))

ccov_dag(g1)
ccov_dag_(g1,old = T)


plot(g1)

ccov_dag_(tg0)
ccov_dag_(tg1)
ccov_dag_(tg2)


g = tg1
root = g$root
ntips = g$n
p = list(length = ntips)
edges = edges.cg_(g)

for(i in 1:ntips){
  p[[i]] = paths(edges, i,root)
}
p

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
