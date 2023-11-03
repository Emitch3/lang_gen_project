

setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project")
source("mixfunctions.R")
source("extendedfunctions.R")
library("Clarity")

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


################################################################################

tg0 <- simCoal_(4,labels=c("A","B","C","O"),outgroup="O")
tg1 <- mixedge_(tg0, 3,2,0.5,1)
tg2 <- removemixedge_(tg1,i = 7)

par(mfrow=c(2,2))
plot.cg_(tg0, main="Original graph", )
plot.cg_(tg1,main="Mixture graph")
plot.cg_(tg2, main = "Mixture edge removed")



## Add to this graph a "mixture edge" from node 2 to node 3
tg1=mixedge_(tg0,2,3,0.5,0.2) ## Firstly with weight 0.2
tg2=mixedge_(tg0,2,3,0.5,0.8) ## Seconldly with weight 0.8


cov(tg0)

ccov_dag(g)

