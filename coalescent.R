
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project")
source("mixfunctions.R")
library("Clarity")
source("extendedfunctions.R")

# include admix in get_depths - get_weightmatrix (at bottom of this page), change format_params tpars in parametrise 

ordermatrix<-function(x,order){
  ## Reorder a matrix x ro the order given
  ## For plotting
  x=x[rownames(x)%in%order,colnames(x)%in%order]
  mat=matrix(NA,nrow=length(order),ncol=length(order))
  rownames(mat)=colnames(mat)=order
  mat[rownames(x),colnames(x)]=x
  mat
}

myscalefun2=function(x)sign(x)*log(abs(1+x))

################################################################################

## Simulate a 4 population model where there is an outgroup
tg0=simCoal_(5,labels=c("A","B","C","D", "O"), alpha = 1,outgroup="O")
## Add to this graph a "mixture edge" from node 1 to node 3
tg1=mixedge_(tg0, 1, 3, 0.5, 0.2)
#tg2=mixedge_(tg1,4,3,0.5,0.8) ## Seconldly with weight 0.8


par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)

# simulate data
data=ccov_dag_(tg1)

## Make a random pair of alternative graphs
trand0=simCoal_(5,labels=c("A","B","C","D","O"),outgroup="O")
trand1=mixedge_(trand0,2,3,0.5,0.5)  ## Careful not to involve the outgroup

plot.cg_(trand0, main="Alternative graph",showedges = T )
plot.cg_(trand1,main="Alternative mixture graph", showedges = T)

## Perform inference
tinf1=infer_graph(trand1, data, maxiter = 50 ,verbose = T, losstol = 0.001, initial_temperature = 0.5, cooling_factor = 1.5) 

tinf1=infer_graph(tinf1[[1]], data, maxiter = 100 ,verbose = T, losstol = 0.001, initial_temperature = 0.01, cooling_factor = 1) 

plot.cg_(tg1,main="Data graph", showedges = T)
plot.cg_(tinf1[[1]], main="Predicted graph")
data
ccov_dag_(tinf1[[1]])
losses = tinf1[[3]]
losses
min(losses)

# Covariance matrices
tinfpred1=ccov_dag_(tinf1[[1]])
tpred1=ccov_dag_(tg1)

## Get the plot that we would make, if we made a plot
pt=plot.cg_(tg1,show=FALSE)

Clarity_Chart(ordermatrix(tpred1,pt$order),scalefun=myscalefun2,text=T )
Clarity_Chart(ordermatrix(tinfpred1,pt$order),scalefun=myscalefun2,text=T)



##################################

## Simulate a 4 population model where there is an outgroup
tg0=simCoal_(5,labels=c("A","B","C","D", "O"), alpha = 1,outgroup="O")
## Add to this graph a "mixture edge" from node 1 to node 3
tg1=mixedge_(tg0, 1, 3, 0.5, 0.2)
#tg2=mixedge_(tg1,4,3,0.5,0.8) ## Seconldly with weight 0.8


par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)

# simulate data
data=ccov_dag_(tg1)

## Make a random pair of alternative graphs
trand0=simCoal_(5,labels=c("A","B","C","D","O"),outgroup="O")
trand1=mixedge_(trand0,2,3,0.5,0.5)  ## Careful not to involve the outgroup

plot.cg_(trand0, main="Alternative graph",showedges = T )
plot.cg_(trand1,main="Alternative mixture graph", showedges = T)

## Perform inference
#tinf1=infer_graph_hc(trand1, data, maxiter = 100 ,verbose = T, losstol = 0.001, patience = 1) 

tinf1=infer_graph(tinf1$g, data, maxiter = 50 ,verbose = T, losstol = 0.001, initial_temperature = 0.25, cooling_factor = 1) 

plot.cg_(tg1,main="Data graph", showedges = T)
plot.cg_(tinf1[[1]], main="Predicted graph")
data
ccov_dag_(tinf1[[1]])
losses = tinf1[[3]]
losses
min(losses)

# Covariance matrices
tinfpred1=ccov_dag_(tinf1[[1]])
tpred1=ccov_dag_(tg1)

## Get the plot that we would make, if we made a plot
pt=plot.cg_(tg1,show=FALSE)

Clarity_Chart(ordermatrix(tpred1,pt$order),scalefun=myscalefun2,text=T )
Clarity_Chart(ordermatrix(tinfpred1,pt$order),scalefun=myscalefun2,text=T)





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


########################################

# 

get_weightmatrix = function(g) {
  edges = edges.cg_(g)
  p = list()
  for(i in g$tips) {p[[i]] = paths(edges, i, g$root)}
  # make dictionary mapping edges to column numbers
  dict = list()
  l = 1
  mixnodes = g$mix
  for(i in 1:length(g$nl)){
    node = g$nl[[i]]
    # if(i %in% mixnodes)
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
      }
    }
    W = rbind(W,unlist(rowvec))
  }
  colnames(W) = sort(names(dict))
  return(W)
}
