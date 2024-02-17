
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project")
#source("mixfunctions.R")
library("Clarity")
source("extendedfunctions.R")


                        ####################
                        #                  #
                        #    TO DO LIST    #         
                        #   ------------   #
                        #                  #
                        ####################


## FIX NN_INTERCHANGE BUG: SWAP NODE WITH ADMIX ON IT, NEEDS TO TAKE ADMIX NODE WITH IT (if using nni)

## implement admix target => scale depth proportional to weight?? probably not

## implement topology record to avoid repeat topologies being tested
## then add pruneregraft and mix regraft - can remove nni after?
## pruneregraft and mix regraft => too many possible moves...
## limit it to only prune/regraft edges source and target edges?


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

plotres <- function(gref,dataref,center=FALSE){
  pred <- ccov_dag_(gref)
  loss <- dataref - pred
  pt=plot.cg_(gref,show=FALSE)
  Clarity_Chart(ordermatrix(loss,pt$order),scalefun=myscalefun2,text=T )
}

################################################################################

## Simulate a 4 population model where there is an outgroup
tg0=simCoal_(5,labels=c("A","B","C","D", "O"), alpha = 1,outgroup="O")
## Add to this graph a "mixture edge" from node 1 to node 3
tg1=mixedge_(tg0, 1, 3, 0.5, 0.2)
tg2=mixedge_(tg1,2,3,0.5,0.8) ## Secondly with weight 0.8

par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)
#plot.cg_(tg2,main="Mixture graph 2", showedges = T)

# simulate data
data = ccov_dag_(tg1)

## Make a random pair of alternative graphs
trand0=simCoal_(5,labels=c("A","B","C","D","O"),outgroup="O")
trand1=mixedge_(trand0,2,3,0.5,0.2)  ## Careful not to involve the outgroup
#trand2=mixedge_(trand1,1,3,0.5,0.2)

plot.cg_(trand0, main="Alternative graph",showedges = T )
plot.cg_(trand1,main="Alternative mixture graph", showedges = T)
#plot.cg_(trand2,main="Alternative mixture graph 2", showedges = T)

## Perform inference
tinf1=infer_graph(trand1, data, maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), verbose = T, losstol = 1e-7)

plot.cg_(tg1,main="Data graph", showedges = T)
plot.cg_(tinf1$g, main="Predicted graph")

trand1 = mixedge_(tinf1$g,2,3,0.5,0.5)
plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(trand1, main="Inferred graph with admix",showedges = T )

tinf2=infer_graph(trand1, data, maxiter = 100, movefreqs = c(0,0,1/3),verbose = T)

plot.cg_(tinf2$g, main = "Inferred graph", showedges = T)
plot.cg_(tg1, main="Actual graph",showedges = T )

# Covariance matrices
tinfpred1=ccov_dag_(tinf2[[1]])
tpred1=ccov_dag_(tg1)

## Get the plot that we would make, if we made a plot
pt=plot.cg_(tg1,show=FALSE)

Clarity_Chart(ordermatrix(tpred1,pt$order),scalefun=myscalefun2,text=T )
Clarity_Chart(ordermatrix(tinfpred1,pt$order),scalefun=myscalefun2,text=T)


##################################
set.seed(123)

## Simulate a 9 population model where there is an outgroup
tg0=simCoal_(9,labels=c("A","B","C","D","E","F","G","H","O"), alpha = 1,outgroup="O")
## Add to this graph a "mixture edge" from node 1 to node 3
tg1=mixedge_(tg0, 1, 3, 0.5, 0.2)
#tg2=mixedge_(tg1,4,3,0.5,0.8) ## Secondly with weight 0.8


par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)

# simulate data
data=ccov_dag_(tg1)

# Make a random proposal graph
trand0=simCoal_(9,labels=c("A","B","C","D","E","F","G","H","O"),outgroup="O")
trand1 = mixedge_(trand0,2,3,0.5,0.5)
plot.cg_(trand0, main="Proposal graph",showedges = T )
plot.cg_(trand1, main="Proposal graph 2",showedges = T )

## Perform inference without mixture edges
tinf1=infer_graph(trand0, data, movefreqs = c(1/2,0,1/2) , maxiter = 200 ,verbose =T)

plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(tg1, main="Actual graph",showedges = T )

## Add mixture edge to inferred graph
trand1 = mixedge_(tinf1$g,2,3,0.5,0.5)
plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(trand1, main="Inferred graph with admix",showedges = T )

## Perform inference again
tinf2=infer_graph(tinf2$g, data, movefreqs = c(0,0,1/3) , maxiter = 300 ,verbose =T)

plot.cg_(tinf2$g, main = "Inferred graph", showedges = T)
plot.cg_(tg1, main="Actual graph",showedges = T )

E=edges.cg_(g)
dag = cbind(E$parent,E$child)

graphrecord <- list()

check = any(unlist(lapply(graphrecord, function(entry) all(dag == entry))))

if(!check) dag_list[[length(dag_list) + 1]] = dag




## prune/regraft and mix regraft swap: 2 10 mixnode: 18 mixswap: 5 3

g = trand1
g = removemixedge_(g,i=18)

g = mypruneregraft_(g, source = 2, target = 10)
plot.cg_(g)

proposal=myregraftmixedge_(g, i = NA, 
                           source = 5,
                           target = 3,
                           alpha = runif(1,0,0.5),w = runif(1,0.005, 0.5))  

plot.cg_(proposal)



dag1 <- matrix(c(10, 10, 9, 9, 8, 8, 7, 7, 6, 6, 3, 1, 5, 8, 10, 7, 1, 6, 4, 2), ncol = 2, byrow = TRUE)
dag2 <- matrix(c( 9, 9,10, 10, 8, 8, 7, 7, 6, 6, 3, 1, 5, 8, 10, 7, 1, 6, 4, 2), ncol = 2, byrow = TRUE)


check_proposal(dagstep_(g,data,movetype = 1),graphrecord)$graphrecord


###########################


tinf1=infer_graph(tinf1$g, data, movefreqs = c(1/3,1/3,1/2) , maxiter = 400 ,verbose =T)



nearestTip <- function(g, i){
  tips = tipsunder_(g, i)
  distancelist = list()
  for (t in tips) {
    distance = 0
    nodepath = nodesabove(g,t)
    for (j in nodepath) {
      distance = distance + g$nl[[j]]$d #*g$nl[[j]]$w
    }
    distancelist[[as.character(t)]] = distance
  }
  return(distancelist)
}

plotres <- function(gref,dataref,center=FALSE){
  pred <- ccov_dag_(gref)
  loss <- dataref - pred
  pt=plot.cg_(gref,show=FALSE)
  Clarity_Chart(ordermatrix(loss,pt$order),scalefun=myscalefun2,text=T )
}

plotcov <- function(gref,center=FALSE){
  pred <- ccov_dag_(gref)

  #loss <- (pred-dataref)^2
  pt=plot.cg_(gref,show=FALSE)
  Clarity_Chart(ordermatrix(pred,pt$order),scalefun=myscalefun2,text=T )
}


# Positive residuals indicate pairs of populations where the model underestimates 
# the observed covariance, and thus populations where the fit might be improved by 
# adding additional edges.

# Negative residuals indicate pairs of populations where the model overestimates 
# the observed covariance; these are a necessary outcome of having positive residuals,
# but can also sometimes be interpreted as populations that are forced too close together
# due to unmodeled migration elsewhere in the graph.




################################################################################

g0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")
g1 = mixedge(g0,2,3,0.5,0.2)
plot(g1)

getp(g1)

gparvec(g1)

printDepths(g1)

gparvec(g0)



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

