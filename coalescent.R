
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project")
#source("mixfunctions.R")
library("Clarity")
source("extendedfunctions.R")


                        ####################
                        #                  #
                        #      NOTES       #
                        #   -----------    #
                        #                  #
                        ####################



ordermatrix<-function(x,order){
  #Reorder a matrix x ro the order given
  #For plotting
  x=x[rownames(x)%in%order,colnames(x)%in%order]
  mat=matrix(NA,nrow=length(order),ncol=length(order))
  rownames(mat)=colnames(mat)=order
  mat[rownames(x),colnames(x)]=x
  mat
}

myscalefun2<-function(x)sign(x)*log(abs(1+x))

plotres <- function(gref,dataref,center=FALSE){
  pred <- ccov_dag_(gref)
  loss <- dataref - pred
  pt=plot.cg_(gref,show=FALSE)
  Clarity_Chart(ordermatrix(loss,pt$order),scalefun=myscalefun2,text=T )
}

################################################################################

set.seed(12)
par(mfrow=c(1,2))


## Simulate a 4 population model where there is an outgroup
tg0=simCoal_(5,labels=c("A","B","C","D", "O"), outgroup="O")
plot.cg_(tg0)

tg1=simCoal_(5,labels=c("A","B","C","D", "O"), outgroup="O")
plot.cg_(tg1)

tg1 = reparameterise(tg0)
plot.cg_(tg1)


tg2=mixedge_(tg0,2,3,0.5,0.2) ## Firstly with weight 0.2
tg3=mixedge_(tg1,2,3,0.5,0.5) ## Seconldly with weight 0.5

tg4=mixedge_(tg2,1,4,0.5,0.3) ## Firstly with weight 0.3
tg5=mixedge_(tg3,1,4,0.5,0.4) ## Seconldly with weight 0.4

# tg6=mixedge_(tg4,2,1,0.5,0.1)
# tg7=mixedge_(tg5,2,1,0.5,0.2)
# 
# tg2=mixedge_(tg0,3,2,0.5,0.2) ## Firstly with weight 0.1
# tg3=mixedge_(tg1,3,2,0.5,03) ## Seconldly with weight 0.5
# 
# tg4=mixedge_(tg2,1,2,0.5,0.4)
# tg5=mixedge_(tg3,1,2,0.5,0.4)


plot.cg_(tg4)
plot.cg_(tg5)
cglist1 = list(tg4, tg5)
class(cglist1) = "cglist"
dataref = lapply(cglist1, ccov_dag_)

cglist0 = list(tg0, tg1)
class(cglist0) = "cglist"

## Add "mixture edges"
cglist1 = add_mixedges(cglist0, n_edges=2, randomlistweights = T)
par(mfrow=c(1,2))
lapply(cglist1, plot.cg_)
# Simulate covariance matrices from these two models
dataref = lapply(cglist1, ccov_dag_)


# data=ccov_dag_(tg4)
# data2=ccov_dag_(tg5)
# dataref = list(data,data2)


trand0=simCoal_(5,labels=c("A","B","C","D","O"), outgroup = "O")
trand1 = add_mixedges(trand0, n_edges = 2)
plot.cg_(trand1)


## And of the two random graphs used for inference
trandlist=list(trand1,trand1)
class(trandlist)="cglist"

plot.new()
tinf1=infer_graph(trandlist, dataref, maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                  verbose = T, plotprogress = T, losstol = 1e-10, MHtol = 0)  # initial_temp = 1,cooling_rate = 0.99


pt=plot.cg_(tinf1$g[[1]],show = FALSE)
tpred1=ccov_dag_(tinf1$g)

Clarity_Chart(ordermatrix(tpred1[[2]],pt$order),scalefun=myscalefun2,text=T,main="Observed 1")


            
            

check_topology(g=tinf1$g, gtrue = cglist1[[1]])

tinf2=infer_graph(trandlist[[1]], dataref[[1]], maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                  verbose = T, plotprogress = T, losstol = 1e-10, MHtol = 0)


check_topology(g=tinf2$g, gtrue = cglist1[[1]])

tinf3=infer_graph(trandlist[[2]], dataref[[2]], maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                  verbose = T, plotprogress = T, losstol = 1e-10, MHtol = 0)



check_topology(g=tinf3$g, gtrue = cglist1[[1]])


plot.cg_(cglist1[[1]])
plot.cg_(cglist1[[2]])


install.packages("mvnfast")
library(mvnfast)



# Plot results
par(mfrow=c(2,2))
lapply(tinf1$g, plot.cg_)
lapply(cglist1, plot.cg_)

adjacency_matrix(tinf1$g[[1]]) == adjacency_matrix(cglist1[[1]])


g = tinf1$g[[1]]

for (o in 1:1000) {
  if(!all(enumerate(tg0)== enumerate(reparameterise(tg0)))) print("Fail")

}








################################################################################



trandlist1 = add_mixedges(tinf1$g,n_edges = 2,randomlistweights = T)
class(trandlist1) = "cglist"

plot.new()
tinf2=infer_graph(trandlist1, dataref, maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                  verbose = T, plotprogress = T, losstol = 1e-10, MHtol = 0.5)  # initial_temp = 1,cooling_rate = 0.99

# "proposal: movetype: reverse mixture swap: 5 7 mixnode: 13"

g = tinf1$g[[1]]
plot.cg_(g)

g1 =mixedge_(g, source = 1, target = 9, weight= 0.07,alpha = 0.5)

g2 =mixedge_(g1, source = 4, target = 7, weight= 0.27,alpha = 0.5)
plot.cg_(g2)

evaluate_proposal(g2)


data = dataref[[1]]

get_depths(g2,data)

get_weightmatrix(g2)


#get_weight_optim(g2, data,maxit=20,factr = 0.01,pgtol = 0.01)

#reversemixture_(g2,mix = 13)

#dagstep2(g2,dataref[[1]],movetype = 4,moveparams = c(NA,NA,13,5,7))

tinf2=infer_graph(g2, dataref[[1]], maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                  verbose = T, plotprogress = F, losstol = 1e-10, MHtol = 0.5)  # initial_temp = 1,cooling_rate = 0.99

#testg = g2


# Check if the correct topology has been inferred
if(all(adjacency_matrix(tinf1$g[[1]]) ==  adjacency_matrix(cglist1[[1]])) == TRUE){
  print("True topology inferred!")
} else{
  print("True topology not inferred")
}




# Plot results
par(mfrow=c(2,2))
lapply(tinf2$g, plot.cg_)
lapply(cglist1, plot.cg_)




# plot.cg_(tg4, main = "Actual graph 1")
# plot.cg_(tg5, main = "Actual graph 2")


trandz = mixedge_(tinf1$g[[1]],3,4,0.5,0.2)
plot.cg_(trandz)

trandlist=list(trandz,trandz)
class(trandlist)="cglist"

tinf2=infer_graph(trandlist, dataref, maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                  verbose = T, plotprogress = T, losstol = 0.001, MHtol = 0.5)


par(mfrow=c(2,2))
lapply(tinf2$g, plot.cg_)
plot.cg_(tg4, main = "Actual graph 1")
plot.cg_(tg5, main = "Actual graph 2")



################################################################################


# install.packages("MVN")
# library(MVN)

install.packages("mvnfast")
library(mvnfast)


infer_graph2 <- function(g, data, n_edges, maxiter=100, movefreqs=c(1/4,1/4,1/4,1/4), 
                         verbose = FALSE, plotprogress = FALSE, losstol = 1e-30, MHtol = 0){
  # Input Tree proposal (without mixture edges) and number of mixture edges to be inferred
  # Mixture edges are added sequentially to help inference
  
  maxiter_i = ceiling(maxiter/(n_edges+1))
  loss_list = NULL
  if(verbose) print("INITIAL RUN: Iteration 1 - Performing inference without mixture edges")
  
  for(i in 0:n_edges){
    max = maxiter_i*(i+1)
    
    if(i>0) {
      if(verbose)  print(paste("NEW RUN: Iteration", i+1, "- Adding one additional mixture edge to the graph"))
      
      ## Use residuals to suggest mixture parent and child nodes
      res = residual(g=g,data=dataref)
      if(is(res,"list")) {
        i = which.max(sapply(res, max))
        res = res[[i]]
      }
      pair = as.vector(which(res == max(res), arr.ind = TRUE)[1,])
      swap = list(rem=NA,source=pair[1],target=pair[2])
      
      if(is(g,"cglist") && isvalidregraftmixture_(g[[1]],rem = NA,ret = pair)) {
        g = add_mixedges(g, n_edges = 1,swap = swap)
        
      }else if(isvalidregraftmixture_(g,rem = NA,ret = pair)){
        g = add_mixedges(g, n_edges = 1,swap = swap)
        
      }else g = add_mixedges(g, n_edges = 1)
    }
    
    tinf = infer_graph(g, data, maxiter=max, movefreqs=movefreqs , verbose = verbose
                       , plotprogress = plotprogress, losstol = losstol, MHtol = MHtol) 
    g = tinf$g
    res = g
    loss = tinf$loss
    loss_list = c(loss_list, tinf$loss_list)
  }
  return(list(g = g, loss = loss, loss_list = loss_list))
}  


residual <- function(g, data){
  if(is(g,"cglist")){ 
    ret = list()
    for(i in 1:length(g))
      ret[[i]] = residual(g[[i]],data[[i]])
    return(ret)
  }else {
    res = data - ccov_dag_(g)
  }
  return(res)
}


trandlist0=list(trand0,trand0)
class(trandlist0)="cglist"

#add_mixedges(tg0,n_edges=1,swap=c(1,4))



infer_graph3 <- function(g, data, n_edges, maxiter=100, movefreqs=c(1/4,1/4,1/4,1/4), 
                         verbose = FALSE, plotprogress = FALSE, losstol = 1e-30, MHtol = 0){
  # Input Tree proposal (without mixture edges) and number of mixture edges to be inferred
  # Mixture edges are added sequentially to help inference
  
  #maxiter_i = ceiling(maxiter/(n_edges+1))
  loss_list = NULL
  if(verbose) print("INITIAL RUN: Iteration 1 - Performing inference without mixture edges")
  
  for(i in 0:1){
    #max = maxiter_i*(i+1)
    
    if(i==1) {
      if(verbose)  print("NEW RUN: Adding mixture edges")#print(paste("NEW RUN: Iteration", i+1, "- Adding one additional mixture edge to the graph"))
      
      for(j in 1:n_edges) {
        ## Use residuals to suggest mixture parent and child nodes
        res = residual(g=g,data=dataref)
        if(is(res,"list")) {
          k = which.max(sapply(res, max))
          res = res[[k]]
        }
        pair = as.vector(which(res == max(res), arr.ind = TRUE)[1,])
        swap = list(rem=NA,source=pair[1],target=pair[2])
        
        if(is(g,"cglist") && isvalidregraftmixture_(g[[1]],rem = NA,ret = pair)) {
          g = add_mixedges(g, n_edges = 1,swap = swap)
          
        }else if(isvalidregraftmixture_(g,rem = NA,ret = pair)) {
          g = add_mixedges(g, n_edges = 1,swap = swap)
          
        }else g = add_mixedges(g, n_edges = 1)
        
      }
    }
    
    tinf = infer_graph(g, data, maxiter=maxiter, movefreqs=movefreqs , verbose = verbose
                       , plotprogress = plotprogress, losstol = losstol, MHtol = MHtol) 
    g = tinf$g
    res = g
    loss = tinf$loss
    loss_list = c(loss_list, tinf$loss_list)
  }
  return(list(g = g, loss = loss, loss_list = loss_list))
}  


tinf2 = infer_graph3(trandlist0, dataref,n_edges=2, maxiter = 200,verbose = T,
             plotprogress = T,movefreqs = c(1/4,1/4,1/4,1/4), losstol = 1e-10, MHtol = 0.5)



res =  ccov_dag_(cglist1[[1]]) - ccov_dag_(tinf2$g[[1]])
diag(res)=0
res

pairs = which(res == max(res), arr.ind = TRUE)

pairs[1,]


rownames(pairs)

dataref[[1]] - ccov_dag_(tinf2$g[[1]]) 
    
res = residual(g=g,data=dataref)
if(is(res,"list")) {
  i = which.max(sapply(res, max))
  res = res[[i]]
}

swap = which(res == max(res), arr.ind = TRUE)


enumerate(g) %in% g$tips == enumerate(tg5) %in% tg5$tips
g = tg0
plot.cg_(tg0)
get_depths(tg0,data)



  
tg = reparameterise(g)
plot.cg_(tg)

printDepths(tg)

printDepths(g)


printTimes(tg)

printTimes(g)

n_features1=200

n_features2=20

results = list()

for(i in 1:1){
  
  simdata = cov(rmvn(n=n_features1, mu = rep(0, 5) ,sigma = data))
  
  simdata2 = cov(rmvn(n=n_features2, mu = rep(0, 5) ,sigma = data2))
  
  dataref = list(simdata,simdata2)
  
  trand0=simCoal_(5,labels=c("A","B","C","D","O"), outgroup = "O")
  
  trand1 = add_mixedges(trand0, n_edges = 2)
  
  class(trandlist)="cglist"
  

  tinf1=infer_graph(trandlist, dataref, maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                    verbose = T, plotprogress = T, losstol = 0.001, MHtol = 0.5 )
  
  results[[i]] = tinf1
  
}





features = seq(from = 10, to = 100, by = 10)

joint_results = list()
single_results_1 = list()
single_results_2 = list()

dataref = lapply(cglist1, ccov_dag_)
j=1

for (i in seq_along(features)) {
  
  # Generate the data with varying noise
  simdata = cov(rmvn(n = features[i], mu = rep(0, 5), sigma = dataref[[1]] ))
  colnames(simdata) = colnames(dataref[[1]])
  rownames(simdata) = rownames(dataref[[1]])
  
  simdata2 = cov(rmvn(n = features[i], mu = rep(0, 5), sigma = dataref[[2]] ))
  colnames(simdata2) = colnames(dataref[[1]])
  rownames(simdata2) = rownames(dataref[[1]])
  
  simdataref = list(simdata, simdata2)

  trand0 = simCoal_(5, labels = c("A", "B", "C", "D", "O"), outgroup = "O")
  trand1 = add_mixedges(trand0, n_edges = 2)
  class(trandlist) = "cglist"
  
  for(k in 1:5){
    # Apply the infer_graph function
    tinf1 = infer_graph(trandlist, simdataref, maxiter = 2, movefreqs = c(1/4, 1/4, 1/4, 1/4), 
                        verbose = T, plotprogress = T, losstol = 0.0001, MHtol = 0)
    
    tinf2 = infer_graph(trandlist[[1]], simdata, maxiter = 2, movefreqs = c(1/4, 1/4, 1/4, 1/4), 
                        verbose = T, plotprogress = T, losstol = 0.0001, MHtol = 0)
    
    tinf3 = infer_graph(trandlist[[2]], simdata2, maxiter = 2, movefreqs = c(1/4, 1/4, 1/4, 1/4), 
                        verbose = T, plotprogress = T, losstol = 0.0001, MHtol = 0)
    
    # Store the results
    joint_results[[j]][[k]] = tinf1
    single_results1[[j]][[k]] = tinf2
    single_results2[[j]][[k]] = tinf3
    j=j+1
  }
}



l = list()

l[[1]] = 3

l[[10]] = 4

vv = list(tinf1)



xg = nn_interchange(tg3, source= 4, target = 1)

enumerate(xg)
enumerate(tg3)

plot.cg_(xg)
plot.cg_(tg3)

adjacency_matrix(xg) ==
  adjacency_matrix(tg3)

g$nl[[11]]$children

plot.cg_(g)
plot.cg_(tg7)

g = tinf1$g[[2]]

enumerate(g)
enumerate(tg7)

all(adjacency_matrix(g) ==
  adjacency_matrix(tg7))




#g = altg


all(adjacency_matrix(tinf1$g[[1]]) == adjacency_matrix(tg5))


enumerate(tg5)
enumerate(altg)

plot.cg_(altg)
plot.cg_(tg3)

adjacency_matrix(altg)
adjacency_matrix(tg5)


plot.cg_(trand2)




################################################################################



## Make a random pair of alternative graphs
trand0=simCoal_(5,labels=c("A","B","C","D","O"), outgroup = "O")
trand1=mixedge_(trand0,2,3,0.5,0.2)  ## Careful not to involve the outgroup
trand2=mixedge_(trand1,1,3,0.5,0.2)
trand3=mixedge_(trand2,3,4,0.5,0.2)

plot.cg_(trand0, main="Alternative graph",showedges = T)
#plot.cg_(trand1,main="Alternative mixture graph", showedges = T)
#plot.cg_(trand2,main="Alternative mixture graph 2", showedges = T)
plot.cg_(trand2,main="Alternative mixture graph 2", showedges = T)

plot.cg_(tg5)


## Perform inference
tinf1=infer_graph(trand2, data, maxiter =200, movefreqs = c(1/4,1/4,1/4,1/4), verbose = T, plotprogress = T, losstol = 1e-31)

trand1=mixedge_(tinf1$g,5,3,0.5,0.2)
trand2=mixedge_(trand1,3,1,0.5,0.2)

plot.cg_(trand3)


## Add to this graph a "mixture edge" from node 1 to node 3
# tg2=mixedge_(tg0, 1, 3, 0.5, 0.2)
# tg3 = mixedge_(tg1, 1, 3, 0.5, 0.4)
# 
# tg2=mixedge_(tg1,2,4,0.5,0.1) ## Secondly with weight 0.8
# tg3=mixedge_(tg2,1,4,0.5,0.2)
# 
# plot.cg_(tg0, main="Original graph",showedges = T )
# #plot.cg_(tg1,main="Mixture graph", showedges = T)
# #plot.cg_(tg2,main="Mixture graph 2", showedges = T)
# plot.cg_(tg2)


trand2=mixedge_(tinf2$g,3,1,0.5,0.2)

tinf3=infer_graph(trand2, data, maxiter = 200, movefreqs = c(0,1/2,0,1/2,0), verbose = T, plotprogress = T, losstol = 1e-31)

trand3=mixedge_(tinf3$g,1,2,0.5,0.2)
# transvector(data)*
plot.cg_(trand3)

tinf4=infer_graph(trand3, data, maxiter = 200, movefreqs = c(1/3,1/3,0,1/3,0), verbose = T, plotprogress = T, losstol = 1e-31)

ng = myregraftmixedge_(tinf4$g,i = 11,source = 10,target = 3,0.4,0.1)

tg = evaluate_proposal(ng,data)

plot.cg_(tg$proposal)

tg = evaluate_proposal(trand3,data)
plot.cg_(tg$proposal)

plot.cg_(tg3,main="Actual graph", showedges = T)
plot.cg_(tinf1$g, main="Predicted graph")

plot.cg_(tinf1$g, main = "Inferred graph",showedges = T)
plot.cg_(trand1, main="Inferred graph with admix",showedges = T )

tinf2=infer_graph(trand3, data, maxiter = 200, movefreqs = c(1/3,1/3,0,1/3),verbose = T, plotprogress = T, losstol = 1e-32)

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
tg0=simCoal_(9,labels=c("A","B","C","D","E","F","G","H","O"), outgroup="O")
## Add to this graph a "mixture edge" from node 1 to node 3
tg1=mixedge_(tg0, 1, 3, 0.5, 0.2)
tg2=mixedge_(tg1,6,4,0.5,0.4) ## Secondly with weight 0.4
tg3=mixedge_(tg2,12,10,0.5,0.3)

par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
#plot.cg_(tg1,main="Mixture graph", showedges = T)
plot.cg_(tg3,main="Mixture graph", showedges = T)

# simulate data
data=ccov_dag_(tg3)

# Make a random proposal graph
trand0=simCoal_(9,labels=c("A","B","C","D","E","F","G","H","O"),outgroup="O")
trand1 = mixedge_(trand0,2,3,0.5,0.5)
trand2=mixedge_(trand1,1,2,0.5,0.4) ## Secondly with weight 0.4
trand3=mixedge_(trand2,4,2,0.5,0.4) ## Secondly with weight 0.4

plot.cg_(trand0, main="Proposal graph",showedges = T )
#plot.cg_(trand1, main="Proposal graph 2",showedges = T )
plot.cg_(trand3, main="Proposal graph 2",showedges = T )

## Perform inference without mixture edges - just using prune/regraft
tinf1=infer_graph(trand0, data, movefreqs = c(1/3,1/3,0,1/3,0) , maxiter = 200 ,verbose =T, plotprogress = T)

plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(tg1, main="Actual graph",showedges = T )

## Add mixture edge to inferred graph
trand1 = mixedge_(tinf1$g,2,3,0.5,0.5)
plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(trand1, main="Inferred graph with admix",showedges = T )

## Perform inference again with mixture
tinf2=infer_graph(trand2, data, movefreqs = c(1/4,1/4,1/4,1/4) , maxiter = 200 ,verbose =T)

plot.cg_(tinf2$g, main = "Inferred graph", showedges = T)
plot.cg_(tg2, main="Actual graph",showedges = T )

plot.cg_(trand2)

##########################



# tg1=mixedge_(tg0, 1, 3, 0.5, 0.2)
# tg2=mixedge_(tg1,2,4,0.5,0.1) ## Secondly with weight 0.8
# tg3=mixedge_(tg2,1,4,0.5,0.2)
g = tinf1$g
plot.cg_(g)



g1=myregraftmixedge_(g,i = 10,source=2,target=4,0.5,0.1)
plot.cg_(g1)



g2=myregraftmixedge_(g1,i = 12,source=1,target=4,0.5,0.1)
plot.cg_(g2)


ng = evaluate_proposal(g,data, factr = 0.1,pgtol=0.1)
tg = ng$proposal
ng$loss
plot.cg_(tg)

g2=myregraftmixedge_(g1,i = 11,source=4,target=3,0.5,0.16)
g3=myregraftmixedge_(g2,i = 12,source=3,target=4,0.5,0.06)
plot.cg_(g3)

## "proposal: movetype  prune/regraft and mix regraft swap: 1 3 mixnode: 11 mixswap: 2 1"
plot.cg_(trand3)
g4 = mypruneandmixregraft(trand3,prsource = 1,prtarget = 3,rem = 11,mixsource = 2,mixtarget = 1)
plot.cg_(g4)
plot.cg_(evaluate_proposal(g4,data)$proposal)

plot.cg_(g4)


ret = randomprunemixpars(g)
print(ret)
g2 = mypruneandmixregraft(g, prsource = ret$prsource,prtarget = ret$prtarget,rem = ret$rem,mixsource = ret$mixsource,mixtarget =ret$mixtarget)
plot.cg_(g2)


g2 = mypruneandmixregraft(g1, prsource = 11,prtarget = 3,rem = 10,mixsource = 2,mixtarget = 1)
plot.cg_(g2)

g3 = mypruneregraft_(g, source=1,target=3)
plot.cg_(g3)


g = trand3
plot.cg_(g)

g1 = myregraftmixedge_(g,i = 18,source = 3,target = 7,alpha = 0.5,w = 0.2)

g2 = myregraftmixedge_(g1 ,i = 20,source = 2,target = 12,alpha = 0.5,w = 0.2)

plot.cg_(g2)

tinf4=infer_graph(g2, data, maxiter = 200, movefreqs = c(0,1/2,0,1/2,0), verbose = T, plotprogress = T, losstol = 1e-31)

# [1] "proposal: movetype mix regraft swap: 3 7 mixnode: 18"
# [1] "Iteration: 2  current loss: 2.13001152721922"
# [1] "proposal: movetype: reverse mixture swap: 7 6 mixnode: 20"
# [1] "proposal: movetype mix regraft swap: 7 12 mixnode: 18"
# [1] "proposal: movetype mix regraft swap: 2 12 mixnode: 20"
# [1] "Iteration: 5  current loss: 0.784875714378374"
# [1] "proposal: movetype: reverse mixture swap: 12 2 mixnode: 20"
# [1] "proposal: movetype mix regraft swap: 5 2 mixnode: 18"
# [1] "proposal: movetype mix regraft swap: 6 13 mixnode: 18"
# Error in rbind(deparse.level, ...) : 
#   numbers of columns of arguments do not match



####

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
# due to unmodelled migration elsewhere in the graph.




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

