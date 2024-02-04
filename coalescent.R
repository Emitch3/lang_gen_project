
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

############################################################################################
### FIX NN_INTERCHANGE BUG: SWAP NODE WITH ADMIX ON IT, NEEDS TO TAKE ADMIX NODE WITH IT ###
############################################################################################



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
tg2=mixedge_(tg1,2,3,0.5,0.8) ## Secondly with weight 0.8

par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)
#plot.cg_(tg2,main="Mixture graph 2", showedges = T)

# simulate data
data=ccov_dag_(tg1)

## Make a random pair of alternative graphs
trand0=simCoal_(5,labels=c("A","B","C","D","O"),outgroup="O")
trand1=mixedge_(trand0,2,3,0.5,0.2)  ## Careful not to involve the outgroup
#trand2=mixedge_(trand1,1,3,0.5,0.2)

plot.cg_(trand0, main="Alternative graph",showedges = T )
plot.cg_(trand1,main="Alternative mixture graph", showedges = T)
#plot.cg_(trand2,main="Alternative mixture graph 2", showedges = T)

## Perform inference
tinf1=infer_graph(trand1, data, maxiter = 50, verbose = T, losstol = 0.0001)
#tinf1=infer_graph(tinf1[[1]], data, maxiter = 100 ,verbose = T, losstol = 0.001)

plot.cg_(tg1,main="Data graph", showedges = T)
plot.cg_(tinf1$g, main="Predicted graph")

g = tg0
trans = transvector(data)

#tg = nn_interchange(tinf1$g,source = 3,target = 4)
tg = mixedge_(g, source = 1, target = 3, alpha =0.3, w =0.5)

plot.cg_(tg)
#tg = myregraftmixedge_(tinf1$g, 10, source=1, target = 3,alpha = 0.5, w = 0.2 )
d = get_depths(tg, data, trans)
ng = parameterise_(tg, pars= d, what = "d")
plot.cg_(ng)

infer_mixparams(ng, data)$par

rg = infer_mixparams(ng, data)$g
plot.cg_(rg)

rg = parameterise_mix(ng, alphas = 0.5, weights = 0.2)
rg

ctree_loss2_(rg, data)
ctree_loss2_(tg1,data)
#data
#ccov_dag_(tinf1[[1]])
#losses = tinf1[[3]]
#losses
#min(losses)

# Covariance matrices
tinfpred1=ccov_dag_(tinf1[[1]])
tpred1=ccov_dag_(tg1)

## Get the plot that we would make, if we made a plot
pt=plot.cg_(tg1,show=FALSE)

Clarity_Chart(ordermatrix(tpred1,pt$order),scalefun=myscalefun2,text=T )
Clarity_Chart(ordermatrix(tinfpred1,pt$order),scalefun=myscalefun2,text=T)



g = trand1
plot.cg_(trand1)
p = get_depths(g, data, 1)
g = parameterise_(g, pars = p, what = "d")
plot.cg_(g)

g = removemixedge_(trand1,i = 10)
p = get_depths(g, data, 1)
g = parameterise_(g, pars = p, what = "d")
plot.cg_(g)


##################################

## Simulate a 4 population model where there is an outgroup
tg0=simCoal_(9,labels=c("A","B","C","D","E","F","G","H","O"), alpha = 1,outgroup="O")
## Add to this graph a "mixture edge" from node 1 to node 3
tg1=mixedge_(tg0, 1, 3, 0.5, 0.2)
#tg2=mixedge_(tg1,4,3,0.5,0.8) ## Seconldly with weight 0.8


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

## Perform inference
tinf1=infer_graph(trand0, data, maxiter = 100 ,verbose = T, losstol = 0.0001)

#trand1 = mixedge_(tinf1$g,2,3,0,0.5)
#tinf2=infer_graph(tinf1$g, data, maxiter = 100 ,verbose = T, movefreqs = c(0,2/3,0,1/3,0),losstol = 0.01)
#tinf1=infer_graph(tinf1$g, data, maxiter = 100 ,verbose = T, losstol = 0.01)


plot.cg_(tinf1$g, main = "Inferred graph", showedges = T)
plot.cg_(tg1, main="Original graph",showedges = T )

#tinf2=infer_graph_hc(trand1, data, maxiter = 50 ,verbose = T, losstol = 0.01) 

plot.cg_(tinf2$g, main = "Predicted graph", showedges = T)
plot.cg_(tg1,main="Data graph", showedges = T)


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

