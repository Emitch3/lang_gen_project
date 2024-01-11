
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
tinf1=infer_topo(trand0, data, maxiter = 30 ,verbose = T)
ccov_dag_(tinf1[[1]])
plot.cg_(tinf1[[1]])
tinf1[[3]]



g0 = tg1
p = get_depths(g0, data)
g = parameterise_(g0,format_params(g0,p) , transform = F)







g0 <- simCoal(4,labels=c("A","B","C","O"), times = "coal",outgroup="O")
g1 <- mixedge(g0, 6, 2, 0.5, 0.2)


tg1 = 

plot.cg_(tg0)
plot.cg_(tg1)

data = ccov_dag_(tg1)
data

p = get_depths(tg1, data)
p
printDepths(tg1)

y = parameterise_(tg1, pars = p, transform = F)

ccov_dag_(y)
data

sum((ccov_dag_(y) - data)^2)

plot.cg_(tg1)
plot.cg_(y)


edges.cg_(tg0)

Wsvd = svd(W)

Wdiag = diag(1/Wsvd$d)

Wsvd$v %*% Wdiag %*% t(Wsvd$u) %*% diag(V)


tg0=simCoal_(4,labels=c("A","B","C","O"),outgroup="O")
## Add to this graph a "mixture edge" from node 2 to node 3
tg1=mixedge_(tg0,2,3,0.5,0.2) ## Firstly with weight 0.2
tg2=mixedge_(tg0,2,3,0.5,0.8) ## Secondly with weight 0.8
## Simulate covariances from these two models
tpred1=ccov_dag_(tg1)
tpred2=ccov_dag_(tg2)

## Make a random pair of alternative graphs
trand0=simCoal_(4,labels=c("A","B","C","O"),outgroup="O")
trand1=mixedge_(trand0,2,3,0.5,0.5)  ## Careful not to involve the outgroup

#tinf1=infer_dag(trand1,tpred1)

tg0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")
tg1=mixedge(tg0,2,3,0.5,0.2) 


plot.cg(tg1)
plot(x)

g = tg0

source = 2
target = 1


################################################################################


reversemixture<-function(g,mixrev){
  ####### INCOMPLETE
  if(is(g,"cglist")){
    ret=lapply(g,reversemixture,mixrev=mixrev)
    class(ret)="cglist"
    return(ret)
  }
  ## Swap the role of mixrev (a mixture node) and the other parent of the target
  target=g$nl[[mixrev]]$cr
  opar= g$nl[[target]]$pl
  if(opar==target) stop("reversemixture problem")
  ## Swap at the parents of the target
  isleftchild=(g$nl[[opar]]$cl==target)
  if(isleftchild){
    ochild=g$nl[[opar]]$cr
  }else{
    ochild=g$nl[[opar]]$cl
  }
  ## Find out the parent of the other child
  paropar=g$nl[[opar]]$pl
  ## Find out the parent of the mixrev
  parmix=g$nl[[mixrev]]$pl
  g$nl[[ochild]]$pl=paropar
  parisleftchild=(g$nl[[paropar]]$cl==opar)
  mixisleftchild=(g$nl[[parmix]]$cl==mixrev)
  ## Reverse the children of the parents
  if(parisleftchild){
    g$nl[[paropar]]$cl=mixrev
  }else{
    g$nl[[paropar]]$cr=mixrev
  }
  g$nl[[mixrev]]$pl=paropar
  if(mixisleftchild){
    g$nl[[parmix]]$cl=opar
  }else{
    g$nl[[parmix]]$cr=opar
  }
  g$nl[[opar]]$pl=paropar
  
  ## tpl=g$nl[[target]]$pl
  ## tpr=g$nl[[target]]$pr
  ## ## Swap at the target
  ## g$nl[[target]]$pl=tpr
  ## g$nl[[target]]$pr=tpl
  g
}



g0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")

g1=mixedge(g0,2,3,0.5,0.2) ## Firstly with weight 0.2
plot(g1)

g2 = reversemixture(g1, mixrev = 8)
plot(g2)


################################################################################


tg0 = simCoal_(4,labels=c("A","B","C","O"), alpha = 0.75,outgroup="O")
plot.cg_(tg0)

t1 = mypruneregraft_(tg0, source = 1, target = 3)
plot.cg_(t1)


g0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")
plot(g0)

t = mypruneregraft(g0, source = 1, target = 3)
plot(t)


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
