
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
tinf1=infer_topo(trand1, data, maxiter = 50 ,verbose = T, patience = 15)
ccov_dag_(tinf1[[1]])
plot.cg_(tg1,main="Data graph", showedges = T)
plot.cg_(tinf1[[1]], main="Predicted graph")
tinf1[[3]]



g0 = tg1
plot.cg_(g0)
pars = get_depths(g0, data)
pars
g = parameterise_(g0,pars)
plot.cg_(g)






parameterise_<-function(g,pars){
  ## Take parameters and put them in their correct place in g
  
  nodes = c(g$tips, g$internal)
  drift = nodes[nodes != g$root]
 # mix = g$mix
  
 # np=npars_(g)
 # ng=length(pars)/np["tot"]
 # if(floor(ng)!=ng) stop("Invalid parameterisation")
  
  tpars = format_params(g,pars)
  tg = g
  for(j in 1:length(drift))  {
    tg$nl[[ drift[j] ]]$d = tpars[j]
    }
  
  #if(length(mix)>0) for(j in mix) tg$nl[[ tg$nl[[mix[j]]]$children[2] ]]$w = tpars[j+np["nd"]]
  return(tg)
}





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
