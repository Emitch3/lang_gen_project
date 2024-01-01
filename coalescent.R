
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

tg0 = simCoal_(6,labels=c("A","B","C", "D", "E","O"), alpha = 0.75,outgroup="O")
plot.cg_(tg0)

tg1 <- mixedge_(tg0, 3, 2, 0.5, 0.2)
tg2 <- mixedge_(tg1, 1, 2, 0.5, 0.2)
#tg3 <- removemixedge_(tg2, 8)
plot.cg_(tg0)
#plot.cg_(tg3)

edges.cg_(tg1)
edges.cg_(tg3)

#tg2<- mixedge_(tg1, 1, 5, 0.5, 0.3)

par(mfrow=c(1,2))
plot.cg_(tg0)
plot.cg_(tg1)
plot.cg_(tg2)

par(mfrow=c(1,2))
plot.cg_(tg0, main="Original graph",showedges = T )
plot.cg_(tg1,main="Mixture graph", showedges = T)
plot.cg_(tg2, main = "Mixture graph 2")

g0 <- simCoal(4,labels=c("A","B","C","O"), times = "coal",outgroup="O")
g1 <- mixedge(g0, 6, 2, 0.5, 0.2)

ccov_dag_(tg0)
ccov_dag_(tg2)

ccov_dag_(g0,old = T)

ccov_dag(g1)
ccov_dag_(g1,old = T)

plot(g1)

c = ccov_dag_(tg0)
ccov_dag_(tg1)
ccov_dag_(tg2)

c_Center(c,d)

ctree_loss2_(tg1,c)

V = ccov_dag_(tg0)

W = get_weightmatrix(tg0)

sol = find_lengths(tg0)


find_lengths(tg2)

system.time(find_lengths(tg1))


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




pruneregraft_ = function(g,source,target,careful=FALSE){
  warning("May need debugging")
  ## source = node to be moved
  
  ## prunes source node branch and attaches it to the target node branch
  
  ## Run this ONLY on type=="split" nodes!
  gstart=g
  
  psourceroot=FALSE
  ptargetroot=FALSE

  ptarget=g$nl[[target]]$parent[1] # non-admix parent of target
  
  psource=g$nl[[source]]$parent[1] # non-admix parent of source. (This node may be removed)
  
  ocsource = g$nl[[psource]]$children 
  ocsource = ocsource[ocsource!=source] # siblings of source node
  
  ppsource=g$nl[[psource]]$parents # grandparent of source (will become parent to other child of this parent)
  
  if(ocsource==target) { # is target sibling of source node?
    print(paste("source",source,"and target",target,"have same parent",psource))
    return(g)
  }
  
  if(is.na(ppsource) ) {
    print(paste("NOTE: parent of source is the root ( node",psource,")"))
    psourceroot=TRUE
    ## psource will become internal
    g$nl[[psource]]$w = 1
    ## ocsource will become the root
    g$nl[[ocsource]]$d = NA
    g$nl[[ocsource]]$parents = NA
    g$root=ocsource
  }else if(length(ocsource) == 1) {
    ## Remove psource since can't have split node with only 1 child
    g$nl[[ocsource]]$d = g$nl[[ocsource]]$d + g$nl[[psource]]$d
    g$nl[[ocsource]]$parents = ppsource
    
    g$nl[[ppsource]]$children = c(g$nl[[ppsource]]$children[g$nl[[ppsource]]$children != psource],ocsource )
    # update graph information
    g$nl[[psource]] =cnode_(id = psource,children = NA, parents = NA, w=1,d=0,type="spare") 
    g$mix=g$mix[-which(g$mix==psource)]
    g$internal=g$internal[-which(g$internal==psource)]
    g$spare=c(psource,g$spare)
  }else{
    g$nl[[psource]]$children = ocsource}  
  
  ## Update the parent of the target
  if(is.na(ptarget)) print(paste("parent of target is the root ( node",ptarget),")")
  
  g$nl[[ptarget]]$children = c(g$nl[[ptarget]]$children, source)
  
  # targetisleftchild = (g$nl[[ptarget]]$cl==target)
  # if(targetisleftchild){
  #   g$nl[[ptarget]]$cl=psource
  # }else{
  #   g$nl[[ptarget]]$cr=psource
  # }
  
  ## Update the split node, the old parent node
  #g$nl[[psource]]$d=g$nl[[target]]$d/2
  #g$nl[[psource]]$cl=target
  #g$nl[[psource]]$cr=source
  #g$nl[[psource]]$pl=ptarget
  ## Update the target and source
  
  #g$nl[[target]]$d=g$nl[[target]]$d/2
  g$nl[[target]]$parents = c( g$nl[[target]]$parents, psource)
  g$nl[[source]]$parents=c(g$nl[[source]]$parents, psource )
  
  ## Done!
  if(careful) {
    if(!checkgraph(g)){
      print(paste("Moving",source,"to",target))
      print(g)
      stop("Created invalid graph!")
    }
  }
  g
}



mypruneregraft_ = function(g,source,target,careful=FALSE){
  ## prunes source node branch and attaches it to the target node branch
  
  ## Run this ONLY on type=="split" nodes!
  gstart=g
  psourceroot=FALSE
  ptargetroot=FALSE
  
  ## Parent of source needs removing
  ptarget=g$nl[[target]]$parent[1] # non-admix parent of target
  psource=g$nl[[source]]$parent[1] # non-admix parent of source. This node is to be moved to the target branch
  
  ocsource = g$nl[[psource]]$children 
  ocsource = ocsource[ocsource!=source] # siblings of source node
  
  ppsource=g$nl[[psource]]$parents # grandparent of source, will become parent to other siblings of source node

  if(ocsource==target) {
    print(paste("source",source,"and target",target,"have same parent",psource))
    return(g)
  }
  if(is.na(ppsource) ) {
    print(paste("NOTE: parent of source is the root ( node",psource,")"))
    psourceroot=TRUE
    ## psource will become internal
    g$nl[[psource]]$w = 1
    ## ocsource will become the root
    g$nl[[ocsource]]$d = NA
    g$nl[[ocsource]]$parents[1] = NA
    g$root=ocsource
  }else{
    ## Update the other child of the parent of the source
    g$nl[[ocsource]]$d = g$nl[[ocsource]]$d + g$nl[[psource]]$d
    g$nl[[ocsource]]$parent[1] = ppsource
    g$nl[[ppsource]]$children = c(g$nl[[ppsource]]$children[g$nl[[ppsource]]$children != psource],ocsource)
  }

  # ## Update the parent of the target
  if(is.na(ptarget)) print(paste("parent of target is the root ( node",ptarget),")")
  
  g$nl[[ptarget]]$children = c(g$nl[[ptarget]]$children[g$nl[[ptarget]]$children != target], psource) 
  
  ## Update the split node, the old parent node
  g$nl[[psource]]$d=g$nl[[target]]$d/2
  g$nl[[psource]]$children = c(target, source)
  g$nl[[psource]]$parents[1] = ptarget
  
  ## Update the target and source
  g$nl[[target]]$d=g$nl[[target]]$d/2
  g$nl[[target]]$parents[1] = psource    
  g$nl[[source]]$parents[1] = psource
  ## Done!
  if(careful) {
    if(!checkgraph(g)){
      print(paste("Moving",source,"to",target))
      print(g)
      stop("Created invalid graph!")
    }
  }
  g
}

tg0 = simCoal_(4,labels=c("A","B","C","O"), alpha = 0.75,outgroup="O")
plot.cg_(tg0)

t1 = mypruneregraft_(tg0, source = 1, target = 3)
plot.cg_(t1)


g0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")
plot(g0)

t = mypruneregraft(g0, source = 1, target = 3)
plot(t)

#pruneregraft_(tg1, source = )


regraftmixedge_<-function(g,i,source,target,alpha,w){
  if(is(g,"cglist")){
    ret=lapply(g,function(x){
      regraftmixedge_(x,i,source,target,alpha,w)
    })
    class(ret)="cglist"
    return(ret)
  }
  ## Does a complete removal and replacement of a mixture edge, with specified parameters
  if(!is.na(i)) g=removemixedge_(g,i)
  g=mixedge_(g,source,target,alpha,w)
  g
}



XX = pruneregraft(tg0, source ,target, careful = F)


edges.cg_(g)
edges.cg_(XX)
plot.cg_(g)
plot.cg_(XX)


g0=simCoal(4,labels=c("A","B","C","O"),outgroup="O")
g1=mixedge(g0,1,3,0.5,0.2) 

gx = mypruneregraft(g = g1, source = 3, target = 1)

plot(g1)
plot(gx)


















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
