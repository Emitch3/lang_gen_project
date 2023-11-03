

setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project")
source("mixfunctions.R")
source("extendedfunctions.R")
library("Clarity")


ccov_dag_<-function(g){
  ## Enumerate all trees implied by a DAG and compute the covariance matrix for each.
  if(is(g,"cglist")){
    ret = lapply(g,ccov_dag_)
    return(ret)
  }
  if(length(g$mix)==0){
    csum=ccov_tree_(g)
    rownames(csum)=colnames(csum)=g$tip.label
    return(csum)
  }
  alltrees=cenumerate_trees(g) # tree ids
  allmatrices=cenumerate_weightmatrices(g,alltrees)
  print("w")
  allcovs=cenumerate_covmatrices(g,alltrees)
  ## allw=sapply(g$mix,cweight,g=g) # weights
  clist=lapply(1:dim(alltrees)[1],function(i){
    c=allcovs[[i]]
    cmult=allmatrices[[i]]
    c * cmult 
  })
  csum=Reduce('+',clist)
  rownames(csum)=colnames(csum)=g$tip.label
  csum
}


cenumerate_trees<-function(g){
  ## Enumerate all trees that can be induced by different switches of mixture edges
  allvals=lapply(g$mix,function(x){c(0,1)})
  names(allvals)=g$mix
  expand.grid(allvals)
}



cenumerate_weights<-function(g,alltrees=NULL){
  ## Enumerate all tree weighting matrices and return them as a list
  if(all(is.null(alltrees))) alltrees = cenumerate_trees(g)
  sapply(1:dim(alltrees)[1],function(i)
    cmultweight(g,alltrees[i,,drop=FALSE])
  )
}


cmultmatrix_ <- function(g,mixsetting){
  ## Takes a mixture DAG g
  ## And a configuration of mixture on/off settings, mixsetting
  ## which is a row vector of length #mixtures with column names giving the node labels of the mixture edges
  ## with values 0 or 1 for whether we follow them or not
  ## We then return the mixture weighting matrix
  cmult=matrix(1,g$n,g$n)
  allw=sapply(g$mix,cweight,g=g) # weights
  for(j in 1:length(allw)){
    ttips=tipsunder(g,g$nl[[g$mix[j] ]]$cr) ## tips under the mixture target
    tw=allw[j]
    if(mixsetting[j]==0) tw = 1 - tw
    cmult[ttips,ttips]= cmult[ttips,ttips] * tw
    cmult= cmult * tw
  }
  cmult 
}

alltrees=cenumerate_trees(tg1)
alltrees

cenumerate_weightmatrices(tg1,alltrees)

tpred1=ccov_dag_(tg1)

tg0=simCoal_(4,labels=c("A","B","C","O"),outgroup="O")

## Add to this graph a "mixture edge" from node 2 to node 3
tg1=mixedge_(tg0,2,3,0.5,0.2) ## Firstly with weight 0.2
tg2=mixedge_(tg0,2,3,0.5,0.8) ## Seconldly with weight 0.8
## Simulate covariances from these two models


ccov_dag<-function(g){
  ## Enumerate all trees implied by a DAG and compute the covariance matrix for each.
  if(is(g,"cglist")){
    ret = lapply(g,ccov_dag)
    return(ret)
  }
  if(length(g$mix)==0){
    csum=ccov_tree(g)
    rownames(csum)=colnames(csum)=g$tip.label
    return(csum)
  }
  alltrees=cenumerate_trees(g) # tree ids
  allmatrices=cenumerate_weightmatrices(g,alltrees)
  allcovs=cenumerate_covmatrices(g,alltrees)
  ## allw=sapply(g$mix,cweight,g=g) # weights
  clist=lapply(1:dim(alltrees)[1],function(i){
    c=allcovs[[i]]
    cmult=allmatrices[[i]]
    c * cmult 
  })
  csum=Reduce('+',clist)
  rownames(csum)=colnames(csum)=g$tip.label
  csum
}

g2 = mixedge(g1,1,2,0.5,0.2) ## Firstly with weight 0.2

cenumerate_trees<-function(g){
  ## Enumerate all trees that can be induced by different switches of mixture edges
  allvals=lapply(g$mix,function(x){c(0,1)})
  names(allvals)=g$mix
  expand.grid(allvals)
}

allvals=lapply(g2$mix,function(x){c(0,1)})
allvals
names(allvals)=g2$mix
cenumerate_trees(g2)
g2$mix


cenumerate_weightmatrices(g2)

tg0 = simCoal_(4,labels=c("A","B","C","O"),outgroup="O")

g0 = simCoal(4,labels=c("A","B","C","O"),outgroup="O")
g1 = mixedge(g0,2,3,0.5,0.2) ## Firstly with weight 0.2
plot(g2)



cenumerate_weightmatrices<-function(g,alltrees=NULL){
  ## Enumerate all tree weighting matrices and return them as a list
  if(all(is.null(alltrees))) alltrees = cenumerate_trees(g)
  lapply(1:dim(alltrees)[1],function(i)
    cmultmatrix(g,alltrees[i,,drop=FALSE])
  )
}


alltrees = cenumerate_trees(g)
alltrees[1,,drop=F]


cmultmatrix<-function(g,mixsetting){
  ## Takes a mixture DAG g
  ## And a configuration of mixture on/off settings, mixsetting
  ## which is a row vector of length #mixtures with column names giving the node labels of the mixture edges
  ## with values 0 or 1 for whether we follow them or not
  ## We then return the mixture weighting matrix
  cmult=matrix(1,g$n,g$n)
  allw=sapply(g$mix,cweight,g=g) # weights
  for(j in 1:length(allw)){
    ttips=tipsunder(g,g$nl[[g$mix[j] ]]$cr) ## tips under the mixture target
    tw=allw[j]
    if(mixsetting[j]==0) tw = 1 - tw
    cmult[ttips,ttips]= cmult[ttips,ttips] * tw
    cmult= cmult * tw
  }
  cmult 
}

mixsetting = alltrees[1,,drop=F]
cmultmatrix(g2, alltrees[1,,drop=F])

cmult=matrix(1,g2$n,g2$n)

allw=sapply(g$mix,cweight,g = g2)
allw

ttips=tipsunder(g2,g2$nl[[g2$mix[1] ]]$cr)

tw = allw[1]

mixsetting[1] == 0

tw = 1 - tw

cmult[ttips,ttips] = tw

alltrees
mixsetting

cmultweight(g, mixsetting)












ccov_dag(g)

printWeight(g2)


cenumerate_weights(g)

cmultweight<-function(g,mixsetting){
  ## Takes a mixture DAG g
  ## And a configuration of mixture on/off settings, mixsetting
  ## which is a row vector of length #mixtures with column names giving the node labels of the mixture edges
  ## with values 0 or 1 for whether we follow them or not
  ## We then return the weight
  cmult=1
  allw=sapply(g$mix,cweight,g=g) # weights
  for(j in 1:length(allw)){
    tw=allw[j]
    if(mixsetting[j]==0) tw = 1 - tw
    cmult= cmult * tw
  }
  cmult 
}

cmultweight(g,alltrees[1])



tipsunder_(tg0,7)


















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



paths = function(edges, i){
  idx = which(edges$child == i)
  path = NULL
  edges$new_weight = 0
  dict_weights = list()
  dict_weights[[as.character(i)]] = 1
  
  while (length(idx) > 0) {
    E = edges[idx,]

    edge = paste0(E$child,"-", E$parent)
    
    for(k in 1:length(edge)) {
      E[k,]$new_weight = E[k,]$weight*dict_weights[[as.character(E[k,]$child)]]
      dict_weights[[as.character(E[k,]$parent)]] = E[k,]$new_weight
    }    

    r = cbind(edge, E)
    path = rbind(path,r)

    idx = which(edges$child %in% E$parent)
  }
  #path <- aggregate(new_weight ~ edge + parent + child + length + mix, data = path, FUN = sum)
  paths = split(path, path$new_weight)
  return(paths)
}

O = function(p_i, p_j){
  overlap_i = p_i[p_i$edge %in% p_j$edge,]
  overlap_j = p_j[p_j$edge %in% p_i$edge,]
  wi = overlap_i$weight
  wj = overlap_j$weight
  ci = overlap_i$length
  res = wi*wj*ci
  return(sum(res))
}

Vij = function(Paths_i, Paths_j){
  v = sum(sapply(Paths_i, function(path_i) sum(sapply(Paths_j, function(path_j) O(path_i, path_j)))))
  return(v)
}

cov = function(g) {
  ntips = g$n
  cov = matrix(0, ntips, ntips)
  edges = edges.cg_(g)
  
  cov = sapply(1:ntips, function(i) {
    sapply(1:ntips, function(j) {
      Vij(paths(edges, i), paths(edges, j))})
  })
  dimnames(cov) = list(g$tip.label, g$tip.label)
  return(cov)
}


cov(tg0)

ccov_dag(g)

system.time(cov(g))
system.time(ccov_dag(g))

ntips = g$n
cov = matrix(0, ntips, ntips,dimnames = list(g$tip.label,g$tip.label))
edges = edges.cg(g)
for(i in 1:ntips){
 for(j in 1:ntips){
   cov[i,j]  = Vij(paths(edges,i),paths(edges,j))
 }
}

cov

ccov_dag(g)


# result = 0
# for(i in 1:length(Paths_i)){
#   for(j in 1:length(Paths_j)){
#     result = result + O(Paths_i[[i]],Paths_j[[j]])
#   }
# }
# result









O(Paths_i[[1]],Paths_i[[1]])


p_i[p_i$edge %in% p_j$edge,]

O = function(p_i, p_j){
  overlap = p_i$edge[p_i$edge %in% p_j$edge]
  w_i = p_i[which(p_i$edge %in% overlap),]$new_weight
  c_i = p_i[which(p_i$edge %in% overlap),]$length
  w_j = p_j[which(p_j$edge %in% overlap),]$new_weight
  return(sum(w_i*w_j*c_i))
}

overlap








d = list()
d[[as.character(E$parent)]] <- E$weight
d

d[["8"]]

O()

#edges
#i = 3
#j = 3

p_i = paths_i(edges,3)
p_j = paths_i(edges,3)
p_i
p_j

overlap = p_i[]$edge[p_i$edge %in% p_j$edge]




split(p_j, p_j$new_weight)

#ccov(p_i, p_j)



O(p_i, p_j) 

#Vij = function
#  O(p_i, p_j) + o


overlap$weight
p_j$weight






overlap = p_j[p_j$edge %in% p_i$edge,]
overlap

aggregated = aggregate(new_weight ~ edge + parent + child + length + mix, data = overlap, FUN = prod)
aggregated
sum(aggregated$new_weight^2*aggregated$length)

ccov_dag(g)[3,3]

overlap_edges <- p_j[ p_j %in% p_i] #names(table(df$edge))[table(df$edge) > 1]
overlap_edges

subset_df <- df[df$edge %in% overlap_edges, ]
subset_df

result <- aggregate(new_weight ~ edge + length, data = subset_df, FUN = function(x) prod(x))

result

sum(result$length*result$new_weight)

overlapW = W[p_j$edge %in% p_i$edge]
overlapLength = p_j$length[p_j$edge %in% p_i$edge]

#overlap = p_j[p_j$edge %in% p_i$edge,]
overlapLength

sum(overlapW*overlapLength)

#sum(overlap$new_weight*overlap$length)

ccov_dag(g)[1,1]




