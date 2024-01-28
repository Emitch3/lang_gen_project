setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project")
source("mixfunctions.R")
library("Clarity")



cnode_<-function(id,parents = NA,children = NA,d=NA,t=NA,p=NA,w=NA,type=NA,mixchild=NA){
  ## This is the generic function to create or update a node
  ## Pass an id to create a new node
  ## Pass a cnode to update it
  ## Returns a cnode
  if(is(id,"cnode")){
    if(!is.na(parents)) id$parents = parents
    if(!is.na(children)) id$children = children
    if(!is.na(d)) id$d=d  # distance from parent node
    if(!is.na(t)) id$t=t # time from tips
    if(!is.na(p)) id$p=p # positions
    if(!is.na(w)) id$w=w 
    if(!is.na(type)) id$type=type
    return(id)
  }else{
    if(is.na(type)) type="split"
    r=list(id=id,
           parents = parents,
           children = children,
           d=d,
           t=t,
           p=p,
           w=w,
           type=type,
           mixchild=mixchild)
    class(r)="cnode"
    return(r)
  }
}

simCoal_<-function(n,times="coal", alpha = 1, labels=paste0("t",1:n),outgroup=numeric()){
  ## Simulate a coalescent tree 
  nodes = list()
  for(i in 1:n) nodes[[i]]=cnode_(i,t=0,d=0,w=1)
  alive = 1:n
  outgroupnums = alive[labels%in%outgroup]
  alive = alive[!labels%in%outgroup]
  time = 0
  fin = FALSE
  for(i in n:2){
    ii <- length(nodes) + 1    
    if(length(alive) < 3) nChildren = 2
      else  nChildren <-  sample(c(2,3), 1, prob = c(alpha, 1 - alpha)) #max(2,sample(1:length(alive), 1))
    clist <- sample(alive,nChildren)
    
    if(times == "coal"){
      t=rexp(1,i)
    }else if(times=="test"){
      t=(n-i+1)  
    }else{
      stop("Invalid times argument!")
    }        
    time=time+t
    if((length(alive)<2)&&(length(outgroup)>0)){ ## Add in the outgroups last
      alive=c(alive,outgroupnums)
      nChildren = 2
      clist <- sample(alive,2)
      fin = TRUE
    }
    
    nodes[[ii]] <- cnode_(ii, t=time, d=0)
    
    for(c in 1:nChildren){
      nodes[[ii]]$children[c] <- clist[c]
      nodes[[clist[c] ]] <- cnode_(nodes[[clist[c] ]], parents = ii, d=time-nodes[[clist[c] ]]$t ,w=1)
      nodes[[clist[c] ]]$parents <- ii
    }
    if(fin == TRUE) {break}
    alive=c(alive[!alive%in%clist],ii)
  }
  g=list(nl=nodes,
         tips=1:n,
         root=length(nodes),
         tip.label=labels,
         internal=(n+1):length(nodes),
         mix=numeric(),
         spare=numeric(),
         mixparmap=numeric(),
         outgroup=outgroup,
         outmix=outgroup,
         n=n)
  class(g)="cg"
  g
}



mixedge_<-function(g,source,target,alpha,weight){
  ## Add a mixture edge to a graph 
  if(g$nl[[target]]$type == "mixture") {stop("ERROR: target cannot be a mixture node")}
  
  # forbid sibling admixture
  #if(g$nl[[target]]$parents[1] == g$nl[[source]]$parents[1]) {stop("ERROR: cannot have admixture between siblings")}
  
  if(length(g$spare)>0){
    ti=tail(g$spare,1)
    g$spare=g$spare[-(g$spare%in%ti)]
  }else{
    ti=length(g$nl)+1 # new node id
  }
  ##    if(target%in%tipsunder(source)) stop("Invalid mixture edge?")
  parsource <- g$nl[[source]]$parents[1] # parent of source node
  ##    if(is.na(parsource)) stop("ERROR: tried to make the root a mixture!")
  ## Add the new node
  g$nl[[ti]] <- cnode_(ti,
                       parents = parsource,
                       children = c(source, target),
                       d=g$nl[[source]]$d*(1-alpha),
                       t=(g$nl[[parsource]]$t - g$nl[[source]]$t)*alpha +  g$nl[[source]]$t,
                       w=1,
                       type="mixture",
                       mixchild = target)

  
  ## Update the target
  l <- length(g$nl[[target]]$parents)
  
  g$nl[[target]]=cnode(g$nl[[target]])
  g$nl[[target]]$w[1] = g$nl[[target]]$w[1] - weight
  g$nl[[target]]$w[length(g$nl[[target]]$w) + 1] = weight

  
  g$nl[[target]]$parents[l+1] = ti
  ## Update the source node (called source)
  g$nl[[source]]=cnode(g$nl[[source]],d=g$nl[[source]]$d*alpha)
  g$nl[[source]]$parents[1] = ti
  ## Update the original source's parent node (called parsource)
  if(!is.na(parsource)) { ## It was not the root.
    idx <- which(g$nl[[parsource]]$children == source)
    g$nl[[parsource]] = cnode(g$nl[[parsource]])
    g$nl[[parsource]]$children[idx]=ti
  }else{ ## It was not the root
    g$root=ti
    g$nl[[source]]=cnode(g$nl[[source]],d=0.5,w=1)
  }
  ## Update the list of nodes
  g$internal=c(g$internal,ti)
  g$mix=c(g$mix,ti)
  g
}


removemixedge_ <- function(g,i,careful=TRUE){
 # warning("check weights")
  ## Remove the node at index i
  ## If we are careful, we expect a proper graph
  ## Otherwise we accept missing right children for "dangling" mixture nodes with no right child
  if(g$nl[[i]]$type!="mixture") stop(paste("ERROR: Invalid request to remove node",i,"which is not a mixture edge"))
  
  csource = g$nl[[i]]$children[1] # child from source node
  ctarget = g$nl[[i]]$children[2] # child from target node
  originalparent = g$nl[[i]]$parents[1] # mixture node's (split) parent
  if(is.na(csource) && is.na(ctarget)) stop(paste("ERROR: node",i,"has no children!")) 
  if(is.na(originalparent) && (i != g$root)) stop(paste("ERROR: node",i,"has no parents!"))
  #if(is.na(cr)&&careful) stop(paste("ERROR: node",i,"has no right child!"))  
  
  if(length(g$nl[[i]]$parents) > 1 ){ #if mixture node has any mixture parents
    print("Warning: admix to admix")
    mixparents = g$nl[[i]]$parents[2:length(g$nl[[i]]$parents)]
    for(p in mixparents){removemixedge_(g,p,careful=TRUE)}
  }
  
  ## Update source child
  g$nl[[csource]]$parents = originalparent
  g$nl[[csource]]$d  = g$nl[[i]]$d + g$nl[[csource]]$d
  
  ## Update target child
  pt <- g$nl[[ctarget]]$parents
  g$nl[[ctarget]]$parents <- pt[pt != i]
  g$nl[[ctarget]]$w = 1
  
  ## Update original parent
  # if(g$root!=i){
  idx <- which(g$nl[[originalparent]]$children == i)
  g$nl[[originalparent]]$children[idx] = csource
  
  # }else{ # Its the root; have to make the left child the new root
  #  g$root=cl
  #}
  
  ## update graph information
  g$nl[[i]] =cnode_(id= i,children = NA, parents = NA, w=1,d=0,type="spare") #id = cr
  g$mix=g$mix[-which(g$mix==i)]
  g$internal=g$internal[-which(g$internal==i)]
  g$spare=c(i,g$spare)
  g
}



tipsunder_<-function(g,i,w=FALSE){
  ## Return the set of tips underneath a specified node
  ## if w (weighted) return this as a matrix of the weight for tip t (rows) under node i
  ## if !w (unweighted) return this as a list of tip ids
  if(w){
    stop("Weighted version not implemented")
    ## Weighted version
    r=rep(0,g$n)
    if(!is.na(g$nl[[i]]$cl)) {
      if(g$nl[[g$nl[[i]]$cl]]$pl==i){ # Is this the left parent?
        w = 1 - g$nl[[g$nl[[i]]$cl]]$w
      }else w = g$nl[[g$nl[[i]]$cl]]$w
      print(paste("node",i,"left child",g$nl[[i]]$cl,"weight",w))
      r = r + w * tipsunder(g,g$nl[[i]]$cl,w=w)
    }
    if(!is.na(g$nl[[i]]$cr)){
      if(g$nl[[g$nl[[i]]$cr]]$pl==i){ # Is this the left parent?
        w = 1 - g$nl[[g$nl[[i]]$cr]]$w
      }else w= g$nl[[g$nl[[i]]$cr]]$w
      print(paste("node",i,"right child",g$nl[[i]]$cr,"weight",w))
      r=r + w * tipsunder(g,g$nl[[i]]$cr,w=w)
    }
    if((is.na(g$nl[[i]]$cl))&(is.na(g$nl[[i]]$cr)) & (i%in%g$tips) ) {
      r[i]=r[i] + 1
    }
    r
  }else{
    ## Unweighted version
    r=c()
    children <- g$nl[[i]]$children
    if(!all(is.na(children))){
      for(k in 1:length(children)){
        r=c(r,tipsunder_(g,children[k],w=w))}
    }
    if((all(is.na(children))) & (i%in%g$tips)) r=c(r,i)
    r
  }
}


ccov_tree_<-function(g){
  ## Compute the induced covariance for a tree-like g
  ## Weightings are ignored
  c=matrix(0,ncol=g$n,nrow=g$n)
  for(i in g$internal){
    if(!is.na(g$nl[[i]]$d)){
      tu=tipsunder_(g,i)
      c[tu,tu]=c[tu,tu]+g$nl[[i]]$d
    }
  }
  for(i in g$tips){
    c[i,i]=c[i,i]+g$nl[[i]]$d
  }
  c
}



assignheights <-function(g,rightalign=FALSE){
  nleft=0
  nnodes=length(g$nl)
  i=g$root
  node=g$root
  for(i in 1:length(g$nl)) g$nl[[i]]$tmp=0 # tracks completed
  
  while(TRUE){
    ## Situations:
    ## 1. We have a left and need to go there
    lefttodo = ((!is.na(g$nl[[node]]$children[1])) && (g$nl[[g$nl[[node]]$children[1]]]$tmp==0))
    ## 2. We place ourself
    selftodo = (g$nl[[node]]$tmp==0)
    ## 3. Right to do
    righttodo = ((!is.na(g$nl[[node]]$children[2])) &&
                   (g$nl[[g$nl[[node]]$children[2]]]$tmp==0) &&
                   (g$nl[[node]]$type!="mixture") )
    ## 4. Return to parent if we have one
    atroot = (g$root==node)
    ## 5. Stop if -back- at root
    if(lefttodo){
      node=g$nl[[node]]$children[1]
      next;
    }else if(selftodo) {
      g$nl[[node]]$tmp=1
      g$nl[[node]]$p=nleft #/nnodes
      if(is.na(g$nl[[node]]$parents[1])){ nleft=nleft + 1
      }else if(g$nl[[ g$nl[[node]]$parents[1] ]]$type=="split") {  nleft=nleft + 1
      }
      next;
    }else if(righttodo){
      node=g$nl[[node]]$children[2]
      next;            
    }else if(!atroot){
      node=g$nl[[node]]$parents[1]
      next;
    }else{
      break;
    }
  }
  g
}

assignlocation_ <- function(g,mindepth=0,maxdepth=Inf){
  ret=matrix(NA,nrow=length(g$nl),ncol=2)
  nnodes=length(g$nl)
  node=g$root
  
  checkAlive <- function(i){
    C1 <- g$nl[[i]]$tmp == 1 
    if(any(is.na(g$nl[[i]]$children))) C2 <- FALSE
    else C2 <- any(sapply(g$nl[[i]]$children, function(c) g$nl[[c]]$tmp == 0))
    return(C1 && C2)}
  
  for(i in 1:length(g$nl)) g$nl[[i]]$tmp=0 # tracks completed
  mydepth=function(x){
    ifelse(is.na(x),0,min(max(mindepth,x),maxdepth))
  }
  i = 0
  fin = FALSE
  while(fin == FALSE){
    
    if(node == g$root) {t0 <- 0
    g$nl[[node]]$t = 0
    g$nl[[node]]$p = 0
    g$nl[[node]]$tmp = 1
    }else{t0 <- g$nl[[node]]$t}
    
    if(any(!is.na(g$nl[[node]]$children))){
      children <- g$nl[[node]]$children
      for(c in children){
        depth = t0 + mydepth(g$nl[[c]]$d)
        if(g$nl[[c]]$tmp == 0) {g$nl[[c]]$t = depth}
        i = i+1
        g$nl[[c]]$p = i
        g$nl[[c]]$tmp = 1
      }
    }
    alive <- (1:nnodes)[sapply(1:nnodes,checkAlive)] 
    if(length(alive) == 0) {fin = TRUE
    break}
    node = alive[1]
    next
    #for(j in alive){
    # print(paste("j:",j, "tmp:", g$nl[[j]]$tmp))
    # if(j==1 && g$nl[[1]]$tmp == 1) {fin = TRUE
    # break}
    #if(g$nl[[j]]$type != "spare" && g$nl[[j]]$tmp == 1){
    # node = j
    #break}
    #}
  }
  g = assignheights(g)
  g
}


c_get_<-function(g,n,what="d"){
  ## Get edge properties from a parent node n to either left or right child
  # if(getleft & is.na(n$cl)) return(NULL)
  #  if((!getleft) & is.na(n$cr)) return(NULL)
  

  if(any(is.na(n$children))) return(NULL)  
  if(what == "d"){
    children <- n$children
    if(n$type == "mixture"){ # Special case: mixture edges have d=0 to target node
      c1 = children[1]
      return(c(g$nl[[c1]]$d ,0))
    }
    d_vector <- sapply(children, function(c){g$nl[[c]]$d})
    return(d_vector)}
  
  if((what == "w") && (n$type != "mixture")){
    children <- n$children
    w_vector <- sapply(children, function(c){g$nl[[c]]$w[1]})
    return(w_vector)}
  
  if(what == "mix"){
    if(n$type=="mixture") {
      return(c(0,1))
    }
    else if(n$type=="split"){
      nchildren <- length(n$children)
      return(rep(0,nchildren))
    }
  }
  if((what == "w") && (n$type == "mixture")){
    m = n$mixchild
    idx = which(g$nl[[m]]$parents == n$id)
    return(c(n$w, g$nl[[m]]$w[idx]))
  }
}


edges.cg_<-function(g){
  ## Extract properties of all edges as a data frame
  edge.length=do.call("c", lapply(rev(g$nl[g$internal]),function(n){c_get_(g,n,"d")}))
  
  edge.w=do.call("c", lapply(rev(g$nl[g$internal]),function(n){c_get_(g,n,"w")}))
  
  edge=do.call("rbind",
               lapply(rev(g$nl[g$internal]),function(n){
                 tr=c()
                 if(any(!is.na(n$children))) tr= cbind(tr,n$id,n$children)
                 tr
               })
  )
  
  edge.ismix=do.call("c", lapply(rev(g$nl[g$internal]),function(n){c_get_(g,n,"mix")}))
  
  colnames(edge)=c("parent","child")
  as.data.frame(cbind(edge,
                      weight=edge.w,
                      length=edge.length,
                      #mix = rep(0,length(g$nl)-1)
                      mix=edge.ismix
  ))
}

locationasmatrix=function(g){
  ## Location of each node 
  tlayout=cbind(t=sapply(g$nl,function(x)x$t),
                p=sapply(g$nl,function(x)x$p),
                index=sapply(g$nl,function(x)x$id))
  rownames(tlayout)=tlayout[,"index"]
  tlayout
}

plot.cg_=function(g,ref=g,arrows.length=0.1,edges=NULL,
                  arrows.col=c("grey","red"),
                  text.col=c("darkgrey","darkred"),
                  digits=1,mindepth=1e-2,maxdepth=Inf,rightalign=FALSE,
                  tips=NULL,cex.edge.text=0.5,
                  label.mixture=TRUE,
                  label.nonmixture=TRUE,
                  label.internal=TRUE,
                  labels.col="black",
                  showedges=TRUE,
                  keeplocation=FALSE,
                  lwd=1,textdelta=0,
                  format=c("triangular","rectangular"),rdelta=0.1,rdelta2=0.1,vadj.edge=0.2,
                  adj.node=0,adj.edge=0,show=TRUE,showaxis=TRUE,cex.labels=1,
                  ...){
  ## Plot without requiring igraph
  if(!keeplocation) g=assignlocation_(g,mindepth,maxdepth)
  if(all(is.null(edges))){
    edges=edges.cg_(g)
    edges[,"weight"]=format(edges[,"weight"],digits=digits)
  }
  tlayout=locationasmatrix(g)
  labels=tlayout[,"index"]
  if(all(is.null(tips))) labels[1:length(g$tip.label)]=g$tip.label
  if(!all(is.null(tips))) labels[1:length(tips)]=tips
  if(!label.internal) labels[(length(g$tip.label)+1):length(labels)]=""
  if(show){
    
    plot(tlayout,type="n",xlab="",ylab="",axes=FALSE,...)
    
    if(is.na(showedges)&& any(edges[,"weight"]<0.99)) showedges=TRUE
    else if(is.na(showedges)) showedges=FALSE
    if(!label.mixture){
      edges[edges[,"mix"]==1,"weight"]=""
    }
    if(!label.nonmixture){
      edges[edges[,"mix"]==0,"weight"]=""
    }
    if(length(lwd)<dim(edges)[1])lwd=rep(lwd,dim(edges)[1])
    for(i in 1:dim(edges)[1]) {
      tt=tlayout[as.numeric(edges[i,1:2]),"t"]
      tp=tlayout[as.numeric(edges[i,1:2]),"p"]
      if(format[1]=="triangular"){
        arrows(tt[1],
               tp[1],
               tt[2],
               tp[2],
               col=arrows.col[edges[i,5]+1],lwd=lwd[i],
               length=arrows.length)
        if(showedges) text(mean(tt),mean(tp),edges[i,"weight"],adj=adj.edge,
                           cex=cex.edge.text,col=text.col[edges[i,5]+1])
      }else{
        arrows(tt[1]-rdelta2*edges[i,5],
               tp[1],
               tt[1]-rdelta2*edges[i,5],
               tp[2]+rdelta*edges[i,5],
               col=arrows.col[edges[i,5]+1],lwd=lwd[i],
               length=0)
        arrows(tt[1]-rdelta2*edges[i,5],
               tp[2]+rdelta*edges[i,5],
               tt[2],
               tp[2]+rdelta*edges[i,5],
               col=arrows.col[edges[i,5]+1],lwd=lwd[i],
               length=arrows.length)
        if(showedges) text(tt[1],tp[2]+rdelta*edges[i,5]+vadj.edge,
                           edges[i,"weight"],adj=adj.edge,
                           cex=cex.edge.text,col=text.col[edges[i,5]+1])
      }
    }
    text(tlayout[,"t"]+textdelta,tlayout[,"p"],labels=labels,col=labels.col,adj=adj.node,cex=cex.labels)
    if(showaxis) axis(1)
  }
  tlayout=cbind(as.data.frame(tlayout),label=labels)
  order=tlayout[order(tlayout[1:g$n,"p"],decreasing=F),"label"]
  invisible(list(layout=tlayout,edges=edges,order=order))
}


paths = function(edges, i, root){
  paths = list()
  path = NULL
  npaths = sum(edges$child == i)
  path_weights = edges$weight[edges$child == i]
  
  for(k in 1:npaths){
    idx = which(edges$child == i)[k]
    while (length(idx) > 0){
      E = edges[idx,]
      edge = paste0(E$child,"-", E$parent)
      r = cbind(edge, E)
      path = rbind(path,r)    
      if(E$parent == root) break
      idx = which(edges$child %in% E$parent)
    }
    path$weight = path_weights[k]
    paths[[as.character(path_weights[k])]] = path
    path = NULL
  }
  return(paths)  
}



paths2 = function(g, i){
  root = g$root
  ni = g$nl[[i]]
  paths = list()
  path = NULL
  npaths = length(ni$w) # length(n$parents)
  path_weights = ni$w
  fin = FALSE
  
  for(k in 1:npaths){
    node = ni
    weight = path_weights[k]
    
    while (fin == FALSE){
      parent = node$parents[k]
      edge =  paste0(node$id,"-",parent) 
      
      if( g$nl[[parent]]$type == "mixture") {length = 0
        }else length = node$d
      
      r = cbind(edge, as.numeric(length), as.numeric(weight))
      path = rbind(path,r)    
      
      if(parent == root){ fin = TRUE
        break}
      
      node = g$nl[[parent]]
    }
    paths[[as.character(path_weights[k])]] = path
    path = NULL
  }
  return(paths)  
}

O = function(p_i, p_j){
  overlap = p_i[p_i$edge %in% p_j$edge,]
  wi = unique(p_i$weight)
  wj = unique(p_j$weight)
  ci = overlap$length
  res = wi*wj*ci
  return(sum(res))
}


Vij = function(Paths_i, Paths_j){
  v = sum(sapply(Paths_i, function(path_i) sum(sapply(Paths_j, function(path_j) O(path_i, path_j)))))
  return(v)
}


ccov_dag_ = function(g, old=FALSE) {
  root = g$root
  ntips = g$n
  p = list(length = ntips)
  
  if(old){edges = edges.cg(g)
  }else edges = edges.cg_(g)
  
  for(i in 1:ntips){
    p[[i]] = paths(edges, i,root)
  }  
  cov = sapply(1:ntips, function(i) {
    sapply(1:ntips, function(j) {
      Vij(p[[i]], p[[j]])})
  })
  dimnames(cov) = list(g$tip.label, g$tip.label)
  return(cov)
}



get_weightmatrix = function(g) {
  edges = edges.cg_(g)
  p = list()
  for(i in g$tips) {p[[i]] = paths(edges, i, g$root)
  }
  # make dictionary mapping edges to column numbers
  dict = list()
  l = 1
  
  for(i in 1:length(g$nl)){
    node = g$nl[[i]]
    if(node$id == g$root) next
    parent = node$parents[1]
    edge =  paste0(node$id,"-",parent)
    dict[as.character(edge)] = l
    l = l + 1
  }
  W = NULL
  for(i in g$tips){
    rowvec = NULL
    for(j in 1:(length(g$nl) - 1 )) rowvec[j] = 0
    for(k in 1:length(p[[i]])){
      for (l in 1:length(p[[i]])){
        p1 = p[[i]][[k]] 
        p2 = p[[i]][[l]]
        # remove mix edges since they have drift 0
        p1 = p1[p1$mix !=1,] 
        p2 = p2[p2$mix !=1,] 
        
        overlap = p1$edge[ p1$edge %in% p2$edge] 
        index = as.numeric(dict[overlap]) 
        rowvec[index] = as.numeric(rowvec[index]) + unique(p1$weight)*unique(p2$weight)
      }
    }
    W = rbind(W,unlist(rowvec))
  }
  colnames(W) = sort(names(dict))
  return(W)
}


get_depths = function(g, V){
  #input topology g and covariance matrix and output depths
  #V = ccov_dag_(g)
  W = get_weightmatrix(g)
  sol = t(W)%*%solve(W%*%t(W))%*%diag(V)
  sol[sol < 0] = 0
  return(sol)
}


mypruneregraft_ = function(g,source,target,careful=FALSE){
  ## prunes source node branch and attaches it to the target node branch
  
  ## Run this ONLY on type=="split" nodes!
  gstart=g
  psourceroot=FALSE
  ptargetroot=FALSE
  
  ## Parent of source needs removing
  ptarget=g$nl[[target]]$parents[1] # non-admix parent of target
  psource=g$nl[[source]]$parents[1] # non-admix parent of source. This node is to be moved to the target branch
  
  ocsource = g$nl[[psource]]$children 
  ocsource = ocsource[ocsource!=source] # siblings of source node
  
  ppsource=g$nl[[psource]]$parents # grandparent of source, will become parent to other siblings of source node
  
  if(ocsource==target) {
    print(paste("source",source,"and target",target,"have same bifurcating parent",psource))
    return(g)}
    # if(length(ocsource) == 1) {print(paste("source",source,"and target",target,"have same bifurcating parent",psource))
    #   return(g)}
    # else if(!is.na(ppsource)){ # Deal with trifucating node case
    #   if(length(g$spare)>0) 
    #     print("multifurcating node to bifurcating node")
    # }    
 # }
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
    g$nl[[ocsource]]$parents[1] = ppsource
    g$nl[[ppsource]]$children[g$nl[[ppsource]]$children == psource] = ocsource
  }
  
  # ## Update the parent of the target
  if(is.na(ptarget)) print(paste("parent of target is the root ( node",ptarget),")")
  #if(ppsource != ptarget)  g$nl[[ptarget]]$children[g$nl[[ptarget]]$children == target] = psource
  g$nl[[ptarget]]$children[g$nl[[ptarget]]$children == target] = psource
  
  ## Update the split node, the old parent node
  g$nl[[psource]]$d=g$nl[[target]]$d/2
  g$nl[[psource]]$children = c(target, source)
  g$nl[[psource]]$parents[1] = ptarget
  
  ## Update the target and source
  g$nl[[target]]$d=g$nl[[target]]$d/2
  g$nl[[target]]$parents[1] = psource    
  g$nl[[source]]$parents[1] = psource # redundant?
  
  ## Done!
  if(careful) {
    if(!checkgraph(g)){
      print(paste("Moving",source,"to",target))
      print(g)
      stop("Created invalid graph!")
    }
  }
  print(g$nl[[7]]$children)
  g
}



getp_<-function(g){
  if(is(g,"cglist")) return(getp_(g[[1]]))
  ## Extracts the indices of parameters of each type
  driftp=(1:length(g$nl))
  driftp=driftp[!driftp%in%g$root]
  if(length(g$mix)>0)
    for(i in g$mix) {
      driftp=driftp[!driftp%in%g$nl[[i]]$children[1]]
    }
  mixp=g$mix
  list(drift=driftp,mix=mixp)
}


npars_<-function(g){
  ## Extract the number of each class of parameter, plus the total
  if(is(g,"cglist")){
    np=npars(g[[1]])
    return(np)
  }
  p<-getp_(g)
  ret=c(nd=length(p$drift),
        nm=length(p$mix),
        tot=0)
  ret["tot"]=ret["nd"]+ret["nm"]
  return(ret)
}


transformpars_<-function(g,pars,inv=F){
  ## Transform parameters from R to their required ranges
  if(is(g,"cglist")){
    np=npars_(g)
    res=lapply(1:length(g),function(i){
      tpars=pars[(i-1)*np["tot"]+(1:np["tot"])]
      transformpars(g[[i]],tpars,inv)
    })
    return(do.call("c",res))
  }
  np=npars_(g)
  if(np["tot"] != length(pars)) stop("Invalid parameterisation")
  for(i in 1:np["nd"]) pars[i] = driftscale(pars[i],inv=inv)
  if(np["nm"]>0) for(i in 1:np["nm"]) pars[np["nd"]+i] = mixscale(pars[np["nd"] + i],inv=inv)
  pars
}



isvalidregraftpair_<-function(g,swap){
  ## Asks: is swap[1] -> swap[2] a valid prune pair?
  done=TRUE
  ## Reject if:
  ## nodes share a parent
  if((!is.na(g$nl[[swap[1]]]$parents[1]))&&(!is.na(g$nl[[swap[2]]]$parents[1]))) if(g$nl[[swap[1]]]$parents[1]==g$nl[[swap[2]]]$parents[1]) done=FALSE
  ## The target is the parent of the source
  if((!is.na(g$nl[[swap[1]]]$parents[1])) && (g$nl[[swap[1]]]$parents[1]==swap[2])) done=FALSE
  ## The target is a descendent of the source
  if(swap[2] %in% nodesunder_(g,swap[[1]])) done=FALSE
  ## The source is a descendent of the target
  if(swap[1] %in% nodesunder_(g,swap[[2]])) done=FALSE
  ## Either is an outgroup
  if(length(g$outgroup)>0){
    tog=which(g$tip.label==g$outgroup)
    if(any(swap %in%tog)) done=FALSE
    ## The root does not have the outgroup
    if(all(g$nl[[g$root]]$children !=tog)) done=FALSE
  }
  return(done)
}


nodesunder_<-function(g,i,visited=numeric()){
  r=c()
  if(i %in% visited) return(numeric()) #stop(paste("ERROR in nodesunder: node",i,"has already been visited!"))
  if(all(!is.na(g$nl[[i]]$children))) {
    for(c in g$nl[[i]]$children){
      j=nodesunder_(g,c,visited=c(visited,i))
      r=c(r,j)
    }
  }
  r=c(r,i)
  return(unique(r))
}


gparvec_<-function(g,invtrans=FALSE){
  ## Extract the parameter vector from g
  ## Optionally transform this into R
  if(is(g,"cglist")){
    r=lapply(g,gparvec,invtrans=invtrans)
    return(do.call("c",r))
  }
  p=getp_(g)
  np=npars_(g)
  pars=numeric(np["tot"])
  for(i in 1:length(p$drift)) pars[i] = g$nl[[ p$drift[i] ]]$d 
  if(np["nm"]>0) for(i in 1:np["nm"])  pars[i+np["nd"]] = g$nl[[ g$nl[[p$mix[i]]]$cr ]]$w
  if(invtrans) pars=transformpars(g,pars,T)
  pars
}


randomregraftpair_<-function(g,maxtries=400,...){
  ## Return random nodes that are not siblings
  ## Careful not to go crazy if there are no valid options...
  done=FALSE
  ntries=0
  while(!done){
    #        print(paste("... random prune pair try",ntries)) ## DEBUG
    ret=randomnodes(g,...)
    #        print(paste(ret,collapse=",")) ## DEBUG
    for(i in 1:length(ret)) {
      while(TRUE){
        ii=ret[i]
        n=g$nl[[ii]]
        targetsparentismixture = (!is.na(n$parents[1])) && (g$nl[[n$parents[1]]]$type=="mixture")
        if(targetsparentismixture) ret[i]=n$parents[1]
        else break;
      }
    }
    done=isvalidregraftpair_(g,ret)
    ntries=ntries+1
    
    if(ntries==maxtries) stop("Error! Reached maximum number of attempts to find valid tree move!")
  }
  return(ret)
}


mypruneregraftstep_<-function(g){
  ## Do a complete prune/regraft step
  if(is(g,"cglist")){
    swap=randomregraftpair_(g[[1]])
    gtest=lapply(g,function(x)mypruneregraft_(x,swap[1],swap[2]))
    class(gtest)="cglist"
  }else{
    swap=randomregraftpair_(g)
    gtest=mypruneregraft_(g,swap[1],swap[2])
  }
  print(swap)
  list(g=gtest,swap=swap,mixswap=NA)
}


ctree_loss2_<-function(gref,dataref,center=FALSE){
  ## Evaluate the loss for a parameterised graph of class cg, NOT a cglist
  pred<-ccov_dag_(gref)
  #if(center){
  #  pred=c_Center(pred)
  #  dataref=c_Center(dataref)
  #}
  dataref=dataref[gref$tip.label,gref$tip.label]
  loss<-sum(na.omit(as.numeric((pred-dataref)^2)))
  return(loss)
}


infergraphpar_<-function(g,ctree_loss,
                         data,
                         lower=-5,
                         method="L-BFGS-B",
                         control=defaultcontrol,
                         ...){
  ## Infer a graph's best parameters
  p = gparvec_(g,invtrans=F)
  
  lower=rep(lower,length(p))
  opt=optim(p,ctree_loss2_,
            gref=g,
            method=method,
            dataref=data,lower=lower,
            control=control,...)
  g=parameterise(g,opt$par)
  list(g=g,par=opt$par,loss=opt$value)
}


myregraftmixedge_ <- function(g,i,source,target,alpha,w){
  if(is(g,"cglist")){
    ret=lapply(g,function(x){
      myregraftmixedge_(x,i,source,target,alpha,w)
    })
    class(ret)="cglist"
    return(ret)
  }
  ## Does a complete removal and replacement of a mixture edge, with specified parameters
  if(!is.na(i)) g=removemixedge_(g,i)
  g=mixedge_(g,source,target,alpha,w)
  g
}


myregraftmixedgestep_<-function(g,add=FALSE){
  if(is(g,"cglist")){
    swap=randomregraftmixture(g[[1]],add=add)
    gtest=lapply(g,function(x){
      myregraftmixedge_(x,swap$rem,swap$source,swap$target,swap$alpha,swap$w)})
    class(gtest)="cglist"
  }else{
    swap=randomregraftmixture_(g,add=add)
    gtest=myregraftmixedge_(g,swap$rem,swap$source,swap$target,swap$alpha,swap$w)
  }
  list(g=gtest,swap=NA,mixswap=swap)
}



isvalidregraftmixture_<-function(g,rem,ret){
  ## Asks if removing a node rem (can be NA for no removal) and then adding a mixture edge from ret[1] to ret[2] is valid?
  valid=TRUE
  ## Is the proposal sound? We need:
  parsource=g$nl[[ret[1]]]$parents[1]
  partarget=g$nl[[ret[2]]]$parents[1]
  ## clsource=g$nl[[ret[1]]]$cl
  ## crsource=g$nl[[ret[1]]]$cr
  ## cltarget=g$nl[[ret[2]]]$cl
  ## crtarget=g$nl[[ret[2]]]$cr
  
  add=all(is.na(rem))
  if(!add){
    if(parsource==rem) parsource=g$nl[[rem]]$parents[1]
    if(partarget==rem) partarget=g$nl[[rem]]$parents[1]
    ## Reject if (after removal of the mixture node):
    if(rem %in% ret) valid=FALSE
    g=removemixedge_(g,rem)
    rem=NA
  }
  ## nodes share a parent
  if((!is.na(parsource)) && (!is.na(partarget)) && (parsource==partarget) ) valid=FALSE
  ## The target is the parent of the source
  if((!is.na(parsource)) && (parsource==ret[2]) ) valid=FALSE
  #######  CARE HERE:
  
  ## The target is the child of the source
  if((!is.na(partarget)) && (partarget==ret[2]) ) valid=FALSE
  #######  CARE HERE:
  
  
  ## The target is the child of the source (via a mixture edge)
  if(all(!is.na(g$nl[[ret[2]]]$parents)) &&
     any(g$nl[[ret[2]]]$parents==ret[1])) valid=FALSE
  #######        
  
  
  ## The target is not already a mixture node target
#  if((!is.na(g$nl[[ret[2]]]$pr))){
#    if(is.na(rem) ||  (g$nl[[ret[2]]]$pr!=rem))
#      valid=FALSE
#  }
  ## The targets parent is not a mixture node
  ##    if((!is.na(partarget)) && (g$nl[[partarget]]$type=="mixture")) valid=FALSE
  ## The target is a descendent of the source
  ##        if(ret[2] %in% nodesunder(g,ret[[1]])) valid=FALSE 
  ## The source is a descendent of the target
  if(ret[1] %in% nodesunder_(g,ret[[2]])) valid=FALSE
  
  ## The source or the target is a forbidden outgroup
  ## Either is an outgroup
  if(length(g$outmix)>0){
    tog=which(g$tip.label==g$outmix)
    if(any(ret%in%tog)) valid=FALSE
  }    
  ## valid...
  return(valid)
}


randomregraftmixture_<-function(g,add=FALSE,maxtries=200,...){
  ntries=0
  done=FALSE
  while(!done){
    if(add){
      rem=NA
    }else{
      ## Which node to remove?
      if(length(g$mix)==1) rem=g$mix
      else rem=sample(g$mix,1)
    }
    ## Which mixture to propose?
    ret= randomnodes(g,internal = FALSE ,...)
    print(c(rem,ret))
    done=isvalidregraftmixture_(g,rem,ret)
    if(ntries==maxtries) stop("Error! Reached maximum number of attempts to find valid tree move!")
  }
  
  alpha= runif(1,0.01,0.5)
  w=runif(1,0.01,0.5)
  ret=list(rem=rem,source=ret[1],target=ret[2],alpha=alpha,w=w)
  return(ret)
}



reversemixture_<-function(g){
  ####### reverse mixture edge direction
  if(is(g,"cglist")){
    ret=lapply(g,reversemixture_)
    class(ret)="cglist"
    return(ret)
  }

  mixnodes = g$mix
  if(length(mixnodes)==0) stop("Error: g has no mix edge")
  mix = mixnodes[sample(1:length(mixnodes),1)]

  g = myregraftmixedge_(g, i= mix, source = g$nl[[mix]]$children[2], 
                        target = g$nl[[mix]]$children[1], 
                        alpha = runif(1,0.01,0.5),
                        w = g$nl[[ g$nl[[mix]]$children[2] ]]$w[2])
  list(g=g,swap=NA,mixswap=NA)
}

 

dagstep_ <- function(g,data,control=defaultcontrol,freqs=c(1/4,1/4,1/4,1/4),verbose=FALSE,...){
  ## Do one iteration of the graph
  movetype=sample(1:4,1,prob=freqs)
  if((is(g,"cglist")&&(length(g[[1]]$mix)==0)) ||
     (is(g,"cg")&&(length(g$mix)==0)))movetype=1 # No mixture edges to worry about
  if(verbose) print(paste("Proposing move of type",movetype,"..."))
  if(movetype==1){
    proposal=mypruneregraftstep_(g)
  }else if(movetype==2){        
    proposal=myregraftmixedgestep_(g)
  }else if(movetype==3){
    proposal0=mypruneregraftstep_(g)
    proposal=myregraftmixedgestep_(proposal0$g)
    proposal$swap=proposal0$swap
  }else if(movetype==4){
    proposal=reversemixture_(g)
  }else stop("Invalid move type in dagstep?!")
  if(verbose){
    print(paste("proposal: movetype",movetype,
                "swap:",paste(proposal$swap,collapse=","),
                "mixswap:",paste(proposal$mixswap,collapse=",")))
  }
  proposal$g
}


format_params <- function(g, pars){
  ## Take parameters and put them in the correct order without admix edges
  ## removes admix edge based on keys of parameter matrix
  if(length(g$mix) == 0) return(pars)
  
  keys = rownames(pars)

  nodes <- sapply(keys,substr,1,1)
  
  duplicates <- duplicated(nodes) | duplicated(nodes, fromLast = TRUE)
  
  duplicatekeys <- keys[duplicates]
  
  mixedge = duplicatekeys[substr(duplicatekeys,3,3) %in% g$mix]
  
  return(pars[rownames(pars) != mixedge , drop = FALSE])
}

parameterise_<-function(g, pars, what){
  ## Take parameters and put them in their correct place in g
  if(what == "d"){
    nodes = c(g$tips, g$internal)
    drift = nodes[nodes != g$root]
  
    #tpars = format_params(g,pars)

    for(j in 1:length(drift))  {
      g$nl[[ drift[j] ]]$d = pars[j] #tpars[j]
    }
  }
  
  if(what == "w"){
    weights = Map(c,pars,1-pars)
    mix = g$mix
    if(length(pars) != length(mix)) stop("Error: invalid weight input")
  
    if(length(mix)>0) for(j in 1:length(mix)){
      g$nl[[ g$nl[[mix[j]]]$children[2] ]]$w = weights[[j]]
    }
  }
  return(g)
}


w_loss <- function(w,g0,data){
  weights = Map(c,w,1-w)
  mix = g0$mix
  if(length(w) != length(mix)) stop("Error: invalid weight input")
  
  if(length(mix)>0) for(j in 1:length(mix)){
    g0$nl[[ g0$nl[[mix[j]]]$children[2] ]]$w = weights[[j]]
  }
  
  p = get_depths(g0, data)
  
  g = parameterise_(g0, pars=p, what ="d")
  
  loss = ctree_loss2_(g, data)
  return(loss)
}


infer_weight <- function(g, data){
  w <- runif(length(g$mix),0.05,0.95)
  opt=optimize(f = function(w) w_loss(w, g, data),
               maximum = FALSE, tol = 0.1,
               interval = c(0.05, 0.95))
  return(opt$minimum)
}


infer_graph_hc <- function(g0, data,maxiter=100,losstol=0.01, patience = 10, verbose = FALSE){
   w = infer_weight(g0, data) 
   g0 = parameterise_(g0, pars=w, what = "w")
   p = get_depths(g0, data) 
   g = parameterise_(g0, pars=p, what = "d")
 
   loss = ctree_loss2_(g,data)
   losses=rep(NA,maxiter)
   losses[1]=loss
   best_loss = loss  # Initialise best loss
   best_solution = g # Initialise best solution
   
   consecutive_loss_count = 0
   
   if(maxiter>1) for(i in 2:maxiter){
     proposal0 = dagstep_(g, data, verbose = verbose)
     w = infer_weight(proposal0, data) 
     proposal0 = parameterise_(proposal0, pars=w, what = "w")
     p = get_depths(proposal0, data) 
     proposal = parameterise_(proposal0, pars=p, what = "d")    
     
     newloss = ctree_loss2_(proposal,data)
     
     if (newloss < loss){
       print("new g")
       g = proposal
       loss = newloss
       
       # Update the best solution if a new best is found
       if (newloss < best_loss) {
         best_loss = newloss
         best_solution = proposal
       }
     }
     losses[i]=loss
     
     if (abs(losses[i] - losses[i - 1]) < losstol) {
       consecutive_loss_count = consecutive_loss_count + 1
     } else {
       consecutive_loss_count = 0
     }
     
     if (consecutive_loss_count >= patience & loss< losstol) {
       cat("Converged at iteration", i, "\n")
       break
     }
     
   }
   print(best_loss)  # Print the best loss
   return(list(best_solution, best_loss, losses))
}



# infer_graph <- function(g0, data, maxiter = 100, losstol = 0.01, patience = 5, verbose = FALSE, initial_temperature = 1, cooling_factor = 0.99) {
#   w = infer_weight(g0, data)
#   print(c(w,w,w))
#   g0 = parameterise_(g0, pars=w, what = "w")
#   p = get_depths(g0, data)
#   g = parameterise_(g0, pars=p, what = "d")
#   
#   loss = ctree_loss2_(g, data)
#   losses = rep(NA, maxiter)
#   losses[1] = loss
#   
#   consecutive_loss_count = 0
#   temperature = initial_temperature
#   
#   if (maxiter > 1) {
#     for (i in 2:maxiter) {
#       proposal0 = dagstep_(g, data, verbose = verbose)
#       w = infer_weight(proposal0, data)
#       proposal0 = parameterise_(proposal0, pars=w, what = "w")
#       p = get_depths(proposal0, data)
#       proposal = parameterise_(proposal0, pars=p, what = "d")
#       
#       newloss = ctree_loss2_(proposal, data)
#       
#       if (newloss < loss || runif(1) < exp((loss - newloss) / temperature)) {
#         print(loss)
#         g = proposal
#         loss = newloss
#       }
#       
#       losses[i] = loss
#       
#       if (abs(losses[i] - losses[i - 1]) < losstol) {
#         consecutive_loss_count = consecutive_loss_count + 1
#       } else {
#         consecutive_loss_count = 0
#       }
#       
#       if (consecutive_loss_count >= patience && loss < losstol) {
#         cat("Converged at iteration", i, "\n")
#         break
#       }
#       
#       temperature = temperature * cooling_factor
#     }
#   }
#   
#   print(loss)
#   return(list(g, loss, losses))
# }




infer_graph <- function(g0, data, maxiter = 100, losstol = 0.01, verbose = FALSE, initial_temperature = 0.5, cooling_factor = 1.5) {
  w = infer_weight(g0, data)
  g0 = parameterise_(g0, pars=w, what = "w")
  p = get_depths(g0, data)
  g = parameterise_(g0, pars=p, what = "d")
  
  loss = ctree_loss2_(g, data)
  best_loss = loss  # Initialise best loss
  best_solution = g # Initialise best solution
  losses = rep(NA, maxiter)
  losses[1] = loss
  
  temperature = initial_temperature
  
  if (maxiter > 1) {
    for (i in 2:maxiter) {
      proposal0 = dagstep_(g, data, verbose = verbose)
      w = infer_weight(proposal0, data)
      proposal0 = parameterise_(proposal0, pars=w, what = "w")
      p = get_depths(proposal0, data)
      proposal = parameterise_(proposal0, pars=p, what = "d")
      
      newloss = ctree_loss2_(proposal, data)
     # print(exp((loss - newloss) / temperature))
      if((newloss < loss ) || runif(n=1, min=0, max = 0.99) < exp((loss - newloss) / temperature)) {
        print(c("Current loss:", newloss))
        g = proposal
        loss = newloss
        
        # Update the best solution if a new best is found
        if (newloss < best_loss) {
          best_loss = newloss
          best_solution = proposal
        }
      }
      
      losses[i] = loss
      
      if (loss < losstol) {
        cat("Converged at iteration", i, "\n")
        break
      }
      
      temperature = temperature * cooling_factor
    }
  }
  
  print(best_loss)  # Print the best loss
  return(list(g = best_solution, loss = best_loss, losslist = losses))
}

