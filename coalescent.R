

setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project")
source("mixfunctions.R")
library("Clarity")


cnode_<-function(id,parents = NA,children = NA,d=NA,t=NA,p=NA,w=NA,type=NA){
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
           type=type)
    class(r)="cnode"
    return(r)
  }
}


simCoal_<-function(n,times="coal",labels=paste0("t",1:n),outgroup=numeric()){
  ## Simulate a coalescent tree in the "clarity graph" framework
  ## Returns a "clarity graph" object
  nodes = list()
  for(i in 1:n) nodes[[i]]=cnode_(i,t=0,d=0,w=1)
  alive = 1:n
  outgroupnums = alive[labels%in%outgroup]
  alive = alive[!labels%in%outgroup]
  time = 0
  fin = FALSE
  for(i in n:2){
    ii <- length(nodes) + 1    
    nChildren <- max(2,sample(1:length(alive), 1))
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


mixedge_<-function(g,source,target,alpha,w){
  ## Add a mixture edge to a clarity graph object
  
  #if(g$nl[[target]]$type == "mixture") {stop("ERROR: target cannot be a mixture node")}
  
  if(length(g$spare)>0){
    ti=tail(g$spare,1)
    g$spare=g$spare[-(g$spare%in%ti)]
  }else{
    ti=length(g$nl)+1 # new node id
  }
  ##    if(target%in%tipsunder(source)) stop("Invalid mixture edge?")
  parsource <- g$nl[[source]]$parent[1] # parent of source node
  ##    if(is.na(parsource)) stop("ERROR: tried to make the root a mixture!")
  ## Add the new node
  g$nl[[ti]] <- cnode_(ti,
                       parents = parsource,
                       children = c(source, target),
                       d=g$nl[[source]]$d*(1-alpha),
                       t=(g$nl[[parsource]]$t - g$nl[[source]]$t)*alpha +  g$nl[[source]]$t,
                       w=1,
                       type="mixture")
  ################
  ## TODO: If the target already has a right parent, we need to create a new node
  
  #if(!is.na(g$nl[[target]]$pr)){
  #  stop("Unimplemented exception: target would have three parents!")
  #}
  
  ## Update the target
  l <- length(g$nl[[target]]$parents)
  
  g$nl[[target]]=cnode(g$nl[[target]],w=g$nl[[target]]$w * w)
  
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
      #print(c)
    }
  }
  for(i in g$tips){
    c[i,i]=c[i,i]+g$nl[[i]]$d
  }
  c
}


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



tg0 <- simCoal_(4, labels=c("A","B","C","O"), times = "coal",outgroup="O")
tg0$nl[[4]]$children



tg1 = mixedge_(tg0,2,3,0.5,0.2)
tg2 = mixedge_(tg1,2,3,0.5,0.2)
tg3 <- removemixedge_(tg1,8)



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

assignlocation_ <- function(g,mindepth=0,maxdepth=Inf){
  ret=matrix(NA,nrow=length(g$nl),ncol=2)
  nnodes=length(g$nl)
  node=g$root
  for(i in 1:length(g$nl)) g$nl[[i]]$tmp=0 # tracks completed
  mydepth=function(x){
    ifelse(is.na(x),0,min(max(mindepth,x),maxdepth))
  }
  fin = FALSE
  while(fin == FALSE){
    nnodes = nnodes - 1
    if(node == g$root) {t0 <- 0
    g$nl[[node]]$t = 0
    }else{t0 <- g$nl[[node]]$t}

    
    
    if(any(!is.na(g$nl[[node]]$children))){
      children <- g$nl[[node]]$children
      for(c in children){
        depth = t0 + mydepth(g$nl[[c]]$d)
        if(g$nl[[c]]$tmp == 0) {g$nl[[c]]$t = depth}
        g$nl[[c]]$tmp = 1
      }
    }
    
    for(j in nnodes:1){
      if(j==1 && g$nl[[1]]$tmp == 1) {fin = TRUE
        break}
      if(g$nl[[j]]$type != "spare" && g$nl[[j]]$tmp == 1){
        node = j
        break}

    }
  }
  g
}


edges.cg<-function(g){
  ## Extract properties of all edges as a data frame
  edge.length=do.call("c",
                      lapply(rev(g$nl[g$internal]),function(n){
                        c(c_get(g,n,TRUE,"d"),
                          c_get(g,n,FALSE,"d"))
                      })
  )
  edge.w=do.call("c",
                 lapply(rev(g$nl[g$internal]),function(n){
                   c(c_get(g,n,TRUE,"w"),
                     c_get(g,n,FALSE,"w"))
                 })
  )
  edge=do.call("rbind",
               lapply(rev(g$nl[g$internal]),function(n){
                 tr=c()
                 if(!is.na(n$cl)) tr= rbind(tr,c(n$id,n$cl))
                 #if(!is.na(n$cr)) tr= rbind(tr,c(n$id,n$cr))
                 tr
               })
  )
  edge.ismix=do.call("c",
                     lapply(rev(g$nl[g$internal]),function(n){
                       pr=c(c_get(g,n,TRUE,"pr"),
                            c_get(g,n,FALSE,"pr"))
                       pr[is.na(pr)]=FALSE
                       as.numeric(pr==n$id)
                     })
  )
  colnames(edge)=c("parent","child")
  as.data.frame(cbind(edge,
                      weight=edge.w,
                      length=edge.length,
                      mix=edge.ismix))
}


tg0 <- simCoal_(4, labels=c("A","B","C","O"), times = "coal",outgroup="O")
x <- assignlocation_(tg0)
plot_.cg(tg0)

printTimes(x)

printFamily(x)

edges

g = simCoal(4, labels=c("A","B","C","O"), times = "coal",outgroup="O")

x <- assignlocation(g)
edges.cg(x)

g$nl[[7]]$d

printDepths(g)

printTimes(g)

printDepths(x)

printTimes(x)

plot(g)



plot_.cg=function(g,ref=g,arrows.length=0.1,edges=NULL,
                  arrows.col=c("grey","red"),
                  text.col=c("darkgrey","darkred"),
                  digits=1,mindepth=1e-2,maxdepth=Inf,rightalign=FALSE,
                  tips=NULL,cex.edge.text=0.5,
                  label.mixture=TRUE,
                  label.nonmixture=TRUE,
                  label.internal=TRUE,
                  labels.col="black",
                  showedges=NA,
                  keeplocation=FALSE,
                  lwd=1,textdelta=0,
                  format=c("triangular","rectangular"),rdelta=0.1,rdelta2=0.1,vadj.edge=0.2,
                  adj.node=0,adj.edge=0,show=TRUE,showaxis=TRUE,cex.labels=1,
                  ...){
  ## Plot without requiring igraph
  if(!keeplocation) g=assignlocation_(g,mindepth,maxdepth)
  if(all(is.null(edges))){
    edges=edges.cg(g)
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

