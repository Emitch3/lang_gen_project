
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
tinf1=infer_graph(trand1, data, maxiter = 50, movefreqs = c(1/2,1/2), verbose = T)
#tinf1=infer_graph(tinf1[[1]], data, maxiter = 100 ,verbose = T, losstol = 0.001)

plot.cg_(tg1,main="Data graph", showedges = T)
plot.cg_(tinf1$g, main="Predicted graph")

trand1 = mixedge_(tinf1$g,2,3,0.5,0.5)
plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(trand1, main="Inferred graph with admix",showedges = T )

tinf2=infer_graph(trand1, data, maxiter = 100, movefreqs = c(1/2,1/2),verbose = T)

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

## Simulate a 4 population model where there is an outgroup
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

## Perform inference
tinf1=infer_graph(trand1, data, movefreqs = c(1/2,1/2) , maxiter = 200 ,verbose =T)

plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(tg1, main="Actual graph",showedges = T )

trand1 = mixedge_(tinf1$g,2,3,0.5,0.5)
plot.cg_(tinf1$g, main = "Inferred graph",showedges = T )
plot.cg_(trand1, main="Inferred graph with admix",showedges = T )



get_depths(tinf2$g, data)
edges.cg_(tinf2$g)




tinf2=infer_graph(trand1, data, maxiter = 100, movefreqs = c(0,1), verbose = TRUE)

plot.cg_(tinf2$g, main = "Inferred graph", showedges = T)
plot.cg_(tg1, main="Actual graph",showedges = T )





data=ccov_dag_(tg1)
#data=transvector(data)*data

g = myregraftmixedge_(tinf1$g,i = NA,source = 1, target = 3,alpha = 0.1,w = 0.7)
edges.cg_(g)
edges.cg_(tinf1$g)
plot.cg_(tinf1$g)
plot.cg_(g)

w=infer_weight(g,mix = 10,data )
ng= parameterise_(g,pars = w, what = "w" )
d=get_depths(ng,V = data)
tg = parameterise_(ng,pars = d, what = "d" )
ctree_loss2_(tg,data)




data

ccov_dag_(tg)

get_weightmatrix_(g)


plot.cg_(tg)

edges.cg_(tinf2$g)

tinf2$g$nl[[18]]$d

plot.cg_(tinf2$g, main = "Inferred graph", showedges = T)
plot.cg_(tg1, main="Actual graph",showedges = T )


ng = myregraftmixedge_(trand0,i = 18, source = 1,target = 3,alpha = 0.4,w = 0.1)
  

get_depths()
infer_weight(ng, mix = ng$mix, data, 1)


ctree_loss2_(ng,dataref = data)

g = myregraftmixedge_(trand1, i = 18,source = 13, target = 12, alpha = 0.4, w = 0.2)  

g$nl[[18]]


# get_paths(edges.cg_(g),i = 1,root = g$root)


# Make a random proposal graph
trand0=simCoal(9,labels=c("A","B","C","D","E","F","G","H","O"),outgroup="O")
trand1 = mixedge(trand0,2,3,0.5,0.5)
plot.cg(trand0)
g = trand1




validPairs = 0
for (i in 1:(length(g$nl) - 1)) {
  for (j in (i + 1):length(g$nl)) {
    swap = c(i,j)
  if(isvalidregraftmixture(g,rem = NA, ret = swap)){
    print(swap)
    validPairs <- validPairs + 1
  }
  }
  print(validPairs)
}


plot(myregraftmixedge(trand1, i=18, source = 1, target = 1, alpha=0.5,w=0.2))



g=trand1
val = validMixturePairs(g)
#isvalidregraftmixture_(trand1, rem = 18, c(1,18))

for (i in 1:val$count) {
  ret = val$mixswaplist[[i]]
  g = myregraftmixedge_(trand1, i = 18,source = ret[1], target = ret[2], alpha = 0.4, w = 0.2)
  print(ret)
  #plot.cg_(g)
  
  print(ccov_dag_(g))
  #print(ret)
}


g = myregraftmixedge_(trand1, i = 18,source = 18, target = 15, alpha = 0.4, w = 0.2)  
g



plot.cg_(trand1)
plot.cg_(g)

ccov_dag_(trand1)


g = tinf1$g
plot.cg_(tg1)
residual_mixpairs(trand1, posmixturepairs(trand1, data))
g = myregraftmixedge_(tg1, i = 18, source = 15, target = 14 ,alpha = 0.4,w = 0.2) 
plot.cg_(g)
ccov_dag_(g)

#isvalidregraftmixture_(tg0, rem=NA, c(12,11))

ng = mixedge_(tg0, source =12, target = 11,alpha = 0.4, weight = 0.2)
plot.cg_(ng)

data = ccov_dag_(ng)
plot.cg_(g)

res = data - ccov_dag_(g)
plotres(g,data)
corr = cor(data) - cor(ccov_dag_(g))

corrplot::corrplot(corr)
corrplot::corrplot(res)

idx = sample(seq_len(36), size = 1, replace = TRUE, prob = correlations)

nodesabove(g, 1)
nodesabove(g,2)
nodesabove(g,7)

tipsunder_(g,16)

nodesunder_(g,16)

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

nearestTip(g,16)

pairs[,idx]

correlations

# Number of elements in your vector
n <- length(correlations)

# Apply Laplace smoothing
smooth_likelihoods <- (correlations + 0.01) / (sum(correlations) + n*0.01)

# Print the new probability distribution
print(smooth_likelihoods)

idx = sample(seq_len(36), size = 1, replace = TRUE, prob = smooth_likelihoods)
pairs[,idx]

sample(df, size = 1, replace = TRUE, prob = df$correlations)

#tinf1=infer_graph(tinf1$g, data, maxiter = 100 ,verbose = T, losstol = 0.01)

plot.cg_(tinf2$g, main = "Inferred graph", showedges = T)
plot.cg_(tg1, main="Actual graph",showedges = T )

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


plot.cg_(tg1)
plot.cg_(trand1)
corrplot::corrplot(res)
plotres(trand1, data)


flattened_matrix <- res[lower.tri(res)]

top_n_indices <- tail(order(flattened_matrix),1)

row_col_indices <- which(lower.tri(res), arr.ind = TRUE)[top_n_indices, ]
row_names <- rownames(res)[row_col_indices[, 1]]
col_names <- colnames(res)[row_col_indices[, 2]]

result <- data.frame(
  row = row_names,
  column = col_names,
  value = flattened_matrix[top_n_indices]
)

result

res = data - ccov_dag_(tinf2$g)

plotres(tinf2$g, data)
plot.cg_(tinf2$g)
flat= res[lower.tri(res)]

low <- median(flat) - 3 * mad(flat, constant = 1)
high <- median(flat) + 3 * mad(flat, constant = 1)

outliers <- subset(flat, flat < low | flat > high)

outliers

sort(outliers,decreasing = T)
max(outliers)


which.max(res)

max_index <- which(res == max(res), arr.ind = TRUE)
rownames(max_index)


row_name <- rownames(res)[max_index[1]]
col_name <- colnames(res)[max_index[2]]


max_index <- which(res == max(res), arr.ind = TRUE)




print(paste("Row Name: ", row_name))
print(paste("Column Name: ", col_name))



corrplot::corrplot(cor(data) -cor(ccov_dag_(trand1)) )


plotcov(tg1)
plotcov(trand1)

Clarity_Chart(ordermatrix(tg1,pt$order),scalefun=myscalefun2,text=T )

g = trand1
pred <- ccov_dag_(g)
res <- dataref - pred

res
qr.default(res)
data

cor

# Positive residuals indicate pairs of populations where the model underestimates 
# the observed covariance, and thus populations where the fit might be improved by 
# adding additional edges.

# Negative residuals indicate pairs of populations where the model overestimates 
# the observed covariance; these are a necessary outcome of having positive residuals,
# but can also sometimes be interpreted as populations that are forced too close together
# due to unmodeled migration elsewhere in the graph.





tr_mean <- function(mat){
  m = nrow(mat)
  sum = 0
  for(i in 1:m){
    for (j in i:m) {
      sum = sum + mat[i,j]
    }
  }
  return(sum/(m*(m-1)/2))
}

# f <- function(g, data){
#   W = ccov_dag_(g)
#   W_hat = tr_mean(W)
#   R = W - data
#   R_hat = tr_mean(R)
#   
#   return(1 - sum((tr_flatten(W) - W_hat)^2)/sum((tr_flatten(R) - R_hat)^2) )
#   
# }

f(trand1, data)
ctree_loss2_(trand1, data)

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

