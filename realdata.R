library("ape")
library("Clarity")
source("mixfunctions.R")

importseparateres=function(x,n,end=".RData"){
    lapply(n,function(i){
        tn=paste0(x,i,end)
        ret=readRDS(tn)
        if(!is(ret,"cglist")){
            if(is(ret[[1]],"cglist")) return(ret[[1]])
            else stop("Invalid data?")
        }else{
            return(ret)
        }
    })
}

ordermatrix<-function(x,order){
    x=x[rownames(x)%in%order,colnames(x)%in%order]
    mat=matrix(NA,nrow=length(order),ncol=length(order))
    rownames(mat)=colnames(mat)=order
    mat[rownames(x),colnames(x)]=x
    mat
}

genemat=as.matrix(read.csv("easttimor/covgene.csv",row.names=1))
lingmat=as.matrix(read.csv("easttimor/ling.csv",row.names=1))
phonmat=as.matrix(read.csv("easttimor/phon.csv",row.names=1))
paintmat=as.matrix(read.csv("easttimor/paintinglab.csv",row.names=1))
genematp=genemat-min(genemat,na.rm=TRUE)
paintmatp=paintmat-min(paintmat,na.rm=TRUE)
lingmatp=lingmat-min(lingmat,na.rm=TRUE)
phonmatp=phonmat-min(phonmat,na.rm=TRUE)
mydatap=list(gene=genematp,paint=paintmatp,ling=lingmatp,phon=phonmatp)
n=dim(genematp)[1]

names=read.table("names.txt")
nm=names[,2]
lg=names[,3]
names(nm)=names(lg)=names[,1]
lg=as.factor(lg)

## genematp=genemat-min(genemat,na.rm=TRUE)
## lingmatp=lingmat-min(lingmat,na.rm=TRUE)
## mydatap=list(gene=genematp,ling=lingmatp)
origorder=c("EAS_CHB","H","G","A","D","C","E","B","F")
mydatap2=lapply(mydatap,function(x)x[origorder,origorder])

## mydatap2=list(gene=mydatap[[1]][origorder,origorder],
##               paint=mydatap[[2]][origorder,origorder],
##               ling=mydatap[[3]][origorder,origorder],
##               phon=mydatap[[4]][origorder,origorder]
##               )

###################
###################
###################
###################
#### Inference with 9 tips and one dataset
runsim=FALSE
if(runsim){
## pruning and regrafting with admixture
set.seed(1)
n=9
g0=simCoal(n,times="test")
truth0=ccov_dag(g0)
noise0=matrix(rnorm(g0$n*g0$n,mean=0,sd=0.01),nrow=g0$n,ncol=g0$n)
data0=truth0 + noise0
truepars0=gparvec(g0)
np0=length(truepars0)

## sapply(g$nl,function(x)x$d)
gm=mixedge(g0,source=1,target=6,1/4,0.4)
allw=sapply(gm$mix,cweight,g=gm) # weights
gmtrees=clistoftrees(gm)

##### Set up the data
g=gm
truth=ccov_dag(g)
set.seed(1)
noise=matrix(rnorm(g$n*g$n,mean=0,sd=0.01),nrow=g$n,ncol=g$n)
for(i in 1:(g$n-1))for(j in (i+1):g$n) noise[i,j]=noise[j,i]
data=truth + noise
truepars=gparvec(g)
np=length(truepars)

#########
set.seed(3)
ginfr0=randomgraphlist(n,0,1)
ginfr1=myregraftmixedgestep(ginfr0,TRUE)
ginfr2=myregraftmixedgestep(ginfr1$g,TRUE)

ginfresr0=infer_dag(ginfr0,data,init="random",maxiter=50,losstol=0.01)
ginfresr1=infer_dag(ginfr1$g,data,init="random",maxiter=50,losstol=0.01)
set.seed(3)
ginfresr2=infer_dag(ginfr2$g,data,init="random",maxiter=50,losstol=0.01)

###################
###################
###################
###################
##Extension to dual datasets


## Truth has zero mixture edges
plist0=c(runif(length(gparvec(g0))),
     runif(length(gparvec(g0))))
glist0=list(g0,g0)
class(glist0)="cglist"
glist0=parameterise(glist0,plist0)

truthlist0=ccov_dag(glist0)
datalist0=lapply(1:length(truthlist0),function(x){
    noise0=matrix(rnorm(n*n,mean=0,sd=0.01),nrow=n,ncol=n)
    for(i in 1:(n-1))for(j in (i+1):n) noise0[i,j]=noise0[j,i]
    truthlist0[[x]] + noise0
})
trueparlist0=gparvec(glist0)

## Truth has one mixture edge
g1=gm
plist1=c(runif(length(gparvec(g1))),
     runif(length(gparvec(g1))))
glist1=list(g1,g1)
class(glist1)="cglist"
glist1=parameterise(glist1,plist1)

truthlist1=ccov_dag(glist1)
datalist1=lapply(1:length(truthlist1),function(x){
    noise1=matrix(rnorm(n*n,mean=0,sd=0.01),nrow=n,ncol=n)
    for(i in 1:(n-1))for(j in (i+1):n) noise1[i,j]=noise1[j,i]
    truthlist1[[x]] + noise1
})
trueparlist1=gparvec(glist1)

## Inference
set.seed(3)
ginfrlist0=randomgraphlist(n,0,2)
class(ginfrlist0)="cglist"
ginfrlist1=list(ginfr1,ginfr1)
class(ginfrlist1)="cglist"
ginfrlist2=list(ginfr2,ginfr2)
class(ginfrlist2)="cglist"

ginfresr0=infer_dag(ginfrlist0,datalist0,init="random",maxiter=20,losstol=0.01)
ginfresr1=infer_dag(ginfr1$g,data,init="random",maxiter=50,losstol=0.01)
set.seed(3)
ginfresr2=infer_dag(ginfr2$g,data,init="random",maxiter=50,losstol=0.01)

par(mfrow=c(2,2))
plot(glist0)
plot(ginfresr0$g)

#############################################
#############################################
#############################################
#############################################
#############################################
#############################################
#############################################
## REAL DATA!

testinit=randomgraphlist(n,0,4)
testinf=infer_dag(testinit,mydatap,init="random",maxiter=2,losstol=0.01)

gmydata0=infer_dag(ginfrlist0,mydatap,init="random",maxiter=20,losstol=0.01)
gmydata1=infer_dag(ginfrlist1,mydatap,init="random",maxiter=20,losstol=0.01)

gmypred0<-ccov_dag(gmydata0$g)

tips=rownames(mydatap[[1]])

library("Clarity")
par(mfrow=c(3,2))
Clarity_Chart(mydatap[[1]],scalefun=I,main="Genetic data")
Clarity_Chart(mydatap[[2]],scalefun=I,main="Linguistic data")
plot(gmydata0$g,main="Inferred DAG (0 mixtures)")
Clarity_Chart(gmypred0[[1]],scalefun=I,main="Genetic predicted")
Clarity_Chart(gmypred0[[2]],scalefun=I,main="Linguistic predicted")

mydatap_m0_N2000_r8=readRDS("old_k8/mydatap_m0_N2000_r8.RData")
mydatap_m0_N100_r8=readRDS("old_k8/mydatap_m0_N100_r8.RData")
mydatap_m1_N1000_r8=readRDS("old_k8/mydatap_m1_N1000_r8.RData")
mydatap_m3_N100_r8=readRDS("old_k8/mydatap_m3_N100_r8.RData")

mydatap_m1_N1000_r8=importseparateres("old_k8/tmp_mydatap_m1_N1000_r8_s",1:8,".Rdata")
mydatap_m2_N1000_r8=importseparateres("old_k8/tmp_mydatap_m2_N1000_r8_s",1:8,".Rdata")
mydatap_m3_N1000_r8=importseparateres("old_k8/tmp_mydatap_m3_N1000_r8_s",1:8,".Rdata")
mydatap_m4_N1000_r8=importseparateres("old_k8/tmp_mydatap_m4_N1000_r8_s",1:8,".Rdata")

g=mydatap_m1_N1000_r8[[1]]
ccov_dag(g)
gg=infer_dag(g,mydatap2,maxiter=1)
par(mfrow=c(2,2))
plot(gg$g)
ctree_loss(gparvec(gg$g,invtrans=FALSE),gg$g,mydatap2,transform=FALSE)



mydatap_m1_N1000_r8=mydatap_m1_N1000_r8[
    order(sapply(mydatap_m1_N1000_r8,function(x)tail(x$loss,1))) ]
mydatap_m2_N1000_r8=mydatap_m2_N1000_r8[
    order(sapply(mydatap_m2_N1000_r8,function(x)tail(x$loss,1))) ]
mydatap_m3_N1000_r8=mydatap_m3_N1000_r8[
    order(sapply(mydatap_m3_N1000_r8,function(x)tail(x$loss,1))) ]
mydatap_m4_N1000_r8=mydatap_m4_N1000_r8[
    order(sapply(mydatap_m4_N1000_r8,function(x)tail(x$loss,1))) ]

mydatap_m3_N100_r8=mydatap_m3_N100_r8[
    order(sapply(mydatap_m3_N100_r8,function(x)tail(x$losses,1))) ]
mydatap_m1_N1000_r8=mydatap_m1_N1000_r8[
    order(sapply(mydatap_m1_N1000_r8,function(x)tail(x$losses,1))) ]
mydatap_m0_N100_r8=mydatap_m0_N100_r8[
    order(sapply(mydatap_m0_N100_r8,function(x)tail(x$losses,1))) ]
mydatap_m0_N2000_r8=mydatap_m0_N2000_r8[
    order(sapply(mydatap_m0_N2000_r8,function(x)tail(x$losses,1))) ]
sapply(mydatap_m0_N100_r8,function(x)tail(x$losses,1))
sapply(mydatap_m0_N2000_r8,function(x)tail(x$losses,1))
sapply(mydatap_m1_N1000_r8,function(x)tail(x$losses,1))
sapply(mydatap_m3_N100_r8,function(x)tail(x$losses,1))

mydatap_m1_N1000_r8glist=lapply(mydatap_m1_N1000_r8,function(x)x)
mydatap_m2_N1000_r8glist=lapply(mydatap_m2_N1000_r8,function(x)x)
mydatap_m3_N1000_r8glist=lapply(mydatap_m3_N1000_r8,function(x)x)
mydatap_m4_N1000_r8glist=lapply(mydatap_m4_N1000_r8,function(x)x)
mydatap_m1_N1000_r8cov=lapply(mydatap_m1_N1000_r8,function(x)ccov_dag(x))
mydatap_m2_N1000_r8cov=lapply(mydatap_m2_N1000_r8,function(x)ccov_dag(x))
mydatap_m3_N1000_r8cov=lapply(mydatap_m3_N1000_r8,function(x)ccov_dag(x))
mydatap_m4_N1000_r8cov=lapply(mydatap_m4_N1000_r8,function(x)ccov_dag(x))


mydatap_m0_N500_r8glist=lapply(mydatap_m0_N2000_r8,function(x)x$g)
mydatap_m1_N1000_r8glist=lapply(mydatap_m1_N1000_r8,function(x)x$g)
mydatap_m3_N100_r8glist=lapply(mydatap_m3_N100_r8,function(x)x$g)
mydatap_m0_N2000_r8cov=lapply(mydatap_m0_N2000_r8,function(x)ccov_dag(x$g))
mydatap_m1_N1000_r8cov=lapply(mydatap_m1_N1000_r8,function(x)ccov_dag(x$g))
mydatap_m3_N100_r8cov=lapply(mydatap_m3_N100_r8,function(x)ccov_dag(x$g))

## Plot the graphs
par(mfrow=c(3,4))
for(i in 1:6){
    plot(mydatap_m0_N500_r8glist[[i]],tips=tips)
}

par(mfrow=c(3,4))
for(i in 1:6){
    plot(mydatap_m1_N1000_r8glist[[i]],tips=tips)
}

par(mfrow=c(3,4))
for(i in 1:6){
    plot(mydatap_m3_N100_r8glist[[i]],tips=tips)
}

## Plot the data and predictions
par(mfrow=c(3,4))
Clarity_Chart(mydatap[[1]],scalefun=I,text=T,main="Genetic data")
Clarity_Chart(mydatap[[2]],scalefun=I,text=T,main="Linguistic data")
for(i in 1:2){
    Clarity_Chart(mydatap_m1_N100_r8cov[[i]][[1]],text=T,scale=I,main=paste("genetics prediction",i))
    Clarity_Chart(mydatap_m1_N100_r8cov[[i]][[2]],text=T,scale=I,main=paste("Ling prediction",i))
    plot(mydatap_m3_N100_r8glist[[i]],tips=tips)
}

## Plot the data and predictions
par(mfrow=c(3,4))
Clarity_Chart(mydatap[[1]],scalefun=I,text=T,main="Genetic data")
Clarity_Chart(mydatap[[2]],scalefun=I,text=T,main="Linguistic data")
for(i in 1:5){
    Clarity_Chart(mydatap_m3_N100_r8cov[[i]][[1]],text=T,scale=I,main=paste("genetics prediction",i))
    Clarity_Chart(mydatap_m3_N100_r8cov[[i]][[2]],text=T,scale=I,main=paste("Ling prediction",i))
}


Clarity_Chart(mydatap_m0_N500_r8cov[[i]][[1]],text=T,scale=I,main=paste("genetics prediction",i))
Clarity_Chart(mydatap_m0_N500_r8cov[[i]][[2]],text=T,scale=I,main=paste("Ling prediction",i))


##############
## Write the files for m=1
numstarts=8
nummixtures=1
numiters=1000
glist=get(paste0("mydatap_m",nummixtures,"_N",numiters,"_r",numstarts,"glist"))
toutfiles=lapply(1:numstarts,function(i){
    paste0("tmp_mydatap_m",nummixtures,"_N",numiters,"_r",numstarts,"_s",i,".RData")})
for(i in 1:length(glist)){
    saveRDS(glist[[i]],file=toutfiles[[i]])
}

#############
## Write the files for m=2
## We're going to take the best 4 from m=1 and add an edge, then the best 4 from m=3 and remove one
numstarts=8
nummixtures=1
numiters=1000
glist1=get(paste0("mydatap_m",nummixtures,"_N",numiters,"_r",numstarts,"glist"))
glist3=get(paste0("mydatap_m",3,"_N",100,"_r",numstarts,"glist"))
glist2a=lapply(glist3,function(x){ ## Remove the same edge from both graphs
    for(i in 1:length(x)) x[[i]]=removemixedge(x[[i]],18)
    x
})
glist2b=lapply(glist3,function(x){ ## Remove the same edge from both graphs
    myregraftmixedgestep(x,TRUE)
})
glist=glist2a
glist[5:8]=glist2b[1:4]
nummixtures=2
toutfiles=lapply(1:numstarts,function(i){
    paste0("tmp_mydatap_m",nummixtures,"_N",numiters,"_r",numstarts,"_s",i,".RData")})
for(i in 1:length(glist)){
    saveRDS(glist[[i]],file=toutfiles[[i]])
}

#############
## Write the files for m=3
numstarts=8
nummixtures=3
numiters=1000
glist=get(paste0("mydatap_m",nummixtures,"_N",100,"_r",numstarts,"glist"))
toutfiles=lapply(1:numstarts,function(i){
    paste0("tmp_mydatap_m",nummixtures,"_N",numiters,"_r",numstarts,"_s",i,".RData")})
for(i in 1:length(glist)){
    saveRDS(glist[[i]],file=toutfiles[[i]])
}

#############
## Write the files for m=4
numstarts=8
nummixtures=3
numiters=1000
glist=get(paste0("mydatap_m",nummixtures,"_N",100,"_r",numstarts,"glist"))
nummixtures=4
glist=lapply(glist,function(x){ ## Remove the same edge from both graphs
    myregraftmixedgestep(x,TRUE)
})
toutfiles=lapply(1:numstarts,function(i){
    paste0("tmp_mydatap_m",nummixtures,"_N",numiters,"_r",numstarts,"_s",i,".RData")})
for(i in 1:length(glist)){
    saveRDS(glist[[i]],file=toutfiles[[i]])
}

    
} ## END ifrunsim

#############
#############
#############
#############
#############
#############
#############
## Starting from treemix
library("Clarity")
source("mixfunctions.R")
td="~/OD/Projects/clarity/genelanguage_2021/EastTimorComplete/EastTimorOutgroup/"
tm2=readtm(paste0(td,"tim_fs_og_pruned.stem.m5"))
g2=tm2cg(tm2)#,outgroup="EAS_CHB")
checkgraph(g2)

g2$outmix=which(g2$tip.label=="EAS_CHB")
g2$outgroup=which(g2$tip.label=="EAS_CHB")
mixtargets=cbind(sapply(g2$mix,function(x)g2$nl[[x]]$cl),
                 sapply(g2$mix,function(x)g2$nl[[x]]$cr))
mixtargets[is.na(mixtargets)]=0
remedges=g2$mix[which(apply(mixtargets==g2$outmix,1,any))]

while(any(!is.na(remedges))){
    remedges=sort(remedges,T)
    for(e in remedges){
        g2=removemixedge(g2,e)
    }
    mixtargets=cbind(sapply(g2$mix,function(x)g2$nl[[x]]$cl),
                     sapply(g2$mix,function(x)g2$nl[[x]]$cr))
    mixtargets[is.na(mixtargets)]=0
    remedges=g2$mix[which(apply(mixtargets==g2$outmix,1,any))]
}
checkgraph(g2)

test2=infer_dag(g2,mydatap$gene,init="graph",maxiter=1,control=list(maxit=200,trace=1,factr=1e12,pgtol=0))

pg2=plot(test2$g,mindepth=0.0002)

g2list=list(g2,g2,g2,g2)
checkgraph(g2)
class(g2list)="cglist"
mydatatm=list(snp=ordermatrix(mydatap$gene,g2$tip.label),
              paint=ordermatrix(mydatap$paint,g2$tip.label),
              lex=ordermatrix(mydatap$ling,g2$tip.label),
              phon=ordermatrix(mydatap$phon,g2$tip.label)
              )
gtm2a=infer_dag(g2list,mydatatm,init="graph",maxiter=1,control=list(maxit=1000,trace=1,factr=1e12,pgtol=0))
gtm2b=infer_dag(gtm2a$g,mydatatm,init="graph",maxiter=1,control=list(maxit=1000,trace=1,factr=1e10,pgtol=0))
##gtm2c=infer_dag(gtm2b$g,mydatatm,init="graph",maxiter=1,control=list(maxit=1000,trace=1,factr=1e8,pgtol=0))
gtm2=gtm2b
## gparvec(g2list[[2]],invtrans=FALSE)
## gparvec(gtm2$g[[2]],invtrans=FALSE)

## g0=gtm2$g
## p=gparvec(g0[[1]],invtrans=TRUE)
## gref<-parameterise(g0,p,transform=FALSE)

## last=infergraphpar(g0,ctree_loss,mydatatm,control=defaultcontrol)


gpred=ccov_dag(gtm2$g)
rownames(gpred[[1]])=colnames(gpred[[1]])=rownames(gpred[[2]])=colnames(gpred[[2]])=gtm2$g[[1]]$tip.label

myscalefun=function(x)sign(x)*sqrt(abs(x))

pdf("ET_GeneLanguageExample.pdf",height=20,width=15)
layout(matrix(1:12,ncol=3,nrow=4))
plot(gtm2$g,
     main=c("Genes (SNPs) Graph","Genes (Painting) Graph","Language (Lex) Graph", "Language (Phon) Graph"),
     mindepth=1,maxdepth=1,digits=2)
Clarity_Chart(ordermatrix(gpred[[1]],pg2$order),
              text=T,scalefun=myscalefun,
              main="Genes (SNPs) Predicted",digits=4)
Clarity_Chart(ordermatrix(gpred[[2]],pg2$order),
              text=T,scalefun=myscalefun,
              main="Genes (Painting) Predicted",digits=4)
Clarity_Chart(ordermatrix(gpred[[3]],pg2$order),
              text=T,scalefun=myscalefun,
              main="Language (Lex) Predicted",digits=3)
Clarity_Chart(ordermatrix(gpred[[4]],pg2$order),
              text=T,scalefun=myscalefun,
              main="Language (Phonetic) Predicted",digits=3)
Clarity_Chart(ordermatrix(genematp,pg2$order),
              text=T,scalefun=myscalefun,
              main="Genes (SNPs)",digits=4)
Clarity_Chart(ordermatrix(paintmatp,pg2$order),
              text=T,scalefun=myscalefun,
              main="Genes (Painting)",digits=4)
Clarity_Chart(ordermatrix(lingmatp,pg2$order),
              text=T,scalefun=myscalefun,
              main="Language (Lex)",digits=3)
Clarity_Chart(ordermatrix(phonmatp,pg2$order),
              text=T,scalefun=myscalefun,
              main="Language (Phon)",digits=3)
dev.off()


par(mfrow=c(1,3))
plot(g2,mindepth=0.0002)
Clarity_Chart(ordermatrix(mydatap$gene,pg2$order),
              text=T,scalefun=myscalefun,
              main="Genes",digits=4)
Clarity_Chart(ordermatrix(mydatap$ling,pg2$order),
              text=T,scalefun=myscalefun,
              main="Language",digits=3)

system("mkdir -p tmpgraphs")
saveRDS(gtm2$g,"tmpgraphs/gtm3.RData")
