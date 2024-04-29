
setwd("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project")
#source("mixfunctions.R")
install.packages("remotes")
remotes::install_github("danjlawson/CLARITY/Clarity")
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


create_chart <- function(mydata, datapred, scalefun, text, main) {
  # Calculate the difference
  d = mydata - datapred #ifelse(mydata != 0, (mydata - datapred)/mydata, NA)
  
  d[1,] = NA
  d[,1] = NA
  # Calculate the percentage and replace NA values with 0
  percentage <- 100 * round(d,3) #replace(d, is.na(d), 0)
  
  # Call the Clarity_Chart function
  Clarity_Chart(
    percentage,
    scalefun = scalefun,
    text = text,
    main = main
  )
}




################################################################################

### Real world data ###

# Load data
#load("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Realdata inference.RData")

genemat=as.matrix(read.csv("easttimor/covgene.csv",row.names=1))
lingmat=as.matrix(read.csv("easttimor/ling.csv",row.names=1))
phonmat=as.matrix(read.csv("easttimor/phon.csv",row.names=1))
paintmat=as.matrix(read.csv("easttimor/paintinglab.csv",row.names=1))

genematp=genemat-min(genemat,na.rm=TRUE)
paintmatp=paintmat-min(paintmat,na.rm=TRUE)
lingmatp=lingmat-min(lingmat,na.rm=TRUE)
phonmatp=phonmat-min(phonmat,na.rm=TRUE)


lingmatp[1,] = 0
lingmatp[,1] = 0
lingmatp[1,1] = lingmatp[2,2]
phonmatp[1,] = 0
phonmatp[,1] = 0
phonmatp[1,1] = phonmatp[2,2]



mydatap=list(gene=genematp,paint=paintmatp,ling=lingmatp,phon=phonmatp)
n=dim(genematp)[1]

names=read.table("names.txt")
nm=names[,2]
lg=names[,3]
names(nm)=names(lg)=names[,1]
lg=as.factor(lg)


origorder=c("EAS_CHB","Mam","Tok","Bun","Kem","Wai","Mak","Fat","Tet")

for(i in 1:4){
  rownames(mydatap[[i]]) = origorder
  colnames(mydatap[[i]]) = origorder
}
mydatap2=lapply(mydatap,function(x) x[origorder,origorder])


# A     B     E     C     D     F     G     H 
# "Bun" "Fat" "Mak" "Wai" "Kem" "Tet" "Tok" "Mam" 




#### Inference with no admixture ###############################################

## Simulate a 9 population model where there is an outgroup - EAS_CHB
tg0=simCoal_(9 ,labels=c("EAS_CHB","Mam","Tok","Bun","Kem","Wai","Mak","Fat","Tet"), outgroup="EAS_CHB")
plot.cg_(tg0)

ginfrlist0 = list(tg0, tg0, tg0, tg0)
class(ginfrlist0) = "cglist"

ginf0 = infer_graph(g = ginfrlist0, data = mydatap2, maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)

par(mfrow=c(2,2))
lapply(ginf0$g, plot.cg_)
predcov0 = ccov_dag_(ginf0$g)

par(mfrow=c(1,2))



plot.cg_(ginf0$g[[1]], main = "Genes (SNPs) Graph")
create_chart(mydatap2[[1]], predcov0[[1]],I, TRUE, "(100 x) Residual covariance matrix (SNPs)")

plot.cg_(ginf0$g[[2]], main = "Genes (Painting) Graph")
create_chart(mydatap2[[2]], predcov0[[2]],I, TRUE, "(100 x) Residual covariance matrix (Painting)")

plot.cg_(ginf0$g[[3]], main = "Language (Lex) Graph")
create_chart(mydatap2[[3]], predcov0[[3]],I, TRUE, "(100 x) Residual covariance matrix (Lex)")

plot.cg_(ginf0$g[[4]], main = "Language (Phonetics) Graph")
create_chart(mydatap2[[4]], predcov0[[4]],I, TRUE, "(100 x) Residual covariance matrix (Phonetics)")


########

# Plot predicted covariance matrices

Clarity_Chart(predcov0[[1]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (SNPs)")

Clarity_Chart(predcov0[[2]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Painting)")

Clarity_Chart(predcov0[[3]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Lexical)")

Clarity_Chart(predcov0[[4]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Phonetics)")


#### Inference with 1 admixture edges ##########################################

# n_mix = 1


tg1 = add_mixedges(g = ginf0$g, n_edges = 1)
plot.cg_(tg1[[1]])

ginf1 = infer_graph(g = ginf1_$g, data = mydatap2, maxiter = 500, losstol = 1e-30, MHtol = 0, verbose = T,plotprogress = T)

plot.cg_(ginf1$g[[1]], main = "Genes (SNPs) Graph")
plot.cg_(ginf1$g[[2]], main = "Genes (Painting) Graph")
plot.cg_(ginf1$g[[3]], main = "Language (Lex) Graph")
plot.cg_(ginf1$g[[4]], main = "Language (Phonetics) Graph")

par(mfrow=c(1,2))
plot.cg_(ginf1$g[[1]], main = "Genes (SNPs) Graph")
create_chart(mydatap2[[1]], predcov1[[1]],I, TRUE, "(100 x) Residual covariance matrix (SNPs)")

plot.cg_(ginf1$g[[2]], main = "Genes (Painting) Graph")
create_chart(mydatap2[[2]], predcov1[[2]],I, TRUE, "(100 x) Residual covariance matrix (Painting)")

plot.cg_(ginf1$g[[3]], main = "Language (Lex) Graph")
create_chart(mydatap2[[3]], predcov1[[3]],I, TRUE, "(100 x) Residual covariance matrix (Lex)")

plot.cg_(ginf1$g[[4]], main = "Language (Phonetics) Graph")
create_chart(mydatap2[[4]], predcov1[[4]],I, TRUE, "(100 x) Residual covariance matrix (Phonetics)")

# Plot predicted covariance matrices

Clarity_Chart(predcov1[[1]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (SNPs)")

Clarity_Chart(predcov1[[2]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Painting)")

Clarity_Chart(predcov1[[3]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Lexical)")

Clarity_Chart(predcov1[[4]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Phonetics)")


#### Inference with 2 admixture edges #####################################

# n_mix = 2

tg2 = add_mixedges(g = ginfrlist0, n_edges = 2)

# tg2 = add_mixedges(g = ginf1$g, n_edges = 1)
plot.cg_(tg2[[1]])

ginf2_2 = infer_graph(g = tg2, data = mydatap2, maxiter = 600, losstol = 1e-30,verbose = T,plotprogress = T)

predcov2 = ccov_dag_(ginf2$g)

# save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Realdata inference.RData")

plot.cg_(ginf2$g[[1]], main = "Genes (SNPs) Graph")
plot.cg_(ginf2$g[[2]], main = "Genes (Painting) Graph")
plot.cg_(ginf2$g[[3]], main = "Language (Lex) Graph")
plot.cg_(ginf2$g[[4]], main = "Language (Phonetics) Graph")

####
par(mfrow=c(1,2))

plot.cg_(ginf2$g[[1]], main = "Genes (SNPs) Graph")
create_chart(mydatap2[[1]], predcov2[[1]],I, TRUE, "(100 x) Residual covariance matrix (SNPs)")

plot.cg_(ginf2$g[[2]], main = "Genes (Painting) Graph")
create_chart(mydatap2[[2]], predcov2[[2]],I, TRUE, "(100 x) Residual covariance matrix (Painting)")

plot.cg_(ginf2$g[[3]], main = "Language (Lex) Graph")
create_chart(mydatap2[[3]], predcov2[[3]],I, TRUE, "(100 x) Residual covariance matrix (Lex)")

plot.cg_(ginf2$g[[4]], main = "Language (Phonetics) Graph")
create_chart(mydatap2[[4]], predcov2[[4]],I, TRUE, "(100 x) Residual covariance matrix (Phonetics)")



#######

# Plot predicted covariance matrices

Clarity_Chart(predcov2[[1]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (SNPs)")

Clarity_Chart(predcov2[[2]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Painting)")

Clarity_Chart(predcov2[[3]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Lexical)")

Clarity_Chart(predcov2[[4]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Phonetics)")



#### Inference with 3 admixture edges ######################################

tg3 = add_mixedges(g = ginfrlist0, n_edges = 3)
# n_mix = 3

tg3 = add_mixedges(g = ginf2$g, n_edges = 1)
plot.cg_(tg3[[1]])



ginf3 = infer_graph(g = tg3, data = mydatap2, maxiter = 600, losstol = 1e-30,verbose = T,plotprogress = T)

predcov3 = ccov_dag_(ginf3$g)

save.image("C:/Users/USER/OneDrive/Bristol Year 4/Project/lang_gen_project/Realdata inference.RData")


plot.cg_(ginf3$g[[1]], main = "Genes (SNPs) Graph")
plot.cg_(ginf3$g[[2]], main = "Genes (Painting) Graph")
plot.cg_(ginf3$g[[3]], main = "Language (Lex) Graph")
plot.cg_(ginf3$g[[4]], main = "Language (Phonetics) Graph")

par(mfrow=c(1,2))

plot.cg_(ginf3$g[[1]], main = "Genes (SNPs) Graph")
create_chart(mydatap2[[1]], predcov3[[1]],I, TRUE, "(100 x) Residual covariance matrix (SNPs)")

plot.cg_(ginf3$g[[2]], main = "Genes (Painting) Graph")
create_chart(mydatap2[[2]], predcov3[[2]],I, TRUE, "(100 x) Residual covariance matrix (Painting)")

plot.cg_(ginf3$g[[3]], main = "Language (Lex) Graph")
create_chart(mydatap2[[3]], predcov3[[3]],I, TRUE, "(100 x) Residual covariance matrix (Lex)")

plot.cg_(ginf3$g[[4]], main = "Language (Phonetics) Graph")
create_chart(mydatap2[[4]], predcov3[[4]],I, TRUE, "(100 x) Residual covariance matrix (Phonetics)")



# Plot predicted covariance matrices

Clarity_Chart(predcov3[[1]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (SNPs)")

Clarity_Chart(predcov3[[2]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Painting)")

Clarity_Chart(predcov3[[3]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Lexical)")

Clarity_Chart(predcov3[[4]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Phonetics)")



#### Inference with 4 admixture edges ########################################


# n_mix = 4

tg4 = add_mixedges(g = ginf3$g, n_edges = 1)
plot.cg_(tg4[[1]])

ginf4 = infer_graph(g = tg4, data = mydatap2, maxiter = 700, losstol = 1e-30,verbose = T,plotprogress = T)

predcov4 = ccov_dag_(ginf4$g)

par(mfrow=c(2,2))
plot.cg_(ginf4$g[[1]], main = "Genes (SNPs) Graph")
plot.cg_(ginf4$g[[2]], main = "Genes (Painting) Graph")
plot.cg_(ginf4$g[[3]], main = "Language (Lex) Graph")
plot.cg_(ginf4$g[[4]], main = "Language (Phonetics) Graph")


par(mfrow=c(1,2))

plot.cg_(ginf4$g[[1]], main = "Genes (SNPs) Graph")
create_chart(mydatap2[[1]], predcov4[[1]],I, TRUE, "(100 x) Residual covariance matrix (SNPs)")

plot.cg_(ginf4$g[[2]], main = "Genes (Painting) Graph")
create_chart(mydatap2[[2]], predcov4[[2]],I, TRUE, "(100 x) Residual covariance matrix (Painting)")

plot.cg_(ginf4$g[[3]], main = "Language (Lex) Graph")
create_chart(mydatap2[[3]], predcov4[[3]],I, TRUE, "(100 x) Residual covariance matrix (Lex)")

plot.cg_(ginf4$g[[4]], main = "Language (Phonetics) Graph")
create_chart(mydatap2[[4]], predcov4[[4]],I, TRUE, "(100 x) Residual covariance matrix (Phonetics)")


# Plot predicted covariance matrices
Clarity_Chart(predcov4[[1]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (SNPs)")

Clarity_Chart(predcov4[[2]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Painting)")

Clarity_Chart(predcov4[[3]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Lexical)")

Clarity_Chart(predcov4[[4]],scalefun=myscalefun2,text=T,main="Predicted covariance matrix (Phonetics)")

#########################################################################

# Plot real data matrices 

pt=plot.cg_(ginf0$g[[1]],show = FALSE)

par(mfrow=c(1,1))
Clarity_Chart(mydatap2[[1]],scalefun=myscalefun2,text=T,main="Genes (SNPs)")

Clarity_Chart(mydatap2[[2]],scalefun=myscalefun2,text=T,main="Genes (Painting)")

Clarity_Chart(mydatap2[[3]],scalefun=myscalefun2,text=T,main="Language (Lex)")

Clarity_Chart(mydatap2[[4]],scalefun=myscalefun2,text=T,main="Language (Phonetics)")




#######################################################################################

# ## SNP inference ##
# 
# ## Simulate a 9 population model where there is an outgroup - EAS_CHB
# tg0=simCoal_(9 ,labels=c("EAS_CHB","Mam","Tok","Bun","Kem","Wai","Mak","Fat","Tet"), outgroup="EAS_CHB")
# plot.cg_(tg0)
# 
# ginf0_snp = infer_graph(g = tg0, data = mydatap2[[1]], maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# 
# ## Phonetics inference ##
# 
# ginf0_phon = infer_graph(g = tg0, data = mydatap2[[4]], maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# ## Joint SNP and phonetics inference ##
# 
# ginfrlist0 = list(tg0, tg0)
# class(ginfrlist0) = "cglist"
# 
# snp_phon_data = list(mydatap2[[1]], mydatap2[[4]])
# 
# ginf0_snp_phon = infer_graph(g = ginfrlist0, data = snp_phon_data, maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# check_topology(ginf0_snp$g, ginf0_phon$g)
# check_topology(ginf0_phon$g, ginf0_snp_phon$g[[1]])
# 
# 
# ## 1 admixture ##
# 
# 
# tg1 = add_mixedges(tg0, n_edges = 1)
# ginfrlist1 = list(tg1, tg1)
# class(ginfrlist1) = "cglist"
# 
# ginf1_snp = infer_graph(g = tg1, data = mydatap2[[1]], maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# ginf1_phon = infer_graph(g = tg1, data = mydatap2[[4]], maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# ginf1_snp_phon = infer_graph(g = ginfrlist1, data = snp_phon_data, maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# plot.cg_(ginf1_snp_phon$g[[1]])
# 
# LRT = 2*(ginf1_snp_phon$loss - (ginf1_snp$loss+ ginf1_phon$loss) )
# 
# pchisq(q=LRT, df = 18, lower.tail = FALSE)
# 
# 
# 
# ## 2 admixtures ##
# 
# 
# tg2 = add_mixedges(tg1, n_edges = 1)
# ginfrlist2 = list(tg2, tg2)
# class(ginfrlist2) = "cglist"
# 
# ginf2_snp = infer_graph(g = tg2, data = mydatap2[[1]], maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# ginf2_phon = infer_graph(g = ginf2_snp$g, data = mydatap2[[4]], maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 
# ginf2_snp_phon = infer_graph(g = ginfrlist2, data = snp_phon_data, maxiter = 500, losstol = 1e-30, MHtol = 0,verbose = T, plotprogress = T)
# 



################################################################################

### Inference example ###

set.seed(12)
par(mfrow=c(1,2))


## Simulate a 4 population model where there is an outgroup
tg0=simCoal_(5,labels=c("A","B","C","D", "O"), outgroup="O")
plot.cg_(tg0)

tg1=simCoal_(5,labels=c("A","B","C","D", "O"), outgroup="O")
plot.cg_(tg1)

tg1 = reparameterise(tg0)
plot.cg_(tg1)


# tg2=mixedge_(tg0,2,3,0.5,0.2) ## Firstly with weight 0.2
# tg3=mixedge_(tg1,2,3,0.5,0.5) ## Seconldly with weight 0.5
# 
# tg4=mixedge_(tg2,1,4,0.5,0.3) ## Firstly with weight 0.3
# tg5=mixedge_(tg3,1,4,0.5,0.4) ## Seconldly with weight 0.4

# cglist1 = list(tg4, tg5)
# class(cglist1) = "cglist"
# dataref = lapply(cglist1, ccov_dag_)

cglist0 = list(tg0, tg1)
class(cglist0) = "cglist"

## Add "mixture edges"
cglist1 = add_mixedges(cglist0, n_edges=2, randomlistweights = T)
par(mfrow=c(1,2))
lapply(cglist1, plot.cg_)
# Simulate covariance matrices from these two models
dataref = lapply(cglist1, ccov_dag_)


# Simulate proposal graph 
trand0=simCoal_(5,labels=c("A","B","C","D","O"), outgroup = "O")
trand1 = add_mixedges(trand0, n_edges = 2)
plot.cg_(trand1)


## And of the two random graphs used for inference
trandlist=list(trand1,trand1)
class(trandlist)="cglist"

tinf1=infer_graph(trandlist, dataref, maxiter = 200, movefreqs = c(1/4,1/4,1/4,1/4), 
                  verbose = T, plotprogress = T, losstol = 1e-10, MHtol = 0)  # initial_temp = 1,cooling_rate = 0.99


pt=plot.cg_(tinf1$g[[1]],show = FALSE)
tpred1=ccov_dag_(tinf1$g)

Clarity_Chart(ordermatrix(tpred1[[2]],pt$order),scalefun=myscalefun2,text=T,main="Observed 1")


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

