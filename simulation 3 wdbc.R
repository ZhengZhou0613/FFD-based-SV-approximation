#simulation 3 sensitivity analysis for the WDBC dataset (Example 2)===========================================================
library(parallel)
library(rpart)
library(dplyr)
library(here)
source(here("shapley functions.R"))
set.seed(3)
#Read and preprocess dataset===========================================================================================
wdbc <- read.csv(here("data","wdbc.csv"),header=F)
wdbc[which(wdbc[,1]=="M"),1]<-0
wdbc[which(wdbc[,1]=="B"),1]<-1
wdbc[,1]<-as.numeric(wdbc[,1])
testindex<-sample(569,114)
traindata<-wdbc[-testindex,]
testdata<-wdbc[testindex,]
traindata$V1<-as.factor(traindata$V1)
d<-dim(traindata)[2]-1
#The modelling method (classification tree)================================================================================
rferr<-function(nowdata){
  itermodel<-1
  result<-rep(0,itermodel)
  for (i in 1:itermodel) {
    a<-rpart(V1~.,data=nowdata,method="class")
    predictresult<-rep(0,114)
    predictresult[which(predict(a,testdata)[,2]>0.5)]<-1
    result[i]<-length(which(predictresult==testdata[,1]))/114
  }
  return(mean(result))
}
#Define the value function of sensitivity analysis, which is provided by Cohen et al. (2005)====================================================================================================
#Unselected features are neutralized by assigning them a constant value across the training data, making them inactive during classification. ==================================
xt<-traindata[1,]
nulldata<-traindata
nulldata[,-1]<-xt[,-1]
val0<-rferr(nulldata)

length(which(testdata[,1]==1))/114
val<-function(sets,xt){
  xnew<-traindata
  xnew[,-c(1,1+sets)]<-xt[-c(1,1+sets)]
  return(rferr(xnew)-val0)
}

#Set the number of parallel cores, d%%corenumber==0 =================================================================================================
corenumber<-33
iter<-100

#RD (Apply recursive designs to approximate Shapley values, Algorithm 1)===========================================
phirec<-matrix(0,iter,d) #The matrix to record the approximations of RD=========================================
nrec<-rep(0,iter)
for (r in 1:iter) {
  time1<-Sys.time()
  #Generate recursive designs=============================================================================================================================
  d<-30
  drecur<-drec(d)
  permutation<-sample(2:drecur)
  design0<-recurdesign(drecur)[,c(1,permutation)][,1:d]
  design<-design0
  for(j in 2:d){
    designnow<-design0
    a<-designnow[,1]
    designnow[,(2:j)-1]<-designnow[,2:j]
    designnow[,j]<-a
    design<-rbind(design,designnow)
  }
  reprowtag<-left_join(data.frame(design),mutate(distinct(data.frame(design)), rn=row_number()))$rn
  design<-design[!duplicated(design),]
  n<-nrow(design)
  nrec[r]<-n
  #Evaluate the coalition values in parallel=====================================================================================================================
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))  
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(lhs),library(rpart),set.seed(NULL)))
  clusterExport(cl,c("permutation","corenumber","recurdesign","rferr","traindata","testdata","ycombine","val","xt","val0","drecur","design","d","rowgroup"),envir = environment())
  resultrec<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-rep(0,length(rowgroup))
    for (i in 1:length(rowgroup)) {
      result[i]<-val(which(design[rowgroup[i],]>0),xt)
    }
    result
  })
  stopCluster(cl)
  a<-resultrec[[1]]
  for (i in 2:corenumber) {
    a<-c(a,resultrec[[i]])
  }
  #Derive the approximations provided by Algorithm 1==================================================================================================
  for(i in 1:d){
    phirec[r,i]<-sum(ycombine(drecur)*a[reprowtag[((i-1)*4*(drecur-1)+1):(i*4*(drecur-1))]])
  }
  time2<-Sys.time()
  print(time2-time1)
  print(phirec[r,])
}

#Bias detection and correction (shrinkage)========================================================================================================================
phirec1star<-matrix(0,iter,d) #The matrix to record the corrected approximations of RD (shrinkage)===========================================================
for(i in 1:iter){
  phirec1star[i,]<-phirec[i,]-(sum(phirec[i,])-val(c(1:d),xt))/sum(phirec[i,])*phirec[i,]
}

#Bias detection and correction (Algorithm 2)========================================================================================================================
nstar<-rep(0,iter)
phirecstar<-matrix(0,iter,d) #The matrix to record the corrected approximations of RD===================================================================================
philsadd<-matrix(0,iter,d)
testresult<-rep(T,d)
#Numbers of key players selected in the two-step correction ================================================================================================================
t_1<-6
t_2<-21
for (r in 1:iter) {
  time1<-Sys.time()
  print(r)
  # Test whether $\widetilde{\bm{\Phi}}$ and $\bm{\Phi}$ are proportional=====================================================================================================
  phirecorder<-order(phirec[r,])
  ud<-round(seq(from=1,to=d,by=(d-1)/(t_1-1)))
  keyplayer1<-c()
  for (j in ud) {
    keyplayer1<-c(keyplayer1,phirecorder[j])
  }
  design<-rbind(RLHD(d),RLHD(d),RLHD(d),RLHD(d))
  n<-nrow(design)
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(rpart),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup","keyplayer1","rferr","traindata","testdata","val0"),envir = environment())
  resultlsadd<-parLapply(cl,1:corenumber,function(k){
    n<-0
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-matrix(0,length(rowgroup),d)
    for(i in 1:length(rowgroup)){
      row<-design[rowgroup[i],]
      prev<-0
      for (j in 1:d) {
        if(row[j]%in%keyplayer1){
          prevnext<-val(row[1:j],xt)
          if(j>=2){
            if(row[j-1]%in%keyplayer1){
              result[i,row[j]]<-result[i,row[j]]+prevnext-prev
              n<-n+1
            }else{
              result[i,row[j]]<-result[i,row[j]]+prevnext-val(row[1:(j-1)],xt)
              n<-n+2
            }
          }else{
            result[i,row[j]]<-result[i,row[j]]+prevnext-0
            n<-n+1
          }
          prev<-prevnext
        }
      }
    }
    list(result,n)
  })
  stopCluster(cl)
  a<-resultlsadd[[1]][[1]]
  nadd<-resultlsadd[[1]][[2]]
  for (i in 2:corenumber) {
    a<-rbind(a,resultlsadd[[i]][[1]])
    nadd<-nadd+resultlsadd[[i]][[2]]
  }
  philsadd[r,keyplayer1]<-colMeans(a[,keyplayer1])
  nstar[r]<-nstar[r]+nadd
  lmresult<-summary(lm(phirec[r,keyplayer1]~philsadd[r,keyplayer1]))
  #If the test passes, correct the bias by the shrinkage method ========================================================================
  if((coef(lmresult)[1,4]>0.05)&(coef(lmresult)[2,4]<0.05)&(lmresult$r.squared>0.75)){
    testresult[r]<-T
    print("linearity test=")
    print(testresult[r])
    phirecstar[r,]<-phirec[r,]-(sum(phirec[r,])-val(c(1:d),xt))/sum(phirec[r,])*phirec[r,]
    print(phirecstar[r,])
    time2<-Sys.time()
    print(time2-time1)
  }else{
    testresult[r]<-F
    print("linearity test=")
    print(testresult[r])
  }
  # If the test fails, we continue to add evaluation coalitions to further explore the potential relationship between $\widetilde{\bm{\Phi}}$ and $\bm{\Phi}$.==================================================
  if(testresult[r]==F){
    ud<-round(seq(from=1,to=d,by=(d-1)/(t_2-1)))
    keyplayer2<-c()
    for (j in ud) {
      keyplayer2<-c(keyplayer2,phirecorder[j])
    }
    addkeyplayer<-setdiff(keyplayer2,keyplayer1)
    cl <- makeCluster(corenumber)
    clusterSetRNGStream(cl)
    clusterEvalQ(cl,c(library(rpart),set.seed(NULL)))
    clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup","addkeyplayer","rferr","traindata","testdata","val0"),envir = environment())
    resultlsadd<-parLapply(cl,1:corenumber,function(k){
      n<-0
      rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
      result<-matrix(0,length(rowgroup),d)
      for(i in 1:length(rowgroup)){
        row<-design[rowgroup[i],]
        prev<-0
        for (j in 1:d) {
          if(row[j]%in% addkeyplayer){
            prevnext<-val(row[1:j],xt)
            if(j>=2){
              if(row[j-1]%in%addkeyplayer){
                result[i,row[j]]<-result[i,row[j]]+prevnext-prev
                n<-n+1
              }else{
                result[i,row[j]]<-result[i,row[j]]+prevnext-val(row[1:(j-1)],xt)
                n<-n+2
              }
            }else{
              result[i,row[j]]<-result[i,row[j]]+prevnext-0
              n<-n+1
            }
            prev<-prevnext
          }
        }
      }
      list(result,n)
    })
    stopCluster(cl)
    a<-resultlsadd[[1]][[1]]
    nadd<-resultlsadd[[1]][[2]]
    for (i in 2:corenumber) {
      a<-rbind(a,resultlsadd[[i]][[1]])
      nadd<-nadd+resultlsadd[[i]][[2]]
    }
    philsadd[r,addkeyplayer]<-colMeans(a[,addkeyplayer])
    nstar[r]<-nstar[r]+nadd
    #Interpolation-based correction==================================================================================================
    model<-ss(phirec[r,keyplayer2],philsadd[r,keyplayer2],0.75)
    phirecstar[r,]<-predict(model,phirec[r,])$y
    print(phirecstar[r,])
    time2<-Sys.time()
    print(time2-time1)
  }
}

#Record the cost of the correction method==========================================================================================
nrecstar<-nrec+nstar


identical(phirec1star,phirecstar)

#identical(phirec1star,phirecstar)=TRUE, which shows that the shrinkage method is consistently selected by Algorithm 2 for bias detection and correction.

#LS========================================================================================================================
phils<-matrix(0,iter,d) #The matrix to record the approximations of LS===========================================================
nls<-rep(0,iter)
for (r in 1:iter) {
  design<-rbind(RLHD(d),RLHD(d),RLHD(d),RLHD(d))
  time1<-Sys.time()
  n<-nrow(design)
  nls[r]<-n*d
  rowindex<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(lhs),library(rpart),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","rferr","val","xt","traindata","testdata","d","rowindex","val0"),envir = environment())
  resultls<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-matrix(0,length(rowindex),d)
    for(i in 1:length(rowindex)){
      row<-design[rowindex[i],]
      prev<-0
      for (j in 1:d) {
        index<-row[1:j]
        prevnext<-val(index,xt)
        result[i,row[j]]<-result[i,row[j]]+prevnext-prev
        prev<-prevnext
      }
    }
    result
  })
  stopCluster(cl)
  time2<-Sys.time()
  a<-resultls[[1]]
  for (i in 2:corenumber) {
    a<-rbind(a,resultls[[i]])
  }
  phils[r,]<-colMeans(a)
  print(time2-time1)
  print(phils[r,])
}

#STRRS========================================================================================================================
phistrrs<-matrix(0,iter,d) #The matrix to record the approximations of StrRS==================================================
#runs should be a multiple of d
runs<-2*d
nstrrs<-rep(runs*2*d,iter)
for (r in 1:iter) {
  design<-strrsdesign(d,runs)
  n<-length(design)
  time1<-Sys.time()
  rowindex<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(lhs),library(rpart),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","rferr","val","xt","traindata","testdata","d","rowindex","val0"),envir = environment())
  resultstrrs<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-rep(0,length(rowindex))
    for(i in 1:length(rowindex)){
      row<-design[[rowindex[i]]]
      result[i]<-val(row,xt)-val(row[-length(row)],xt)
    }
    result
  })
  stopCluster(cl)
  time2<-Sys.time()
  a<-resultstrrs[[1]]
  for (i in 2:corenumber) {
    a<-c(a,resultstrrs[[i]])
  }
  for(i in 1:d){
    phistrrs[r,i]<-mean(a[((i-1)*runs+1):(i*runs)])
  }
  print(time2-time1)
  print(phistrrs[r,])
}

#SRS==============================================================================================================================
phirs<-matrix(0,iter,d) #The matrix to record the approximations of SRS===========================================================
nrs<-rep(0,iter)
for (r in 1:iter) {
  design<-rsdesign(d,ceiling(nrecstar[r]/d))
  time1<-Sys.time()
  n<-nrow(design)
  nrs[r]<-n*d
  rowindex<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(lhs),library(rpart),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","rferr","val","xt","traindata","testdata","d","rowindex","val0"),envir = environment())
  resultrs<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-matrix(0,length(rowindex),d)
    for(i in 1:length(rowindex)){
      row<-design[rowindex[i],]
      prev<-0
      for (j in 1:d) {
        index<-row[1:j]
        prevnext<-val(index,xt)
        result[i,row[j]]<-result[i,row[j]]+prevnext-prev
        prev<-prevnext
      }
    }
    result
  })
  stopCluster(cl)
  time2<-Sys.time()
  a<-resultrs[[1]]
  for (i in 2:corenumber) {
    a<-rbind(a,resultrs[[i]])
  }
  phirs[r,]<-colMeans(a)
  print(time2-time1)
  print(phirs[r,])
}

#True (Obtain highly accurate approximations of Shapley values through intensive computations, serving as substitutes for true Shapley values)========================================================================================================================
phirstrue<-matrix(0,iter,d)
for (r in 1:iter) {
  print(r)
  design<-rsdesign(d,d^2)
  time1<-Sys.time()
  n<-nrow(design)
  rowindex<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(lhs),library(rpart),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","rferr","val","xt","traindata","testdata","d","rowindex","val0"),envir = environment())
  resultrs<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-matrix(0,length(rowindex),d)
    for(i in 1:length(rowindex)){
      row<-design[rowindex[i],]
      prev<-0
      for (j in 1:d) {
        index<-row[1:j]
        prevnext<-val(index,xt)
        result[i,row[j]]<-result[i,row[j]]+prevnext-prev
        prev<-prevnext
      }
    }
    result
  })
  stopCluster(cl)
  time2<-Sys.time()
  a<-resultrs[[1]]
  for (i in 2:corenumber) {
    a<-rbind(a,resultrs[[i]])
  }
  phirstrue[r,]<-colMeans(a)
  print(time2-time1)
  print(phirstrue[r,])
}

true<-colMeans(phirstrue) # The substitutes for true Shapley values.

#Compute the squared loss of each method==============================================================================================================================
mserecstar<-rep(0,iter)
msels<-rep(0,iter)
msestrrs<-rep(0,iter)
msers<-rep(0,iter)

for(i in 1:iter){
  mserecstar[i]<-sum((true-phirecstar[i,])^2)
  msels[i]<-sum((true-phils[i,])^2)
  msestrrs[i]<-sum((true-phistrrs[i,])^2)
  msers[i]<-sum((true-phirs[i,])^2)
}

save.image(here("results","simulation 3.RData"))


load(here("results","simulation 3.RData"))

#Cost Comparision (Print \bm{r} and \bm{I})=========================================
r<-c(median(nrecstar),median(nls),median(nstrrs),median(nrs))
i<-c(r[2]/d^2,r[3]/d^2/2,r[4]/d)
print(r)
print(i)

#Plot (Figure 7 (a))=============================================================================================
library(ggplot2)
library(ggformula)
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 15, base_line_size = 2/lwd_pt, base_rect_size = 3/lwd_pt)+	theme(panel.grid.major.x=element_line(colour=NA),panel.grid.minor = element_blank())+theme(plot.title=element_text(hjust=0.5)))

mse4d<-data.frame(c(mserecstar,msels,msestrrs,msers),rep(c("RD*","LS","StrRS","SRS"),rep(100,4)))
colnames(mse4d)<-c("values","methods")
mse4d$methods<-factor(mse4d$methods,levels=c("RD*","LS","StrRS","SRS"),ordered = TRUE)
ggplot(mse4d, aes(x=methods, y=values,color=methods))+geom_boxplot()+labs(x="Methods",y="Square Loss")+scale_color_manual(values=c("#db6968","#4d97cd","#459943","#e8c559"))+theme(legend.position = 'none')+scale_y_continuous()+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text=element_text(size=22))
ggsave(here("figures","wdbc.pdf"),plot = last_plot(),width = 7,height = 6)

#Correct rate of identifying top 10 important features===============================================================
ifrrecstar<-rep(0,iter)
ifrrs<-rep(0,iter)
ifrstrrs<-rep(0,iter)
ifrls<-rep(0,iter)
iftrue<-order(true)[21:30]
for (i in 1:iter) {
  ifrrecstar[i]<-length(intersect(iftrue,order(phirecstar[i,])[21:30]))/10
  ifrls[i]<-length(intersect(iftrue,order(phils[i,])[21:30]))/10
  ifrstrrs[i]<-length(intersect(iftrue,order(phistrrs[i,])[21:30]))/10
  ifrrs[i]<-length(intersect(iftrue,order(phirs[i,])[21:30]))/10
}

mean(ifrrecstar)
mean(ifrls)
mean(ifrstrrs)
mean(ifrrs)