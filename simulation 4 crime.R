#simulation 4 local model explanation for the CC dataset (Example 6)===========================================================
library(parallel)
library(dplyr)
library(npreg)
library(randomForest)
library(here)
source(here("shapley functions.R"))
set.seed(4)
#Read and preprocess dataset=========================================================================================================
crime<-as.matrix(read.csv(here("data","crime.csv"),header=F))
crime<-apply(crime,2,as.numeric)
crime<-data.frame(na.omit(crime))
d<-dim(crime)[2]-1
n<-dim(crime)[1]
#Define the value function of the baseline Shapley for local model explanation, which is provided by Sundararajan and Najmi (2020)====================================================================================================
xt<-crime[1,]
xb<-colMeans(crime)
f<-randomForest(V100~.,data=crime,ntree=20)
val0<-predict(f,newdata=data.frame(t(xb)))
val<-function(sets){
  if(length(sets)>0){
    xnow<-xt
    xnow[-sets]<-xb[-sets]
    return(predict(f,xnow)-val0)
  }else{
    return(0)
  }
}

#Set the number of parallel cores, d%%corenumber==0 =================================================================================================
corenumber<-33
iter<-100
#RD (Apply recursive designs to approximate Shapley values, Algorithm 1)===========================================
phirec<-matrix(0,iter,d) #The matrix to record the approximations of RD=========================================
nrec<-rep(0,iter)
for (r in 1:iter) {
  time1<-Sys.time()
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
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))  
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(randomForest),set.seed(NULL)))
  clusterExport(cl,c("design","corenumber","recurdesign","ycombine","drec","d","rowgroup","val","val0","xt","xb","f"),envir = environment())
  resultrec<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-rep(0,length(rowgroup))
    for (i in 1:length(rowgroup)) {
      result[i]<-val(which(design[rowgroup[i],]>0))
    }
    result
  })
  stopCluster(cl)
  a<-resultrec[[1]]
  for (i in 2:corenumber) {
    a<-c(a,resultrec[[i]])
  }
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
  phirec1star[i,]<-phirec[i,]-(sum(phirec[i,])-val(c(1:d)))/sum(phirec[i,])*phirec[i,]
}

#Bias detection and correction (Algorithm 2)========================================================================================================================
nstar<-rep(0,iter)
phirecstar<-matrix(0,iter,d) #The matrix to record the corrected approximations of RD (Algorithm 2)=================================================================
philsadd<-matrix(0,iter,d)
testresult<-rep(T,d)
t_1<-6
t_2<-21
for (r in 1:iter) {
  time1<-Sys.time()
  print(r)
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
  clusterEvalQ(cl,c(library(randomForest),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","d","rowgroup","val","val0","xt","xb","f","keyplayer1"),envir = environment())
  resultlsadd<-parLapply(cl,1:corenumber,function(k){
    n<-0
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-matrix(0,length(rowgroup),d)
    for(i in 1:length(rowgroup)){
      row<-design[rowgroup[i],]
      prev<-0
      for (j in 1:d) {
        if(row[j]%in%keyplayer1){
          prevnext<-val(row[1:j])
          if(j>=2){
            if(row[j-1]%in%keyplayer1){
              result[i,row[j]]<-result[i,row[j]]+prevnext-prev
              n<-n+1
            }else{
              result[i,row[j]]<-result[i,row[j]]+prevnext-val(row[1:(j-1)])
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
  print(coef(lmresult)[,4])
  print(lmresult$r.squared)
  if((coef(lmresult)[1,4]>0.05)&(coef(lmresult)[2,4]<0.05)&(lmresult$r.squared>0.75)){
    testresult[r]<-T
    print("linearity test=")
    print(testresult[r])
    phirecstar[r,]<-phirec[r,]-(sum(phirec[r,])-val(c(1:d)))/sum(phirec[r,])*phirec[r,]
    print(phirecstar[r,])
    time2<-Sys.time()
    print(time2-time1)
  }else{
    testresult[r]<-F
    print("linearity test=")
    print(testresult[r])
  }
  if(testresult[r]==F){
    ud<-round(seq(from=1,to=d,by=(d-1)/(t_2-1)))
    keyplayer2<-c()
    for (j in ud) {
      keyplayer2<-c(keyplayer2,phirecorder[j])
    }
    length(keyplayer2)
    addkeyplayer<-setdiff(keyplayer2,keyplayer1)
    cl <- makeCluster(corenumber)
    clusterSetRNGStream(cl)
    clusterEvalQ(cl,c(library(randomForest),set.seed(NULL)))
    clusterExport(cl,c("corenumber","design","d","rowgroup","val","val0","xt","xb","f","addkeyplayer"),envir = environment())
    resultlsadd<-parLapply(cl,1:corenumber,function(k){
      n<-0
      rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
      result<-matrix(0,length(rowgroup),d)
      for(i in 1:length(rowgroup)){
        row<-design[rowgroup[i],]
        prev<-0
        for (j in 1:d) {
          if(row[j]%in% addkeyplayer){
            prevnext<-val(row[1:j])
            if(j>=2){
              if(row[j-1]%in%addkeyplayer){
                result[i,row[j]]<-result[i,row[j]]+prevnext-prev
                n<-n+1
              }else{
                result[i,row[j]]<-result[i,row[j]]+prevnext-val(row[1:(j-1)])
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
    model<-ss(phirec[r,keyplayer2],philsadd[r,keyplayer2])
    phirecstar[r,]<-predict(model,phirec[r,])$y
    print(phirecstar[r,])
    time2<-Sys.time()
    print(time2-time1)
  }
}

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
  clusterEvalQ(cl,c(library(randomForest),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","d","rowindex","val","val0","xt","xb","f"),envir = environment())
  resultls<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-matrix(0,length(rowindex),d)
    for(i in 1:length(rowindex)){
      row<-design[rowindex[i],]
      prev<-0
      for (j in 1:d) {
        index<-row[1:j]
        prevnext<-val(index)
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
  clusterEvalQ(cl,c(library(randomForest),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","d","rowindex","val","val0","xt","xb","f"),envir = environment())
  resultstrrs<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-rep(0,length(rowindex))
    for(i in 1:length(rowindex)){
      row<-design[[rowindex[i]]]
      result[i]<-val(row)-val(row[-length(row)])
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
  clusterEvalQ(cl,c(library(randomForest),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","d","rowindex","val","val0","xt","xb","f"),envir = environment())
  resultrs<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-matrix(0,length(rowindex),d)
    for(i in 1:length(rowindex)){
      row<-design[rowindex[i],]
      prev<-0
      for (j in 1:d) {
        index<-row[1:j]
        prevnext<-val(index)
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
  design<-rsdesign(d,20*d)
  time1<-Sys.time()
  n<-nrow(design)
  rowindex<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(randomForest),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","d","rowindex","val","val0","xt","xb","f"),envir = environment())
  resultrs<-parLapply(cl,1:corenumber,function(k){
    rowindex<-(sum(rowindex[0:(k-1)])+1):sum(rowindex[0:k])
    result<-matrix(0,length(rowindex),d)
    for(i in 1:length(rowindex)){
      row<-design[rowindex[i],]
      prev<-0
      for (j in 1:d) {
        index<-row[1:j]
        prevnext<-val(index)
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

save.image(here("results","simulation 4.RData"))

load(here("results","simulation 4.RData"))

#Cost Comparision (Print \bm{r} and \bm{I})=========================================
r<-c(median(nrecstar),median(nls),median(nstrrs),median(nrs))
i<-c(r[2]/d^2,r[3]/d^2/2,r[4]/d)
print(r)
print(i)

#Plot (Figure 7 (b))=============================================================================================
library(ggplot2)
library(ggformula)
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 15, base_line_size = 2/lwd_pt, base_rect_size = 3/lwd_pt)+	theme(panel.grid.major.x=element_line(colour=NA),panel.grid.minor = element_blank())+theme(plot.title=element_text(hjust=0.5)))

mse2d<-data.frame(c(mserecstar,msels,msestrrs,msers),rep(c("RD*","LS","StrRS","SRS"),rep(100,4)))
colnames(mse2d)<-c("values","methods")
mse2d$methods<-factor(mse2d$methods,levels=c("RD*","LS","StrRS","SRS"),ordered = TRUE)
ggplot(mse2d, aes(x=methods, y=values,color=methods))+geom_boxplot()+labs(x="Methods",y="Square Loss")+scale_color_manual(values=c("#db6968","#4d97cd","#459943","#e8c559"))+theme(legend.position = 'none')+scale_y_continuous()+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text=element_text(size=22))
ggsave(here("figures","crime.pdf"),plot = last_plot(),width = 7,height = 6)

#Correct rate of identifying top 10 important features===============================================================
ifrrecstar<-rep(0,iter)
ifrrs<-rep(0,iter)
ifrstrrs<-rep(0,iter)
ifrls<-rep(0,iter)

iftrue<-order(abs(true))[90:99]
for (i in 1:iter) {
  ifrrecstar[i]<-length(intersect(iftrue,order(abs(phirecstar[i,]))[90:99]))/10
  ifrls[i]<-length(intersect(iftrue,order(abs(phils[i,]))[90:99]))/10
  ifrstrrs[i]<-length(intersect(iftrue,order(abs(phistrrs[i,]))[90:99]))/10
  ifrrs[i]<-length(intersect(iftrue,order(abs(phirs[i,]))[90:99]))/10
}

mean(ifrrecstar)
mean(ifrls)
mean(ifrstrrs)
mean(ifrrs)