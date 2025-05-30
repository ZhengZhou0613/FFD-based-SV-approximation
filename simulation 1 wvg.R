#simulation 1 weighted voting game (Example 3)=================================================
library(parallel)
library(dplyr)
library(npreg)
library(here)
source(here("shapley functions.R"))
set.seed(1)
#Initialize a weighted voting game==========================================================================================
#d=129, {w_1,...,w_d} randomly sampled from {1,\ldots,500}=============================================================================================================
xt<-sample(1:500,129)
d<-length(xt)
#The value function of the weighted voting game=============================================================================
val<-function(sets,xt){
  if(sum(xt[sets])>=(sum(xt)/2)){
    return(1) 
  }else{
    return(0)
  }
}

#Set the number of parallel cores, d%%corenumber==0 =================================================================================================
corenumber<-33
iter<-100
#RD (Apply recursive designs to approximate Shapley values, Algorithm 1)===============================================================================
phirec<-matrix(0,iter,d) #The matrix to record the approximations of RD=================================================================================
nrec<-rep(0,iter)
for (r in 1:iter) {
  print(r)
  time1<-Sys.time()
  #Generate recursive designs=============================================================================================================================
  drecur<-drec(d)
  permutation<-sample(2:drecur)
  design0<-recurdesign(drecur)[,c(1,permutation)][,1:d]
  design<-do.call(rbind, replicate(d,design0, simplify=FALSE))
  for(j in 2:d){
    index<-((j-1)*4*(drecur-1)+1):(j*4*(drecur-1))
    a<-design[index,1]
    design[index,(2:j)-1]<-design[index,2:j]
    design[index,j]<-a
  }
  reprowtag<-left_join(data.frame(design),mutate(distinct(data.frame(design)), rn=row_number()))$rn
  design<-design[!duplicated(design),]
  n<-nrow(design)
  nrec[r]<-n
  #Evaluate the coalition values in parallel=====================================================================================================================
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))  
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("permutation","corenumber","recurdesign","ycombine","val","xt","drecur","design","d","rowgroup"),envir = environment())
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

#Bias detection and correction (Algorithm 2)========================================================================================================================
nstar<-rep(0,iter)
phirecstar<-matrix(0,iter,d) #The matrix to record the corrected approximations of RD===================================================================================
philsadd<-matrix(0,iter,d)
testresult<-rep(T,iter)
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
  print(keyplayer1)
  design<-rbind(RLHD(d),RLHD(d),RLHD(d),RLHD(d),RLHD(d),RLHD(d))
  n<-nrow(design)
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(library(lhs),set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup","keyplayer1"),envir = environment())
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
    clusterEvalQ(cl,c(library(lhs),set.seed(NULL)))
    clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup","addkeyplayer"),envir = environment())
    resultlsadd<-parLapply(cl,1:corenumber,function(k){
      n<-0
      rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
      result<-matrix(0,length(rowgroup),d)
      for(i in 1:length(rowgroup)){
        row<-design[rowgroup[i],]
        prev<-0
        for (j in 1:d) {
          if(row[j]%in%addkeyplayer){
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
    which(testresult==F)
    #Interpolation-based correction==================================================================================================
    model<-ss(phirec[r,keyplayer2],philsadd[r,keyplayer2],spar=0.75)
    phirecstar[r,]<-predict(model,phirec[r,])$y
    plot(model)
    print(phirecstar[r,])
    time2<-Sys.time()
    print(time2-time1)
  }
}

#Record the cost of the correction method==========================================================================================
nrecstar<-nrec+nstar

#LS========================================================================================================================
phils<-matrix(0,iter,d) #The matrix to record the approximations of LS===========================================================
nls<-rep(0,iter)
for (r in 1:iter) {
  if(testresult[r]==T){
    design<-rbind(RLHD(d),RLHD(d),RLHD(d))
  }else{
    design<-rbind(RLHD(d),RLHD(d),RLHD(d),RLHD(d))
  }
  time1<-Sys.time()
  n<-nrow(design)
  nls[r]<-n*d
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup"),envir = environment())
  resultls<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-matrix(0,length(rowgroup),d)
    for(i in 1:length(rowgroup)){
      row<-design[rowgroup[i],]
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
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup"),envir = environment())
  resultstrrs<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-rep(0,length(rowgroup))
    for(i in 1:length(rowgroup)){
      row<-design[[rowgroup[i]]]
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
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup"),envir = environment())
  resultrs<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-matrix(0,length(rowgroup),d)
    for(i in 1:length(rowgroup)){
      row<-design[rowgroup[i],]
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
phitrue<-matrix(0,iter,d)
for(r in 1:iter){
  design<-rsdesign(d,800000)
  time1<-Sys.time()
  n<-nrow(design)
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","xt","d","rowgroup"),envir = environment())
  resulttrue<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-matrix(0,length(rowgroup),d)
    for(i in 1:length(rowgroup)){
      row<-design[rowgroup[i],]
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
  a<-resulttrue[[1]]
  for (i in 2:corenumber) {
    a<-rbind(a,resulttrue[[i]])
  }
  phitrue[r,]<-colMeans(a)
  print(time2-time1)
}
true<-colMeans(phitrue) # The substitutes for true Shapley values.

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

save.image(here("results","simulation 1.RData"))

load(here("results","simulation 1.RData"))

#Index1=the indices of replications outputting $\widetilde{\bm{\Phi}^*}$ as the final approximation======================================================================================
#Index2=the indices of replications outputting $\widetilde{\bm{\Phi}^**}$ as the final approximation=====================================================================================
index1<-which(testresult==T)
index2<-which(testresult==F)

#Cost Comparision (Print \bm{r} and \bm{I})=========================================
r<-c(median(nrecstar[index1]),median(nls[index1]),median(nstrrs[index1]),median(nrs[index1]))
i<-c(r[2]/d^2,r[3]/d^2/2,r[4]/d)
print(r)
print(i)

r<-c(median(nrecstar[index2]),median(nls[index2]),median(nstrrs[index2]),median(nrs[index2]))
i<-c(r[2]/d^2,r[3]/d^2/2,r[4]/d)
print(r)
print(i)

#Plots==========================================================================================================================================================
library(ggplot2)
library(ggformula)
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 15, base_line_size = 2/lwd_pt, base_rect_size = 3/lwd_pt)+	theme(panel.grid.major.x=element_line(colour=NA),panel.grid.minor = element_blank())+theme(plot.title=element_text(hjust=0.5)))

#Plot 1 (Fig 2 (a))==============================================================================================================================================
phirecorder<-order(phirec[10,])
ud<-round(seq(from=1,to=d,by=(d-1)/19))
keyplayer<-c()
for (j in ud) {
  keyplayer<-c(keyplayer,phirecorder[j])
}
relation<-matrix(0,d,3)
relation[,1]<-true
relation[,2]<-phirec[10,]
relation[,3]<-rep(0,d)
relation[keyplayer,3]<-1
relation<-data.frame(relation)
relation[,3]<-as.character(relation[,3])
colnames(relation)<-c("true","approximation","group")
ggplot(relation, aes(y=true,x=approximation,color=group,shape=group,linetype=group))+geom_point()+scale_shape_manual(values=c(19,8),guide="none")+scale_color_manual(values=c("black","red"),guide="none")+geom_smooth(method="lm",se=F)+geom_smooth(data = subset(relation, group>0),method = "lm",color = "red",se=F)+labs(y="True Shapley Values",x="Approximations obtained by Algorithm 1")+guides(linetype=F)
ggsave(here("figures","wvg1.pdf"),plot = last_plot(),width = 7,height = 6)

#Plot 2 (Fig 3 (a))===========================================================================================================================================================================
mse4d<-data.frame(c(mserecstar[index1],msels[index1],msestrrs[index1],msers[index1]),rep(c("RD*","LS","StrRS","SRS"),rep(length(index1),4)))
colnames(mse4d)<-c("values","methods")
mse4d$methods<-factor(mse4d$methods,levels=c("RD*","LS","StrRS","SRS"),ordered = TRUE)
ggplot(mse4d, aes(x=methods, y=values,color=methods))+geom_boxplot()+labs(x="Methods",y="Squared Loss")+scale_color_manual(values=c("#db6968","#4d97cd","#459943","#e8c559"))+theme(legend.position = 'none')+scale_y_continuous(limits = c(0,0.012))+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text=element_text(size=22))
ggsave(here("figures","wvg2.pdf"),plot = last_plot(),width = 7,height = 6)

#Plot 3 (Fig 3 (b))=============================================================================================================================================================================
mse4d<-data.frame(c(mserecstar[index2],msels[index2],msestrrs[index2],msers[index2]),rep(c("RD**","LS","StrRS","SRS"),rep(length(index2),4)))
colnames(mse4d)<-c("values","methods")
mse4d$methods<-factor(mse4d$methods,levels=c("RD**","LS","StrRS","SRS"),ordered = TRUE)
ggplot(mse4d, aes(x=methods, y=values,color=methods))+geom_boxplot()+labs(x="Methods",y="Squared Loss")+scale_color_manual(values=c("#db6968","#4d97cd","#459943","#e8c559"))+theme(legend.position = 'none')+scale_y_continuous(limits = c(0,0.012))+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text=element_text(size=22))
ggsave(here("figures","wvg3.pdf"),plot = last_plot(),width = 7,height = 6)


#Plot 4 (Fig 4)============================================================================================================================================================
phicompare<-matrix(0,d*length(index1)*5,3)
phicompare<-data.frame(phicompare)
phicompare[,1]<-as.character(rep(xt,length(index1)*5))
phicompare[,2]<-c(rep(true,length(index1)),as.vector(t(phirecstar[index1,])),as.vector(t(phils[index1,])),as.vector(t(phistrrs[index1,])),as.vector(t(phirs[index1,])))
phicompare[,3]<-rep(c("True","RD*","LS","StrRs","SRS"),each=d*length(index1))
colnames(phicompare)<-c("Player","Values","Method")
key<-sort(xt)[seq(from=5,to=125,by=10)]
phicompare<-phicompare[which(phicompare[,1]%in%key),]
phicompare$Player<-factor(phicompare$Player,levels=key)
phicompare$Method<-factor(phicompare$Method,levels=c("True","RD*","LS","StrRs","SRS"),ordered = T)
ggplot(phicompare, aes(x=Player,y=Values,color=Method))+geom_boxplot()+labs(x="Weights of Players",y="Values")+scale_color_manual(values=c("black","#db6968","#4d97cd","#459943","#e8c559"))+scale_y_continuous()+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text=element_text(size=20),legend.text = element_text(size=20))
ggsave(here("figures","wvg4.pdf"),plot = last_plot(),width = 12,height = 6)

