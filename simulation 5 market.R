#simulation 5 evaluating values of companies for the market dataset (Example 5)===========================================================
library(parallel)
library(dplyr)
library(here)
source(here("main","shapley functions.R"))
set.seed(5)
#Read the dataset=========================================================================================================================
W<-as.matrix(read.csv(here("data","crime.csv"),header = F))
d<-17
#Define the value function================================================================================================================
val<-function(sets){
  sum(diag(W)[sets])-sum(W[-sets,sets])
}

#Set the number of parallel cores, d%%corenumber==0 =================================================================================================
corenumber<-17
iter<-100
#RD (Apply recursive designs to approximate Shapley values, Algorithm 1)=============================================================================
#Since the RD method has been proved to provide the true Shapley values for this example, we do not repeat the RD method and designate its results as the true Shapley values.
true<-rep(0,d)
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
rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))  
cl <- makeCluster(corenumber)
clusterSetRNGStream(cl)
clusterEvalQ(cl,c(set.seed(NULL)))
clusterExport(cl,c("permutation","corenumber","recurdesign","ycombine","val","W","drecur","design","d","rowgroup"),envir = environment())
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
  true[i]<-sum(ycombine(drecur)*a[reprowtag[((i-1)*4*(drecur-1)+1):(i*4*(drecur-1))]])
}

#LS========================================================================================================================
phils<-matrix(0,iter,d) #The matrix to record the approximations of LS===========================================================
for (r in 1:iter) {
  design<-rbind(RLHD(d),RLHD(d))
  time1<-Sys.time()
  n<-nrow(design)
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","d","rowgroup","W"),envir = environment())
  resultls<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-matrix(0,length(rowgroup),d)
    for(i in 1:length(rowgroup)){
      row<-design[rowgroup[i],]
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
runs<-d
for (r in 1:iter) {
  design<-strrsdesign(d,runs)
  n<-length(design)
  time1<-Sys.time()
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","d","rowgroup","W"),envir = environment())
  resultstrrs<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-rep(0,length(rowgroup))
    for(i in 1:length(rowgroup)){
      row<-design[[rowgroup[i]]]
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
for (r in 1:iter) {
  design<-rsdesign(d,2*d)
  time1<-Sys.time()
  n<-nrow(design)
  rowgroup<-c(rep(n%/%corenumber+1,n%%corenumber),rep(n%/%corenumber,corenumber-n%%corenumber))
  cl <- makeCluster(corenumber)
  clusterSetRNGStream(cl)
  clusterEvalQ(cl,c(set.seed(NULL)))
  clusterExport(cl,c("corenumber","design","val","d","rowgroup","W"),envir = environment())
  resultrs<-parLapply(cl,1:corenumber,function(k){
    rowgroup<-(sum(rowgroup[0:(k-1)])+1):sum(rowgroup[0:k])
    result<-matrix(0,length(rowgroup),d)
    for(i in 1:length(rowgroup)){
      row<-design[rowgroup[i],]
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

save.image(here("results","simulation 5.RData"))

load(here("results","simulation 5.RData"))

#Compute the squared loss of each method==============================================================================================================================
msers<-rep(0,iter)
msestrrs<-rep(0,iter)
msels<-rep(0,iter)

for(i in 1:iter){
  msers[i]<-sum((true-phirs[i,])^2)
  msestrrs[i]<-sum((true-phistrrs[i,])^2)
  msels[i]<-sum((true-phils[i,])^2)
}

#Cost Comparision (Print \bm{r} and \bm{I})=========================================
r<-c(2*d^2-18,2*d^2,2*d^2,2*d^2)
i<-c(2,1,2*d)
print(r)
print(i)

#Plot (Figure 5 (b))===================================================================================================================
library(ggplot2)
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 15, base_line_size = 2/lwd_pt, base_rect_size = 3/lwd_pt)+	theme(panel.grid.major.x=element_line(colour=NA),panel.grid.minor = element_blank())+theme(plot.title=element_text(hjust=0.5)))

mse2d<-data.frame(c(rep(0,iter),msels,msestrrs,msers),rep(c("RD","LS","StrRS","SRS"),rep(iter,4)))
colnames(mse2d)<-c("values","methods")
mse2d$methods<-factor(mse2d$methods,levels=c("RD","LS","StrRS","SRS"),ordered = TRUE)
ggplot(mse2d, aes(x=methods, y=values,color=methods))+geom_boxplot()+labs(x="Methods",y="Square Loss")+scale_color_manual(values=c("#db6968","#4d97cd","#459943","#e8c559"))+theme(legend.position = 'none')+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text=element_text(size=22))
ggsave(here("figures","market.pdf"),plot = last_plot(),width = 7,height = 6)
