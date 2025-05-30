library(MASS)
library(mvtnorm)
library(parallel)
library(here)
source(here("main","shapley functions.R"))

# Compute the coefficients $B_i^j$ defined in Theorem 6, $i=2,...,8$, $j=2,...,2^i$. 
hmax<-8
B<-matrix(0,hmax,2^hmax)
for (i in 2:hmax) {
  n<-2^(i-1)
  for (j in 2:(2^hmax)) {
    if(j%%2==1){
      B[i,j]<-B[i-1,j]*2^(j-1)
      for (s in 1:min((j-1)/2,n)) {
        B[i,j]<-B[i,j]+B[i-1,j-2*s]*choose(n-j+2*s,s)*2^(j-2*s-1)
      }
    }
    if(j%%4==2){
      B[i,j]<-B[i-1,j]*2^(j-1)
      if(j>2){
        for (s in 1:min(j/2-1,n)) {
          B[i,j]<-B[i,j]+B[i-1,j-2*s]*choose(n-j+2*s,s)*2^(j-2*s-1)
        }
      }
    }
    if(j%%4==0){
      B[i,j]<-B[i-1,j]*2^(j-1)+choose(n,j/2)
      for (s in 1:min(j/2-1,n)) {
        B[i,j]<-B[i,j]+B[i-1,j-2*s]*choose(n-j+2*s,s)*2^(j-2*s-1)
      }
    }
  }
}
for (i in (1:128)*2-1) {
  B[,i]<-0
}

#Define $\psi^*(x)=e^{-x}$====================================================================================================
dist<-function(d){
  return(exp(-d))
}

#The function to compute  f(k,m) defined in Theorem 6==========================================================================
#Given $h$ the function retures a matrix $F=(f_{i,j})$ with $f_{i,j}=f(i,j)$===================================================
fm<-function(h){
  n<-2^h+1
  result<-matrix(0,n,n)
  for (k in 1:n) {
    for (m in 1:k) {
      indices <- (k+m-n):min(k, m)
      a <- sum(unlist(mclapply(indices, function(i) {
        choose(k-1,i-1) * choose(n-k,m-i) * dist(2*(k+m-2*i))
      })))
      a <- choose(n-1,k-1)*a
      result[k,m]<-a
      result[m,k]<-a
    }
  }
  return(result)
}

#The function to compute $g(k,m)$ defined in Theorem 6==========================================================================
#Given $h$ the function retures a matrix $G=(g_{i,j})$ with $g_{i,j}=g(i,j)$====================================================
gm<-function(h){
  n<-2^h+1
  result<-matrix(0,n,n)
  for (k in 1:n) {
    for (m in 1:k) {
      indices <- (k+m-n):min(k, m)
      a <- sum(unlist(mclapply(indices, function(i) {
        choose(k-1,i-1) * choose(n-k,m-i) * dist(2*(k+m-2*i))
      })))
      a <- choose(n-2,k-2)*a
      result[k,m]<-a
      result[m,k]<-a
    }
  }
  for (k in 1:(n-1)) {
    for (m in 1:k) {
      indices <- (k+m-n):min(k-1, m-1)
      a <- sum(unlist(mclapply(indices, function(i) {
        choose(k,i) * choose(n-k-1,m-i-1) * dist(2*(k+m-2*i))
      })))
      a <- choose(n-2,k-1)*a
      result[k,m]<-result[k,m]+a
      result[m,k]<-result[k,m]+a
    }
  }
  return(result)
}

#sigmaxx, sigmaxy and sigmayy are three functions to compute the terms in $Cov(\Phi_i,\Phi_j)$, $Cov(\Phi_i,\widetilde{\Phi}_j)$ and $Cov(\widetilde{\Phi}_i,\widetilde{\Phi}_j)$, respectively, which are defined in Theorem 6.====================================
sigmaxx<-function(h,f=fm(h),g=gm(h)){
  n<-2^h+1
  result<-matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:i) {
      if(i==j){
        covmatrix<-f
      }else{
        covmatrix<-g
      }
      if(n%%2==0){
        for (k in ((1:(n/2))*2-1)) {
          for (m in ((1:(n/2))*2-1)){
            result[i,j]<-result[i,j]+covmatrix[k,m]/(k*m)
          }
        }
      }else{
        for (k in ((1:((n+1)/2))*2-1)) {
          for (m in ((1:((n+1)/2))*2-1)){
            result[i,j]<-result[i,j]+covmatrix[k,m]/(k*m)
          }
        }
      }
      result[j,i]<-result[i,j]
    }
  }
  return(result)
}

sigmaxy<-function(h,f=fm(h),g=gm(h),C){
  n<-2^h+1
  result<-matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:i) {
      if(i==j){
        covmatrix<-f
      }else{
        covmatrix<-g
      }
      if(n%%2==0){
        for (k in ((1:(n/2))*2-1)) {
          for (m in ((1:(n/2))*2-1)){
            result[i,j]<-result[i,j]+covmatrix[k,m]*C[h,k]/m
          }
        }
      }else{
        for (k in ((1:((n+1)/2))*2-1)) {
          for (m in ((1:((n+1)/2))*2-1)){
            result[i,j]<-result[i,j]+covmatrix[k,m]*C[h,k]/m
          }
        }
      }
      result[j,i]<-result[i,j]
    }
  }
  return(result)
}

sigmayy<-function(h,f=fm(h),g=gm(h),C){
  n<-2^h+1
  result<-matrix(0,n,n)
  for (i in 1:n) {
    for (j in 1:i) {
      if(i==j){
        covmatrix<-f
      }else{
        covmatrix<-g
      }
      if(n%%2==0){
        for (k in ((1:(n/2))*2-1)) {
          for (m in ((1:(n/2))*2-1)){
            result[i,j]<-result[i,j]+covmatrix[k,m]*C[h,k]*C[h,m]
          }
        }
      }else{
        for (k in ((1:((n+1)/2))*2-1)) {
          for (m in ((1:((n+1)/2))*2-1)){
            result[i,j]<-result[i,j]+covmatrix[k,m]*C[h,k]*C[h,m]
          }
        }
      }
      result[j,i]<-result[i,j]
    }
  }
  return(result)
}

#Compute the Spearman correlation coefficient between $\widetilde{\bm{\Phi}}$ and $\bm{\Phi}$ for $n=2^h+1$, $h=2,...,8$.=============================================================================================
iter<-100
corresult<-matrix(0,8,iter)
for (h in 2:8) {
  print(h)
  C<-matrix(0,h,(2^h+1))
  C[,1]<-1
  C[,3]<-1/3
  for (j in 1:h) {
    for (i in ((2:(2^(h-1)))*2+1)) {
      if(i<(2^j+1)){
        C[j,i]<-2*B[j,i-1]/(choose(2^j,i-1)*3)+1/3
      }
    }
  }
  f<-fm(h)
  g<-gm(h)
  xx<-sigmaxx(h=h,f=f,g=g)
  xy<-sigmaxy(h=h,f=f,g=g,C)
  yy<-sigmayy(h=h,f=f,g=g,C)
  
  sigma<-rbind(cbind(xx,xy),cbind(t(xy),yy))
  samples <- rmvnorm(n = iter, mean = rep(0,dim(sigma)[2]), sigma = sigma)
  
  for (i in 1:dim(samples)[1]) {
    true<-samples[i,1:(dim(samples)[2]/2)]
    approx<-samples[i,(dim(samples)[2]/2+1):dim(samples)[2]]
    corresult[h,i]<-cor(x=approx,y=true,method = "spearman")
  }
}

save.image(here("results","Figure 8.RData"))

load(here("results","Figure 8.RData"))

#Plot (Figure 8)==========================================================================================================================================================================================================================================================
library(ggplot2)
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 15, base_line_size = 2/lwd_pt, base_rect_size = 3/lwd_pt)+	theme(panel.grid.major.x=element_line(colour=NA),panel.grid.minor = element_blank())+theme(plot.title=element_text(hjust=0.5)))

sccfig<-data.frame(rep(c("9","17","33","65","129","257"),100),as.vector(corresult[3:8,]))
colnames(sccfig)<-c("n","scc")
sccfig$n<-factor(sccfig$n,levels=c("9","17","33","65","129","257"),ordered = TRUE)
ggplot(sccfig, aes(x=n, y=scc, group=n))+geom_boxplot()+labs(x="Number of players",y="Spearman coefficients")+theme(legend.position = 'none')+scale_y_continuous()+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text=element_text(size=22))
ggsave(here("figures","sccfig.pdf"),plot = last_plot(),width = 8,height = 5)