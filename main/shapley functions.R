#Construct two-levle full factorial designs wit p factors
fulldesign<-function(p){
  n<-2^p
  design<-matrix(0,n,p)
  for (i in 1:n) {
    i_0<-i-1
    for (j in 1:p) {
      x<-floor(i_0/2^(p-j))
      if(x==1){
        design[i,j]<-1
      }else{
        design[i,j]<--1
      }
      i_0<-i_0%%2^(p-j)
    }
  }
  return(design)
}

#Construct the regression matrix of the aimdesign under model (3)
regdesign<-function(aimdesign){
  N<-dim(aimdesign)[1]
  p<-dim(aimdesign)[2]
  colintername<-list()
  for (i in 1:(p+1)) {
    colintername[[i]]<-i-1
  }
  r<-p+1
  for (i in 2:p) {
    print(i)
    combine<-combn(p,i)
    for(j in 1:dim(combine)[2]){
      print(j)
      colindex<-as.vector(combine[,j])
      r<-r+1
      colintername[[r]]<-colindex
      col<-rep(1,N)
      for (k in colindex){
        col<-col*aimdesign[,k]
      }
      aimdesign<-cbind(aimdesign,col)
    }
  }
  regdesign<-cbind(rep(1,N),aimdesign)
  colnames(regdesign)<-NULL
  return(list(design=regdesign,colname=colintername))
}

#Construct the recursive design with p factors, p=2^h+1
recurdesign<-function(p){
  h<-log2(p-1)
  H<-matrix(1,1,1)
  for (i in 1:h) {
    H<-cbind(rbind(H,H),rbind(H,-H))
  }
  H<-cbind(c(rep(1,2^(h+1)),rep(-1,2^(h+1))),rbind(H,-H,-H,H))
  return(H)
}

drec<-function(d){
  for(i in 2:15){
    if(d<=(2^i+1)){
      break
    }
  }
  return((2^i+1))
}

#Return the l defined in Theorem 3
ycombine<-function(d){
  n<-4*(d-1)
  l<-c((d+1)/(3*n),rep(2/(3*n),(d-2)))
  l<-c(l,l,-l,-l)
  return(2*l)
}

#Randomly generate an LHD of order d
RLHD<-function(d){
  design<-matrix(0,d,d)
  design[1,]<-1:d
  for (i in 2:d) {
    design[i,]<-c(design[i-1,2:d],design[i-1,1])
  }
  permutation<-sample(1:d)
  design<-design[,permutation]
  permutation<-sample(1:d)
  design<-design[permutation,]
  return(design)
}

#Randomly generate a design, each row of the design is a permutation of 1:d
rsdesign<-function(d,runs){
  design<-matrix(0,runs,d)
  for (i in 1:runs) {
    design[i,]<-sample(1:d)
  }
  return(design)
}

swap<-function(a,b,x){
  save<-x[a]
  x[a]<-x[b]
  x[b]<-save
  return(x)
}

#Randomly generate a design used for StrRS, number of runs should be a multiple of 2*d^2

strrsdesign<-function(d,runs){
  result<-vector("list",0)
  iter<-runs/d
  indesign<-vector("list",runs)
  for (i in 1:runs) {
    indesign[[i]]<-sample(1:d)
  }
  for (i in 1:d) {
    indesignnow<-indesign
    for (j in 1:runs) {
      mainpos<-floor((j-1)/iter)+1
      posnow<-which(indesignnow[[j]]==i)
      if(posnow!=mainpos){
        indesignnow[[j]]<-swap(mainpos,posnow,indesignnow[[j]])[1:mainpos]
      }else{
        indesignnow[[j]]<-indesignnow[[j]][1:mainpos]
      }
    }
    result<-c(result,indesignnow)
  }
  return(result)
}
