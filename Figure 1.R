library(here)
source(here("shapley functions.R"))

#The vector to record the coefficient n/d^2 ===========================================================================================================================================
recordco<-rep(0,257)

#Generate recursive designs for $d=5,...,257$, and compute $n/d^2$, where $n$ is the number of rows of the recursive design.===========================================================  
for (d in 5:257) {
  print(d)
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
  design<-design[!duplicated(design),]
  n<-nrow(design)
  recordco[d]<-n/d^2
  print(recordco[d])
}

rm(design)
rm(design0)

save.image(here("results","Figure 1.Rdata"))

load(here("results","Figure 1.Rdata"))

#Plot (Figure 1)================================================================================================================================================================
library(ggplot2)
lwd_pt <- .pt*72.27/96
theme_set(theme_bw(base_size = 15, base_line_size = 2/lwd_pt, base_rect_size = 3/lwd_pt)+	theme(panel.grid = element_blank())+theme(plot.title=element_text(hjust=0.5)))
plotdata<-matrix(0,253,2)
plotdata[,1]<-as.numeric(5:257)
plotdata[,2]<-recordco[5:257]
plotdata<-data.frame(plotdata)
colnames(plotdata)<-c("d","co")
ggplot(plotdata, aes(x=d, y=co))+geom_line(color="red",size=1)+geom_hline(aes(yintercept=4),linetype="dashed",colour="purple",size=1)+geom_hline(aes(yintercept=2),linetype="dashed",colour="blue",size=1)+scale_y_continuous(limits=c(1,5))+labs(y="r(C)/n^2",x="n")+theme(axis.title.x=element_text(size=22),axis.title.y=element_text(size=22),axis.text=element_text(size=22))
ggsave(here("figures","range.pdf"),plot = last_plot(),width = 12,height = 6)
