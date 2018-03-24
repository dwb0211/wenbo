setwd("F:/New_Genome_data/ZZA/bio_open")
RNA_seq=read.table("RNA-Seq.txt",sep="\t",  header=TRUE,stringsAsFactors = FALSE)
ma=read.table("ma1.txt",sep="\t",  header=TRUE,stringsAsFactors = FALSE)
colnames(RNA_seq)=c( "Gene.Symbol","End","Mes","Epi")

ma_a=ma[,c(1,2,3,8)]
ma_b=ma[,c(1,4,5,8)]
ma_c=ma[,c(1,6,7,8)]


#ma_a=ma_a[ma_a$A_call=="P",]
#ma_b=ma_b[ma_b$B_call=="P",]
#ma_c=ma_c[ma_c$C_call=="P",]

#########     End 2-2-2018   ############
ma_a_o=ma_c[order(ma_c$Gene.Symbol,decreasing = FALSE),]
library(dplyr)

gr=group_by(ma_a_o,Gene.Symbol)
Gr_index=group_indices(gr)  


for (i in length(Gr_index):1) {
  #print("aa")
  if(i!=length(Gr_index) && ma_a_o[i,4] == ma_a_o[i+1,4]){
    ma_a_o[i,2]=ma_a_o[i,2]+ma_a_o[i+1,2]
    # ma_a_o[i,2]=ma_a_o[i,2]+ma_a_o[i+1,2]
    #ma_a_o[i,1]=paste(ma_a_o[i+1,1],",",ma_a_o[i,1])
    #print(paste0("aa",i,"_",b[i-1,21],",",b[i,16])) 
  }else  {
    ma_a_o[i,2]=ma_a_o[i,2]
  }
}



ma_a_p=ma_a_o[!duplicated(ma_a_o$Gene.Symbol),]


a=merge(RNA_seq,ma_a_p,by="Gene.Symbol",all=TRUE)
a[is.na(a)]=0
#mcor = cor(a[,c(2,6)],method = c("pearson"))


a_s=a[order(a[,2],decreasing = FALSE),]
library(dplyr)
gr=group_by(a_s,Epi.x)
Gr_index=group_indices(gr)
a_s$RNA_rank=Gr_index

a_s2=a_s[order(a_s[,6],decreasing = FALSE),]
library(dplyr)
gr=group_by(a_s2,Epi.y)
Gr_index=group_indices(gr)
a_s2$microarray=Gr_index

write.table(a_s2,"Epi_list.txt",sep="\t",quote=FALSE,row.names = FALSE)




library(ggplot2)
p <- ggplot(a_s2,aes(a_s2[,8],a_s2[,9]))
pic=p+geom_point(colour="#0000ff") 
pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))   +  theme(panel.grid=element_blank())   # no grid


mcor = cor(a_s2[,c(8,9)],method = c("pearson"))
mcor




######



















































ma_a=ma_a[!duplicated(ma_a$Gene.Symbol),]

a=merge(RNA_seq,ma_a,by="Gene.Symbol",all=TRUE)
a[is.na(a)]=0
#mcor = cor(a[,c(2,6)],method = c("pearson"))


a_s=a[order(a[,2],decreasing = FALSE),]
library(dplyr)
gr=group_by(a_s,End.x)
Gr_index=group_indices(gr)
a_s$RNA_rank=Gr_index

a_s2=a_s[order(a_s[,6],decreasing = FALSE),]
library(dplyr)
gr=group_by(a_s2,End.y)
Gr_index=group_indices(gr)
a_s2$microarray=Gr_index




b=merge(RNA_seq,ma_b,by="Gene.Symbol",all=TRUE)
b[is.na(b)]=0
mcor = cor(b[,c(3,6)],method = c("pearson"))



c=merge(RNA_seq,ma_c,by="Gene.Symbol",all=TRUE)
c[is.na(c)]=0
mcor = cor(c[,c(4,6)],method = c("pearson"))



library(ggplot2)


p <- ggplot(a_s2,aes(a_s2[,8],a_s2[,9]))
pic=p+geom_point(aes(color="a"),alpha=1)
pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid


mcor = cor(a_s2[,c(8,9)],method = c("pearson"))

p <- ggplot(a_s2,aes(a_s2[,8],a_s2[,9]))
pic=p+geom_point(aes(color="a"),alpha=1)+coord_cartesian(xlim=c(-5,20000),ylim=c(-2,20000))
pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid


mcor = cor(a_s2[,c(8,9)],method = c("pearson"))













ma_a_o=ma_a[order(ma_a$Gene.Symbol,decreasing = FALSE),]
library(dplyr)

gr=group_by(ma_a_o,Gene.Symbol)
Gr_index=group_indices(gr)  


for (i in length(Gr_index):1) {
  #print("aa")
  if(i!=length(Gr_index) && ma_a_o[i,4] == ma_a_o[i+1,4]){
    ma_a_o[i,2]=(ma_a_o[i,2]+ma_a_o[i+1,2])/2
   # ma_a_o[i,2]=ma_a_o[i,2]+ma_a_o[i+1,2]
    #ma_a_o[i,1]=paste(ma_a_o[i+1,1],",",ma_a_o[i,1])
    #print(paste0("aa",i,"_",b[i-1,21],",",b[i,16])) 
  }else  {
    ma_a_o[i,2]=ma_a_o[i,2]
  }
}



ma_a_p=ma_a_o[!duplicated(ma_a_o$Gene.Symbol),]


a=merge(RNA_seq,ma_a_p,by="Gene.Symbol",all=TRUE)
a[is.na(a)]=0
#mcor = cor(a[,c(2,6)],method = c("pearson"))


a_s=a[order(a[,2],decreasing = FALSE),]
library(dplyr)
gr=group_by(a_s,End.x)
Gr_index=group_indices(gr)
a_s$RNA_rank=Gr_index

a_s2=a_s[order(a_s[,6],decreasing = FALSE),]
library(dplyr)
gr=group_by(a_s2,End.y)
Gr_index=group_indices(gr)
a_s2$microarray=Gr_index


library(ggplot2)
p <- ggplot(a_s2,aes(a_s2[,8],a_s2[,9]))
pic=p+geom_point(colour="#0000ff") 
pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid


mcor = cor(a_s2[,c(8,9)],method = c("pearson"))
