#########  ########             Feng PKM2    Mac    10-1-2017
setwd("/Volumes/EXT/Feng_xuhui/mRNA_seq")

library(Rsubread)
library(edgeR)



asd

asd




Exp=featureCounts(files=c("IonXpress_RNA_001.bam",
                          "IonXpress_RNA_002.bam",
                          "IonXpress_RNA_003.bam",
                          "IonXpress_RNA_004.bam",
                          "IonXpress_RNA_005.bam",
                          "IonXpress_RNA_006.bam",
                          "IonXpress_RNA_007.bam",
                          "IonXpress_RNA_008.bam",
                          "IonXpress_RNA_009.bam",
                          "IonXpress_RNA_010.bam",
                          "IonXpress_RNA_011.bam",
                          "IonXpress_RNA_012.bam" ),
                  annot.ext = "/Volumes/DATA/Mus_musculus_UCSC_mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf",
                  isGTFAnnotationFile = TRUE,
                  nthreads = 8)
#write.table(Exp$counts,"PKM2_mRNA",sep="\t",quote = FALSE,col.names =FALSE)
#write.table(cbind(Exp$counts,Exp$annotation$Length),"WT_PR_AHMR_length.txt",sep="\t",quote = FALSE,col.names =FALSE)
D=DGEList(counts=Exp$counts)
d=calcNormFactors(D)
colnames(d$counts)=c("WT1","WT2","WT3","PKM2KO1","PKM2KO2","PKM2KO3","WT_TAC1","WT_TAC2","WT_TAC3","PKM2KO_TAC1","PKM2KO_TAC2","PKM2KO_TAC3")
rpkm_exp=rpkm(d,gene.length = Exp$annotation$Length)
d$counts=rpkm_exp
#d$samples$group=c("WT1","WT2","WT3","PKM2KO1","PKM2KO2","PKM2KO3","WT_TAC1","WT_TAC2","WT_TAC3","PKM2KO_TAC1","PKM2KO_TAC2","PKM2KO_TAC3")
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)
write.table(d$pseudo.counts,"All_group_rpkm.txt",sep="\t",row.names=TRUE)











#################   windows   ##############
########   headmap   ###########

setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.table("All_group_rpkm_ori.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

keep=grep("^Mir",rownames(exp))
data_value2=exp[-keep,]
#colnames(data_value2)=c("Ctrl_Sham_1","Ctrl_Sham_2","Ctrl_Sham_3","Ctrl_TAC_1","Ctrl_TAC_2","Ctrl_TAC_3",
#                        "Mut-Sham-1","Mut-Sham-2","Mut-Sham-3","Mut-TAC-1","Mut-TAC-2","Mut-TAC-3")



#keep_a = rownames(data_value2) %in% c("Myl2","Mb","Myh6","Tnni3","Actc1","Tnnt2")
#data_value2=data_value2[!keep_a,]



exp_mtr=as.matrix(data_value2)
b=as.matrix(exp_mtr)

keep1=rowSums(b)>12
data_value3=b[keep1,]



data_value5=data_value3[,c(1:3,7:9,4:6,10:12)]
colnames(data_value5)=c("Ctrl_Sham_1","Ctrl_Sham_2","Ctrl_Sham_3","Ctrl_TAC_1","Ctrl_TAC_2","Ctrl_TAC_3",
                        "Mut-Sham-1","Mut-Sham-2","Mut-Sham-3","Mut-TAC-1","Mut-TAC-2","Mut-TAC-3")
library(gplots); 
library("RColorBrewer")
myheatcol <- bluered(100)

#col=rev(colorRampPalette(c("red", "yellow", "blue"))(256))

col = colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255)

hr <- hclust(as.dist(1-cor(t(data_value5), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(data_value5, method="spearman")), method="complete")
heatmap.2(data_value5, Rowv=as.dendrogram(hr), Colv=FALSE, col=col,keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))





heatmap.2(data_value3, Rowv=as.dendrogram(hr), Colv=FALSE, col=myheatcol,keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))



heatmap.2(data_value3, Rowv=as.dendrogram(hr), Colv=FALSE, col=rev(redblue(256)), keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))


heatmap.2(data_value5, Rowv=as.dendrogram(hr), Colv=FALSE, col=rev(redblue(256)), keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))




#####################  scatter_plot 10-1-2017  #####
setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.csv("All_group_rpkm.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

keep=grep("^Mir",rownames(exp))
data_value2=exp[-keep,]




data_value2$WT_mean=rowSums(data_value2[,1:3])/3
data_value2$KO_mean=rowSums(data_value2[,4:6])/3
data_value2$WT_TAC_mean=rowSums(data_value2[,7:9])/3
data_value2$KO_TAC_mean=rowSums(data_value2[,10:12])/3

data_value2$KO_TAC_vs_WT_TAC=(data_value2$KO_TAC_mean+0.1)/(data_value2$WT_TAC_mean+0.1)
data_value2$KO_WT_TAC_flag=cut(data_value2$KO_TAC_vs_WT_TAC,breaks=c(-Inf,0.5,2,Inf),labels=c("-1","0","1"))


library(ggplot2)
library(ggrepel)

p <- ggplot(data_value2, aes(log(WT_TAC_mean,2),log(KO_TAC_mean,2)))
pic=p+geom_point(aes(color=factor(KO_WT_TAC_flag)),alpha=0.8)+
  coord_cartesian(xlim=c(0,15),ylim=c(0,15))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)+
  xlab("WT TAC")+ylab("KO TAC")
pic+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) +  theme(panel.grid=element_blank())



## WT vs WT_TAC

data_value2$WT_vs_WT_TAC=(data_value2$WT_mean+0.1)/(data_value2$WT_TAC_mean+0.1)
data_value2$WT_WT_TAC_flag=cut(data_value2$WT_vs_WT_TAC,breaks=c(-Inf,0.5,2,Inf),labels=c("-1","0","1"))


p <- ggplot(data_value2, aes(log(WT_mean,2),log(WT_TAC_mean,2)))
pic=p+geom_point(aes(color=factor(WT_WT_TAC_flag)),alpha=0.8)+
  coord_cartesian(xlim=c(0,15),ylim=c(0,15))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)+
  xlab("WT")+ylab("WT TAC")
pic+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) +  theme(panel.grid=element_blank())



## WT vs KO_TAC

data_value2$WT_vs_KO_TAC=(data_value2$WT_mean+0.1)/(data_value2$KO_TAC_mean+0.1)
data_value2$WT_KO_TAC_flag=cut(data_value2$WT_vs_KO_TAC,breaks=c(-Inf,0.5,2,Inf),labels=c("-1","0","1"))


p <- ggplot(data_value2, aes(log(WT_mean,2),log(KO_TAC_mean,2)))
pic=p+geom_point(aes(color=factor(WT_KO_TAC_flag)),alpha=0.8)+
  coord_cartesian(xlim=c(0,15),ylim=c(0,15))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)+
  xlab("WT")+ylab("KO TAC")
pic+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) +  theme(panel.grid=element_blank())

#+
#  geom_text_repel(
#    data = SC_LE[SC_LE$Gene_name %in% c("Ajuba","Pttg1","Prom1","Ttbk2","Mc4r","Dpysl3","Enpp3"),],
#    aes(label = Gene_name),
#    size = 5,
#    box.padding = unit(1, "lines"), 
#    point.padding = unit(0.1, "lines"))


   






###############   PCA plot   9-7-2017   ##############

#install.packages('ggfortify')
#library(ggfortify)
#library(ggplot2)
#library(cluster)
#aa=SC_LE[,2:7]
#bb=t(aa)
#autoplot(prcomp(log(aa+1,2)))
#autoplot(bb)
#library(ggfortify)
#df <- iris[c(1, 2, 3, 4)]
#autoplot(prcomp(df))













######################    coorelationship efficience    10-1-2017   ###########
setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.csv("All_group_rpkm.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

keep=grep("^Mir",rownames(exp))
data_value2=exp[-keep,]



library(ggplot2)
library(ggrepel)

mcor = cor(data_value2,method = c("pearson"))
col=colorRampPalette(c("red", "yellow", "blue"))(200)
heatmap.2(mcor,col=col,density.info="none", trace="none", margins=c(12,12),symm = TRUE)



dists = dist( t(data_value2 ) )
mat = as.matrix( dists )
library("RColorBrewer")
hmcol=colorRampPalette(c("red", "yellow", "blue"))(200)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13), density.info="none")





p <- ggplot(SC_LE, aes(log(WT_LE3,2),log(KO_LE1,2)))
pic=p+geom_point()+
  coord_cartesian(xlim=c(0,15),ylim=c(0,15))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)+
  xlab("WT")+ylab("V2Ltf")
pic



###########   headmap for Energy metabolism   11-19-2017  ############

setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.table("All_group_rpkm_ori.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

keep=grep("^Mir",rownames(exp))
data_value2=exp[-keep,]




keep = rownames(data_value2) %in% c("Acot2","Mlycd","Sod2","Acacb","Acot1","Acaa2","Idh2","Dgat2","Acadm","Idh3b","Gnmt",
                                    "Nr1d1","Adra1b","Coq9","Ndufs6","Pink1","Slc25a13","Dhtkd1","Ogdhl","Ces1d","Phkg1",
                                    "Pdk2","Cox8b","Acot3","Pdk4")

sample_genes=data_value2[keep,]
b=as.matrix(sample_genes)

a=b[,c(1:3,7:9,4:6,10:12)]
colnames(a)=c("Ctrl_Sham_1","Ctrl_Sham_2","Ctrl_Sham_3","Ctrl_TAC_1","Ctrl_TAC_2","Ctrl_TAC_3",
                        "Mut-Sham-1","Mut-Sham-2","Mut-Sham-3","Mut-TAC-1","Mut-TAC-2","Mut-TAC-3")



library(gplots); 
library("RColorBrewer")

col = colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255)

hr <- hclust(as.dist(1-cor(t(a), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(a, method="spearman")), method="complete")
heatmap.2(a, Rowv=as.dendrogram(hr), Colv=FALSE, col=col,keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))



###########   headmap for TCA  11-19-2017  ############

setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.table("All_group_rpkm_ori.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

keep=grep("^Mir",rownames(exp))
data_value2=exp[-keep,]




keep = rownames(data_value2) %in% c("Acly","Pck2","Idh2","Suclg1","Sdhb","Pdhb","Sdhd","Aco2","Dlst","Suclg2","Sucla2",
                                    "Idh3b","Sdha","Fh1","Mdh1","Ogdh","Sdhc","Ogdhl","Dld","Pdha1","Dlat","Idh3a","Cs",
                                    "Mdh2","Idh3g","Pcx","Aco1","Idh1")

sample_genes=data_value2[keep,]
b=as.matrix(sample_genes)

a=b[,c(1:3,7:9,4:6,10:12)]
colnames(a)=c("Ctrl_Sham_1","Ctrl_Sham_2","Ctrl_Sham_3","Ctrl_TAC_1","Ctrl_TAC_2","Ctrl_TAC_3",
              "Mut-Sham-1","Mut-Sham-2","Mut-Sham-3","Mut-TAC-1","Mut-TAC-2","Mut-TAC-3")



library(gplots); 
library("RColorBrewer")

col = colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255)

hr <- hclust(as.dist(1-cor(t(a), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(a, method="spearman")), method="complete")
heatmap.2(a, Rowv=as.dendrogram(hr), Colv=FALSE, col=col,keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))



###########   headmap for Inflammation  11-19-2017  ############

setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.table("All_group_rpkm_ori.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

keep=grep("^Mir",rownames(exp))
data_value2=exp[-keep,]




keep = rownames(data_value2) %in% c("Nt5e","S100a9","Ccr1","Camk4","Nfkbid","Ccl12","Ncf1","Ccr7","Aoah","Cxcr6","Ccl5","Ccl8","Siglece",
                                    "Cxcl10","Btk","Ackr1","Apod","Casp4","Sdc1","Serpinf1","Ccr2","Ada","Pik3ap1","Tlr2","Nupr1","Fcgr1",
                                    "Fcer1g","Thbs1","Fn1","Slc7a2","Csf1r","Apoe","Sphk1","Pstpip1","Cd300a","Il17ra","Ly86","Cybb","Sbno2",
                                    "Tlr1","Ctss","Acp5","Ephb6","Ccr5","Cd14","Nfkbiz","Itgam","Cd180","Hck","C1qtnf3","Spn","Hmox1","Adam8",
                                    "Fcgr3","Serpina3n","Tnfaip8l2","Il34","Cd24a","Tril","Il6","Pla2g7","Tlr13","Il33","Tgm2","Fcgr2b","Clcf1",
                                    "Tlr6","Ptafr","Ager","Saa3","Pik3cg","Ccl3","Themis2","Cysltr1","Mefv","Tlr9","Rasgrp1","Selp","Cxcr3","Nlrp3",
                                    "Clec7a","Cx3cr1","Egfr","Tlr8","Lat","Cd44","Tlr7","Tnfrsf11a","Serpine1","Tnfrsf4","Il1rn","Hp","Aim2","Chil1",
                                    "Ighg2b","Ccl2","Il1b","Ccl4","Nlrc4","Ankrd42")

sample_genes=data_value2[keep,]
b=as.matrix(sample_genes)

a=b[,c(1:3,7:9,4:6,10:12)]
colnames(a)=c("Ctrl_Sham_1","Ctrl_Sham_2","Ctrl_Sham_3","Ctrl_TAC_1","Ctrl_TAC_2","Ctrl_TAC_3",
              "Mut-Sham-1","Mut-Sham-2","Mut-Sham-3","Mut-TAC-1","Mut-TAC-2","Mut-TAC-3")



library(gplots); 
library("RColorBrewer")
library(heatmap3)

col = colorRampPalette(rev(brewer.pal(11, "RdBu")) )(255)
heatmap.2(a, Rowv=as.dendrogram(hr), Colv=FALSE, col=col,keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))

a_s=scale(a,center=F,scale=T)

hr <- hclust(as.dist(1-cor(t(a), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(a, method="spearman")), method="complete")

heatmap.2(a_s, trace="none", col = col, Colv=FALSE, margin=c(1, 13), density.info="none",keysize = 1.1)

heatmap.3(a, trace="none", scale="row", zlim=c(-3,3), reorder=FALSE,
          distfun=distCor, hclustfun=hclustAvg, col=col, symbreak=FALSE) 


heatmap.2(a, Rowv=as.dendrogram(hr), Colv=FALSE, col=col,keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8),breaks=c(-1,1))



#c=log2(a)
#c[is.infinite(c)]=0
#c_mean=(max(c)+min(c))/2
#c_s=c-c_mean
#library(gplots); 
#library("RColorBrewer")
#hmcol=colorRampPalette(c("red", "yellow", "blue"))(200)
#heatmap.2(c_s, trace="none", col = col, Colv=FALSE, margin=c(1, 13), density.info="none",keysize = 1.1)









###########   headmap of MI and MII genes   ############

setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.csv("All_group_rpkm.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)


keep = rownames(exp) %in% c("Cd86","Cd80","Tlr2","Nos2","Ccr7","Cd163","Mrc1","Retnla","Arg1")

sample_genes=exp[keep,]
b=as.matrix(sample_genes)

a=b[,-c(4:6)]
a=a[c(4,5,8,2,3,9,7,6,1),]


myheatcol <- bluered(100)
hr <- hclust(as.dist(1-cor(t(a), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(a, method="spearman")), method="complete")
heatmap.2(a, Rowv=FALSE, Colv=FALSE, col=hmcol,keysize = 1.2,
          scale="row", density.info="none", trace="none",margins=c(10,8))


#library("RColorBrewer")
#hmcol=rev(colorRampPalette(c("red", "white", "blue"))(200))
#heatmap.2(b, trace="none", col = rev(hmcol), margin=c(13, 13),Rowv=FALSE,Colv=FALSE, density.info="none")











































################         old code     2016           #########################
D=DGEList(counts=Exp$counts[,c(1:3,7:9)])
d=calcNormFactors(D)
colnames(d$counts)=c("WT1","WT2","WT3","WT_TAC1","WT_TAC2","WT_TAC3")
rpkm_exp=rpkm(d,gene.length = Exp$annotation$Length)
d$counts=rpkm_exp

d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)
d$samples$group=c("1","1","1","2","2","2")
et=exactTest(d)
p_val=cbind(2^(et$table$logFC),et$table$PValue)
colnames(p_val) = c("Fold","pValue")
data_value=cbind(d$pseudo.counts,p_val)
write.table(data_value,"Six_group_with_pvalue",sep="\t")






#write.table(d$pseudo.counts,"All_group_rpkm.txt",sep="\t",row.names=TRUE)


#####   WT and WT_TAC
D=DGEList(counts=Exp$counts)
d=calcNormFactors(D)
colnames(d$counts)=c("WT1","WT2","WT3","PKM2KO1","PKM2KO2","PKM2KO3","WT_TAC1","WT_TAC2","WT_TAC3","PKM2KO_TAC1","PKM2KO_TAC2","PKM2KO_TAC3")
rpkm_exp=rpkm(d,gene.length = Exp$annotation$Length)
d$counts=rpkm_exp
#d$samples$group=c("WT1","WT2","WT3","PKM2KO1","PKM2KO2","PKM2KO3","WT_TAC1","WT_TAC2","WT_TAC3","PKM2KO_TAC1","PKM2KO_TAC2","PKM2KO_TAC3")
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)




#et=exactTest(d)
#p_val=cbind(et$table$logFC,et$table$PValue)
#colnames(p_val) = c("log_Fold","pValue")
#data_value=cbind(d$counts,d$pseudo.counts,a)
#write.table(data_value,"D6_WT_PR_AM_all_rpkm.txt",sep="\t")





######### GO_Gene_list 2-2-2018    ######################


setwd("F:/New_Genome_data/Feng xuhui/mRNA_seq")
exp=read.table("All_group_rpkm_ori.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

keep=grep("^Mir",rownames(exp))
data_value2=exp[-keep,]
data_value2$symbol=row.names(data_value2)

target_gene=read.table("F:/New_Genome_data/Feng xuhui/mRNA_seq/pics/GO_gene_list.txt",sep="\t",header=FALSE,stringsAsFactors = FALSE)
colnames(target_gene)=c("symbol")


#small_list=merge(target_gene, data_value2,by="symbol",all.x = TRUE)


a=data_value2[target_gene$symbol,]
a$fold_change=rowSums(a[,c(10:12)])/rowSums(a[,c(7:9)])
write.table(a,"GO_target_list.txt",sep="\t",quote=FALSE,row.names = FALSE)


