#########        edgeR for RPKM calculation

setwd("/Volumes/EXT/fenghua/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")

source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")
biocLite("edgeR")
library(Rsubread)
library(edgeR)

Exp=featureCounts(files=c("D6WT1.sam","D6WT2.sam","D6PR1.sam","D6PR2.sam","D6AM1.sam","D6AM2.sam"),
                  annot.ext = "/Volumes/DATA/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2013-03-06-15-01-24/Genes/genes.gtf",
                  isGTFAnnotationFile = TRUE,
                  nthreads = 8)
write.table(Exp$counts,"WT_PR_AHMR.txt",sep="\t",quote = FALSE,col.names =FALSE)
write.table(cbind(Exp$counts,Exp$annotation$Length),"WT_PR_AHMR_length.txt",sep="\t",quote = FALSE,col.names =FALSE)
D=DGEList(counts=Exp$counts)
d=calcNormFactors(D)
rpkm_exp=rpkm(d,gene.length = Exp$annotation$Length)
d$counts=rpkm_exp
d$samples$group=c("WT","WT","PR","PR","WT","WT")
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)
et=exactTest(d)
p_val=cbind(et$table$logFC,et$table$PValue)
colnames(p_val) = c("log_Fold","pValue")
data_value=cbind(d$counts,d$pseudo.counts,a)
write.table(data_value,"D6_WT_PR_AM_all_rpkm.txt",sep="\t")

#########        end

# -----------------------------       WT Vs PR  -----------------------
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_exp=read.table("WT_PR_AHMR_length.txt",sep="\t")
WT_PR_exp=all_exp[,c(2:5)]
#keep=rowSums(WT_PR_exp)>=2
#WT_PR_exp=WT_PR_exp[keep,]
row.names(WT_PR_exp)=all_exp$V1
colnames(WT_PR_exp)=c("WT1","WT2","PR1","PR2")
d=DGEList(counts=WT_PR_exp)
rpkm_exp=rpkm(d,gene.length = all_exp$V8)
d$counts=rpkm_exp
d=calcNormFactors(d)
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)
d$samples$group=c("WT","WT","PR","PR")
et=exactTest(d)
p_val=cbind(2^-et$table$logFC,et$table$PValue)
data_value=cbind(d$counts,d$pseudo.counts,p_val)
colnames(data_value)=c("WT1_ori","WT2_ori","PR1_ori","PR2_ori","WT1","WT2","PR1","PR2","fold","pValue")
write.table(data_value,"D6_WT_PR_rpkm.txt",sep="\t")

keep=rowSums(data_value[,5:8])>=4
data_value=data_value[keep,]
keep=data_value[,9] <= 0.5 | data_value[,9]>=2
data_value=data_value[keep,]
#keep=data_value[,9]>2 
#data_value=data_value[keep,]
b=as.matrix(data_value[,5:8])
#b=as.matrix(b)
library(gplots); 
myheatcol <- bluered(100)
hr <- hclust(as.dist(1-cor(t(b), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(b, method="spearman")), method="complete")
#mycl <- cutree(hr, h=max(hr$height)/1.2);
#mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=1); 
#mycolhc <- mycolhc[as.vector(mycl)]
heatmap.2(b, Rowv=as.dendrogram(hr), Colv=FALSE,col=myheatcol,keysize=1.6,
          scale="row", density.info="none", trace="none",margins=c(8,8))



# -----------------------------       WT Vs AM  -----------------------
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_exp=read.table("WT_PR_AHMR_length.txt",sep="\t")
WT_AM_exp=all_exp[,c(2:3,6:7)]
#keep=rowSums(WT_PR_exp)>=2
#WT_PR_exp=WT_PR_exp[keep,]
row.names(WT_AM_exp)=all_exp$V1
colnames(WT_AM_exp)=c("WT1","WT2","AM1","AM2")
d=DGEList(counts=WT_AM_exp)
rpkm_exp=rpkm(d,gene.length = all_exp$V8)
d$counts=rpkm_exp
d=calcNormFactors(d)
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)
d$samples$group=c("WT","WT","AM","AM")
et=exactTest(d)
p_val=cbind(2^-et$table$logFC,et$table$PValue)
data_value=cbind(d$counts,d$pseudo.counts,p_val)
colnames(data_value)=c("WT1_ori","WT2_ori","AM1_ori","AM2_ori","WT1","WT2","AM1","AM2","fold","pValue")
write.table(data_value,"D6_WT_AM_rpkm.txt",sep="\t")
keep=rowSums(data_value[,5:8])>=4
data_value=data_value[keep,]
keep=data_value[,9] <= 0.5 | data_value[,9]>=2
data_value=data_value[keep,]
#keep=data_value[,9]>2 
#data_value=data_value[keep,]
b=as.matrix(data_value[,5:8])
#b=as.matrix(b)
library(gplots); 
myheatcol <- bluered(100)
hr <- hclust(as.dist(1-cor(t(b), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(b, method="spearman")), method="complete")
#mycl <- cutree(hr, h=max(hr$height)/1.2);
#mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=1); 
#mycolhc <- mycolhc[as.vector(mycl)]

heatmap.2(b, Rowv=as.dendrogram(hr), Colv=FALSE,col=myheatcol,keysize=1.6,
          scale="row", density.info="none", trace="none",margins=c(8,8))


# -----------------------------       PR Vs AM  -----------------------
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")

all_exp=read.table("WT_PR_AHMR_length.txt",sep="\t")
PR_AM_exp=all_exp[,c(4:7)]

#keep=rowSums(WT_PR_exp)>=2
#WT_PR_exp=WT_PR_exp[keep,]

row.names(PR_AM_exp)=all_exp$V1
colnames(PR_AM_exp)=c("PR1","PR2","AM1","AM2")
d=DGEList(counts=PR_AM_exp)
rpkm_exp=rpkm(d,gene.length = all_exp$V8)
d$counts=rpkm_exp
d=calcNormFactors(d)
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)
d$samples$group=c("PR","PR","AM","AM")
et=exactTest(d)
p_val=cbind(2^-et$table$logFC,et$table$PValue)
data_value=cbind(d$counts,d$pseudo.counts,p_val)
colnames(data_value)=c("PR1_ori","PR2_ori","AM1_ori","AM2_ori","PR1","PR2","AM1","AM2","fold","pValue")
write.table(data_value,"D6_PR_AM_rpkm.txt",sep="\t")

keep=rowSums(data_value[,5:8])>=4
data_value=data_value[keep,]
keep=data_value[,9] <= 0.5 | data_value[,9]>=2
data_value=data_value[keep,]
#keep=data_value[,9]>2 
#data_value=data_value[keep,]
b=as.matrix(data_value[,5:8])
#b=as.matrix(b)
library(gplots); 
myheatcol <- bluered(100)
hr <- hclust(as.dist(1-cor(t(b), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(b, method="spearman")), method="complete")
#mycl <- cutree(hr, h=max(hr$height)/1.2);
#mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=1); 
#mycolhc <- mycolhc[as.vector(mycl)]

heatmap.2(b, Rowv=as.dendrogram(hr), Colv=FALSE,col=myheatcol,keysize=1.6,
          scale="row", density.info="none", trace="none",margins=c(8,8))

# -----------------------------       Xchoromesome genes  -----------------------
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
X_genes=read.table("X_genes.txt",sep="\t")
colnames(X_genes)=c("gene_name")
all_exp=read.table("WT_PR_AHMR_length.txt",sep="\t")
PR_AM_exp=all_exp[,c(2:7)]
#keep=rowSums(WT_PR_exp)>=2
#WT_PR_exp=WT_PR_exp[keep,]
row.names(PR_AM_exp)=all_exp$V1
colnames(PR_AM_exp)=c("WT1","WT2","PR1","PR2","AM1","AM2")
d=DGEList(counts=PR_AM_exp)
rpkm_exp=rpkm(d,gene.length = all_exp$V8)
d$counts=rpkm_exp
d=calcNormFactors(d)
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)
data_value=as.data.frame(d$pseudo.counts)
colnames(data_value)=c("WT1","WT2","PR1","PR2","AM1","AM2")
data_value$gene_name=row.names(data_value)

X_gene_exp=merge(x = data_value, y = X_genes, by = "gene_name", all.y  = TRUE)
write.table(X_gene_exp,"X_gene_exp.txt",sep="\t")
b=as.matrix(X_gene_exp[,2:5])
row.names(b)=X_gene_exp$gene_name
keep=rowSums(b)>=6
b=b[keep,]
#b=as.matrix(b)
library(gplots); 
myheatcol <- bluered(100)

hr <- hclust(as.dist(1-cor(t(b), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(b, method="spearman")), method="complete")
#mycl <- cutree(hr, h=max(hr$height)/1.2);
#mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=1); 
#mycolhc <- mycolhc[as.vector(mycl)]

heatmap.2(b, Rowv=as.dendrogram(hr), Colv=NULL,col=myheatcol,keysize=1.5,
          scale="row", density.info="none", trace="none",margins=c(8,8))





# -----------------------------       gene_anno_with_chromosome_info  11-19-2017 -----------------------

setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_gene=read.table("New_Ring1b_rpkm.txt",sep="\t",header = TRUE,stringsAsFactors=FALSE)
keep1=!grepl("Mir",row.names(all_gene))
all_exp_1=all_gene[keep1,]
small_exp=all_exp_1[,c(3:4,7:10)]
colnames(small_exp)=c("WT1","WT2","AM1","AM2","PR1","PR2")
small_exp$symbol=row.names(small_exp)
anno_file=read.table("F:/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2013-03-06-15-01-24/Genes/refFlat.txt",
                     sep="\t",header=FALSE,stringsAsFactors=FALSE)
colnames(anno_file)=c("symbol","ID","chr","strand","start","end","ID2","chr2","start2","end3","end4")
small_anno_file=subset(anno_file,select = c(symbol,ID,chr,start,end,strand))
annoed=merge(x=small_exp,y=small_anno_file,by="symbol",all.y =TRUE)

annoed1=annoed[!duplicated(annoed$ID),]

all_exp=annoed1[,c(1,9:12,2:8)]


#################                   remove and aggreate all "ID" in a same column    Descending 
c=all_exp[order(annoed1$symbol,decreasing = FALSE),]

library(dplyr)
gr=group_by(c,symbol)
Gr_index=group_indices(gr)


for (i in length(Gr_index):1) {
  # print("aa")
  if(i!=length(Gr_index) && Gr_index[i]==Gr_index[i+1]){
    #print("aa")
    c[i,13]=paste(c[i+1,12],",",c[i,12])
    #print(paste0("aa",i,"_",b[i-1,21],",",b[i,16])) 
  }else  {
    c[i,13]=c[i,12]
  }
}

##################                  Ascending 
#for (i in Gr_index) {
#  # print("aa")
#  if (i!=1 && Gr_index[i]==Gr_index[i-1]){
#    #print("aa")
#    b[i,21]=paste(b[i-1,21],",",b[i,16])
#    #print(paste0("aa",i,"_",b[i-1,21],",",b[i,16])) 
#  }else {
#    b[i,21]=b[i,16]
#  }
#}



d=c[!duplicated(c$symbol),] 
d[,12]=d[,13]
d=d[,-13]

###  remove miRNA and switch column
keep1=!grepl("Mir",d$symbol)
d=d[keep1,]


d=d[,c(1:5,12,6:11)]
colnames(d)=c("symbol","chr","start","end", "strand","ID","D6WT1",  "D6WT2" ,"D6AM1", "D6AM2","D6PR1","D6PR2")
write.table(d,"D6_WT_PR_AM_all_rpkm_annoed_with_chrom_v3.txt",sep="\t",row.names = FALSE)





################   calculate DEG between WT and PR   11-20-2017 ######
library(edgeR)
DL=DGEList(counts=d[,c(7:8,11:12)])
DL=calcNormFactors(DL)
DL=estimateCommonDisp(DL)
DL=estimateTagwiseDisp(DL)
DL$samples$group=c("WT","WT","PR","PR")
et=exactTest(DL)
p_val=cbind(2^-et$table$logFC,et$table$PValue)
data_value_PR=cbind(d[,1:6], DL$counts,DL$pseudo.counts,p_val)
colnames(data_value_PR)=c("Symbol","chr","start","end","strand","transcript_ID","WT1_ori","WT2_ori","PR1_ori","PR2_ori","WT1","WT2","PR1","PR2","fold","pValue")
write.table(data_value_PR,"data_value_PR_with_Fold.txt",sep="\t",row.names = FALSE)






################   calculate DEG between WT and AM   11-20-2017 ######
library(edgeR)
DL2=DGEList(counts=d[,c(7:10)])
DL2=calcNormFactors(DL2)
DL2=estimateCommonDisp(DL2)
DL2=estimateTagwiseDisp(DL2)
DL$samples$group=c("2","2","1","1")
et=exactTest(DL2)
p_val=cbind(2^et$table$logFC,et$table$PValue)
data_value_AM=cbind(d[,1:6], DL2$counts,DL2$pseudo.counts,p_val)
colnames(data_value_AM)=c("Symbol","chr","start","end","strand","transcript_ID","WT1_ori","WT2_ori","AM1_ori","AM2_ori","WT1","WT2","AM1","AM2","fold","pValue")
write.table(data_value_AM,"data_value_AM_with_Fold.txt",sep="\t",row.names = FALSE)

















# -----------------------------       gene_anno_with_chromosome_info  12-8-2016 -----------------------
####

setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_gene=read.table("D6_WT_PR_AM_all_rpkm.txt",sep="\t",header = TRUE,stringsAsFactors=FALSE)
anno_file=read.table("F:/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2013-03-06-15-01-24/Genes/refFlat.txt",
                     sep="\t",header=FALSE,stringsAsFactors=FALSE)
all_gene$symbol=row.names(all_gene)
colnames(anno_file)=c("symbol","ID","chr","strand","start","end","ID2","chr2","start2","end3","end4")
small_anno_file=subset(anno_file,select = c(symbol,ID,chr,start,end,strand))
annoed=merge(x=all_gene,y=small_anno_file,by="symbol",all.y =TRUE)
## remove p value
annoed=annoed[,-c(14:15)]   
## remove duplicated b$ID
annoed1=annoed[!duplicated(annoed$ID),]

###  remove miRNA and switch column
keep1=!grepl("Mir",annoed1$symbol)
all_exp=annoed1[keep1,]
all_exp=all_exp[,c(1,15:18,2:14)]


#################                   remove and aggreate all "ID" in a same column    Descending 
c=all_exp[order(all_exp$symbol,decreasing = FALSE),]

library(dplyr)
gr=group_by(c,symbol)
Gr_index=group_indices(gr)


for (i in length(Gr_index):1) {
  # print("aa")
  if(i!=length(Gr_index) && Gr_index[i]==Gr_index[i+1]){
    #print("aa")
    c[i,19]=paste(c[i+1,19],",",c[i,18])
    #print(paste0("aa",i,"_",b[i-1,21],",",b[i,16])) 
  }else  {
    c[i,19]=c[i,18]
  }
}

##################                  Ascending 
#for (i in Gr_index) {
#  # print("aa")
#  if (i!=1 && Gr_index[i]==Gr_index[i-1]){
#    #print("aa")
#    b[i,21]=paste(b[i-1,21],",",b[i,16])
#    #print(paste0("aa",i,"_",b[i-1,21],",",b[i,16])) 
#  }else {
#    b[i,21]=b[i,16]
#  }
#}
####

d=c[!duplicated(c$symbol),] 
d[,18]=d[,19]
d=d[,-19]
d=d[,-c(6:11)]
d=d[,c(1:5,12,6:11)]
colnames(d)=c("symbol","chr","start","end", "strand","ID","D6WT1",  "D6WT2" ,"D6PR1", "D6PR2"," D6AM1","D6AM2")
write.table(d,"D6_WT_PR_AM_all_rpkm_annoed_with_chrom_v2.txt",sep="\t",row.names = FALSE)


########          WT vs PR
####
library(edgeR)
D=DGEList(counts=d[,7:8])
D=calcNormFactors(D)
D=estimateCommonDisp(D)
D=estimateTagwiseDisp(D)
d1=D
d1$samples$group=c("1","1","2","2")
et=exactTest(d1)
WT_VS_PR=cbind(d[,1:10],2^et$table$logFC,et$table$PValue)
colnames(WT_VS_PR)=c("symbol","chr","start","end", "strand","ID","D6WT1",  "D6WT2" ,"D6PR1", "D6PR2","fold","pValue")
write.table(WT_VS_PR,"WT_VS_PR_with_Fold.txt",sep="\t",row.names = FALSE)

########          WT vs AM
####
library(edgeR)
D=DGEList(counts=d[,c(7:8,11:12)])
D=calcNormFactors(D)
D=estimateCommonDisp(D)
D=estimateTagwiseDisp(D)
d1=D
d1$samples$group=c("1","1","2","2")
et=exactTest(d1)
WT_VS_AM=cbind(d[,c(1:8,11:12)],2^et$table$logFC,et$table$PValue)
colnames(WT_VS_AM)=c("symbol","chr","start","end", "strand","ID","D6WT1",  "D6WT2" ,"D6AM1", "D6AM2","fold","pValue")
write.table(WT_VS_AM,"WT_VS_AM_with_Fold.txt",sep="\t",row.names = FALSE)

###############                   12-8-2016   All_down & all_up in PR and Amhr

#d=c[!duplicated(c$symbol),] 
#d$ID=d$V18
#d=d[,-18]
WT_PR_AM=cbind(WT_VS_PR,WT_VS_AM[,9:12])
colnames(WT_PR_AM)=c("symbol","chr","start","end", "strand","ID","D6WT1",  "D6WT2" ,"D6PR1", "D6PR2","PR_fold","PR_pValue","D6AM1", "D6AM2","AMHR_fold","AMHR_pValue")

keep=rowSums(WT_PR_AM[,c(7:10,13:14)])>5
WT_PR_AM_keep=WT_PR_AM[keep,]

#keep=rowSums(WT_VS_AM[,7:10])>5
#AM_keep=WT_VS_AM[keep,]

#write.table(d_keep,"D6_WT_PR_AM_all_rpkm_annoed_with_chrom_above_5.txt",sep="\t",row.names = FALSE)
#d_keep$PR=(rowMeans(d_keep[,11:12])+0.01)/(rowMeans(d_keep[,7:8])+0.01)
#d_keep$Amhr=(rowMeans(d_keep[,9:10])+0.01)/(rowMeans(d_keep[,7:8])+0.01)

keep_all=WT_PR_AM_keep$PR_fold>=2 & WT_PR_AM_keep$AMHR_fold>=2 & WT_PR_AM_keep$PR_pValue<0.05 &  WT_PR_AM_keep$AMHR_pValue<0.05
All_up=WT_PR_AM_keep[keep_all,]

keep_down=WT_PR_AM_keep$PR_fold<=0.5 & WT_PR_AM_keep$AMHR_fold<0.5 & WT_PR_AM_keep$PR_pValue<0.05 &  WT_PR_AM_keep$AMHR_pValue<0.05
All_down=WT_PR_AM_keep[keep_down,]


write.table(All_up[,c(1:10,13:14)],"All_up_anno_v3.txt",sep="\t",row.names = FALSE)
write.table(All_down[,c(1:10,13:14)],"All_down_anno_v3.txt",sep="\t",row.names = FALSE)



















########                       Scatter_plot  &&   MA plot       ############
########                       12-4-2016  
########                                                        =========>>>>
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_gene=read.table("D6_WT_PR_AM_all_rpkm_annoed_with_chrom.txt",sep="\t",header = TRUE,stringsAsFactors=FALSE)
data_for_sc=all_gene[,c(1,12:17)]
head(data_for_sc)
keep=rowSums(data_for_sc[,2:7])>2
data_for_sc=data_for_sc[keep,]

mcor = cor(data_for_sc[,2:7],method = c("pearson"))

library(gplots)

col=colorRampPalette(c("blue", "yellow", "red"))(500)
col <- bluered(100)
heatmap(x = mcor, col = col, symm = TRUE)
heatmap.2(mcor,col=col,density.info="none", trace="none",
          margins=c(10,10),symm = TRUE)

####  sample_distance
dists = dist( t(data_for_sc[,2:7] ) )
mat = as.matrix( dists )
library("RColorBrewer")
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace="none",col=colorRampPalette(c("blue", "yellow", "red"))(500), margin=c(13, 13), density.info="none",key=TRUE)

heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13), density.info="none")
heatmap.2(mat, trace="none", col =bluered(75), margin=c(13, 13), density.info="none",key=TRUE)
heatmap.2(mat, trace="none", col =redgreen(75), margin=c(13, 13), density.info="none",key=FALSE)



library(ggplot2)
p <- ggplot(data_for_sc,aes(log(data_for_sc[,2],2),log(data_for_sc[,3],2)))
pic=p+geom_point(aes(color="a"),alpha=0.8)+coord_cartesian(xlim=c(0,15),ylim=c(0,15))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)
pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid





library(reshape2)

Exp=melt(data_for_sc,id.vars = c( "symbol","D6WT1.sam.1"))
p <- ggplot(Exp, aes(log(D6WT1.sam.1,2), log(value,2))) + 
  geom_point(aes(color="a"),alpha=0.8)+coord_cartesian(xlim=c(0,15),ylim=c(0,15))
p + facet_grid(. ~ variable)
ggsave(p + facet_grid(. ~ variable), file="WT1.pdf")


Exp2=melt(data_for_sc,id.vars = c( "symbol","D6WT2.sam.1"))
p <- ggplot(Exp2, aes(log(D6WT2.sam.1,2), log(value,2))) + 
  geom_point(aes(color="a"),alpha=0.8)+coord_cartesian(xlim=c(0,15),ylim=c(0,15))
ggsave(p + facet_grid(. ~ variable), file="WT2.pdf")

Exp3=melt(data_for_sc,id.vars = c( "symbol","D6PR1.sam.1"))
p <- ggplot(Exp3, aes(log(D6PR1.sam.1,2), log(value,2))) + 
  geom_point(aes(color="a"),alpha=0.8)+coord_cartesian(xlim=c(0,15),ylim=c(0,15))
ggsave(p + facet_grid(. ~ variable), file="PR1.pdf")


Exp4=melt(data_for_sc,id.vars = c( "symbol","D6PR2.sam.1"))
p <- ggplot(Exp4, aes(log(D6PR2.sam.1,2), log(value,2))) + 
  geom_point(aes(color="a"),alpha=0.8)+coord_cartesian(xlim=c(0,15),ylim=c(0,15))
ggsave(p + facet_grid(. ~ variable), file="PR2.pdf")


Exp5=melt(data_for_sc,id.vars = c( "symbol","D6AM1.sam.1"))
p <- ggplot(Exp5, aes(log(D6AM1.sam.1,2), log(value,2))) + 
  geom_point(aes(color="a"),alpha=0.8)+coord_cartesian(xlim=c(0,15),ylim=c(0,15))
ggsave(p + facet_grid(. ~ variable), file="AM1.pdf")

Exp6=melt(data_for_sc,id.vars = c( "symbol","D6AM2.sam.1"))
p <- ggplot(Exp6, aes(log(D6AM2.sam.1,2), log(value,2))) + 
  geom_point(aes(color="a"),alpha=0.8)+coord_cartesian(xlim=c(0,15),ylim=c(0,15))
ggsave(p + facet_grid(. ~ variable), file="AM2.pdf")




all_gene=read.table("WT_PR_AHMR.txt",sep="\t",header = FALSE,stringsAsFactors=FALSE)
data_for_sc=all_gene[,c(7:12)]


library(DESeq2)
a=as.matrix(all_gene[,2:7])
rld <- rlog(a, blind=FALSE)
plot(rld[,1:2], pch=16, cex=0.3)
plot(rld[,2:3], pch=16, cex=0.3)
rld=as.data.frame(rld)


library(ggplot2)
p <- ggplot(rld,aes(log(rld[,2],2),log(rld[,3],2)))
pic=p+geom_point()
  coord_cartesian(xlim=c(0,15),ylim=c(0,15))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)
pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid








##    Scatter_plot
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_gene=read.table("D6_WT_PR_AM_all_rpkm_annoed_with_chrom.txt",sep="\t",header = TRUE,stringsAsFactors=FALSE)
data_for_sc=all_gene[,c(1,12:17)]
head(data_for_sc)
keep=rowSums(data_for_sc[,2:7])>2
data_for_sc=data_for_sc[keep,]
library(ggplot2)

#con=log2((data_value[,5]+data_value[,6])/2)
#Deci=log2((data_value[,7]+data_value[,8])/2)

###  add mean of each genetype
data_for_sc=transform(data_for_sc,
                     con_mean=(data_for_sc[,2]+data_for_sc[,3])/2+1,
                     PR_mean=(data_for_sc[,5]+data_for_sc[,4])/2+1,
                     Amhr_mean=(data_for_sc[,7]+data_for_sc[,6])/2+1
                    )

data_for_sc$PR_fold=(data_for_sc[,9]+0.1)/(data_for_sc[,8]+0.1)
data_for_sc$Amhr_fold=(data_for_sc[,10]+0.1)/(data_for_sc[,8]+0.1)
data_for_sc$PR_flag=cut(data_for_sc$PR_fold,breaks=c(-Inf,0.5,2,Inf),labels=c("-1","0","1"))
data_for_sc$Amhr_flag=cut(data_for_sc$Amhr_fold,breaks=c(-Inf,0.5,2,Inf),labels=c("-1","0","1"))


library(ggplot2)
library(ggrepel)

p <- ggplot(data_for_sc, aes(log(con_mean,2), log(PR_mean,2)))
pic=p+geom_point(aes(color=factor(PR_flag)))+
  coord_cartesian(xlim=c(0,12),ylim=c(0,12))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)+
  geom_text_repel(
    data = data_for_sc[data_for_sc$symbol %in% c("Nr5a1","Hsd3b1","Cyp11a1","Akr1c18","Aldh1a1","Star","Hsd17b7"),],
    aes(label = symbol),
    size = 5,
    box.padding = unit(1, "lines"), 
    point.padding = unit(0.1, "lines") )


pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid


p <- ggplot(data_for_sc, aes(log(con_mean,2), log(Amhr_mean,2)))
pic=p+geom_point(aes(color=factor(Amhr_flag)))+
  coord_cartesian(xlim=c(0,12),ylim=c(0,12))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)+
  geom_text_repel(
    data = data_for_sc[data_for_sc$symbol %in% c("Nr5a1","Hsd3b1","Cyp11a1","Akr1c18","Aldh1a1","Star","Hsd17b7"),],
    aes(label = symbol),
    size = 5,
    box.padding = unit(1, "lines"), 
    point.padding = unit(0.1, "lines") )

pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid









########                       Venn                             ############
########                       12-5-2016  
########                                                        =========>>>>
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_gene=read.table("D6_WT_PR_AM_all_rpkm_annoed_with_chrom.txt",sep="\t",header = TRUE,stringsAsFactors=FALSE)
install.packages('venneuler')
install.packages("ggplot2")
library(venneuler)
library(ggplot2)
keep=!grepl("random",all_gene$chr)
gene_rem=all_gene[keep,]
gene_rem$PR_fold=(gene_rem[,14]+gene_rem[,15]+0.1)/(gene_rem[,12]+gene_rem[,13]+0.1)
gene_rem$AMHR_fold=(gene_rem[,16]+gene_rem[,17]+0.1)/(gene_rem[,12]+gene_rem[,13]+0.1)

keep_pr=(gene_rem[,14]+gene_rem[,15]+gene_rem[,12]+gene_rem[,13])>5
gene_rem_pr=gene_rem[keep_pr,]

keep_am=(gene_rem[,16]+gene_rem[,17]+gene_rem[,12]+gene_rem[,13])>5
gene_rem_am=gene_rem[keep_am,]




p=ggplot(gene_rem_pr, aes(chr, log(PR_fold,2))) + geom_jitter(aes(color=chr))+ylim(0,7)+geom_text(aes(symbol=="Nr5a1"))
p
b=ggplot(gene_rem_am, aes(chr, log(AMHR_fold,2))) + geom_jitter(aes(color=chr))+ylim(0,7)
b

p





C2$gene_name=paste(C2$chrom,":",C2$start,"-",C2$end,sep="")
F1_EiPSC$gene_name=paste(F1_EiPSC$chrom,":",F1_EiPSC$start,"-",F1_EiPSC$end,sep="")
H1$gene_name=paste(H1$chrom,":",H1$start,"-",H1$end,sep="")
iPSC_CM$gene_name=paste(iPSC_CM$chrom,":",iPSC_CM$start,"-",iPSC_CM$end,sep="")

#install.packages( "gplots" )

venn_in=  list(C2=C2$gene_name, F1_EiPSC=F1_EiPSC$gene_name,  iPSC_CM=iPSC_CM$gene_name)

library(gplots)

venn( 
  venn_in
)


library(VennDiagram)
venn.diagram(venn_in,"VennDiagram.tiff",col = "transparent",
             fill = c("skyblue", "pink1","mediumorchid"  ),
             alpha = 0.50
             #label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
             # cat.col = c("skyblue", "pink1","mediumorchid", "orange")
             #c("darkblue", "darkgreen", "orange", "darkorchid4")
)













########       Pathway_analysis_KEGG      1-27-2018   =========>>>>
                                                     
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
source("https://bioconductor.org/biocLite.R")
biocLite("gage")
biocLite("pathview")


library(gage)
kg.mouse<- kegg.gsets("mouse")
kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
lapply(kegg.gs[1:3],head)
library(pathview)

all_gene=read.table("diff_exp_2018.txt",sep="\t",header = TRUE,stringsAsFactors=FALSE)
gene.symbol.eg= id2eg(ids=all_gene$Symbol, category='SYMBOL', org='Mm') 
all_gene$eg = gene.symbol.eg[,2]
all_gene=all_gene[!duplicated(all_gene$eg),]
keep1=is.na(all_gene$eg)!=TRUE
all_gene=all_gene[keep1,]

data_for_gage=''
data_for_gage=as.data.frame(all_gene$fold)


row.names(data_for_gage)=all_gene$eg

#data_for_gage=cbind(all_gene$fold,all_gene$eg)
fc.kegg.p<- gage(data_for_gage, gsets= kegg.gs, ref=NULL, samp=NULL)
write.table(fc.kegg.p,"b.txt",sep="\t")
rownames(fc.kegg.p)


###  up-regulate genes
sel<- fc.kegg.p$greater[,"p.val"] < 0.8 & !is.na(fc.kegg.p$greater[,"p.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
path.ids.great <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids.great, function(pid) pathview(gene.data = data_for_gage, pathway.id = pid,species = "mmu"))
###  down-regulate genes
sel.l<- fc.kegg.p$less[,"p.val"] < 0.8 & !is.na(fc.kegg.p$less[,"p.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids.l2 <- substr(path.ids.l, 1, 8)
path.ids.l2=path.ids.l2[-1]
pv.out.list <- sapply(path.ids.l2 [1:8], function(pid) pathview(gene.data = data_for_gage, pathway.id = pid,species = "mmu"))
###   Both
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
#path.ids3=path.ids2[-1]
pv.out.list <- sapply(path.ids2[1:20], function(pid) pathview(gene.data = data_for_gage, pathway.id = pid,species = "mmu"))



















setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
all_gene=read.table("D6_WT_PR_AM_all_rpkm_annoed_with_chrom.txt",sep="\t",header = TRUE,stringsAsFactors=FALSE)
all_gene$PR_fold=(all_gene[,14]+all_gene[,15]+0.1)/ (all_gene[,12]+all_gene[,13]+0.1)
keep=(all_gene[,14]+all_gene[,15]+0.1)/ (all_gene[,12]+all_gene[,13]+0.1)<= 0.5 | (all_gene[,14]+all_gene[,15]+0.1)/ (all_gene[,12]+all_gene[,13]+0.1)>=2
PR_diff=all_gene[keep,]
a=matrix(diff_gene$log10Fold)
rownames(a)=diff_gene$X
source("https://bioconductor.org/biocLite.R")
biocLite("gage")
biocLite("pathview")
library(gage)
kg.mouse<- kegg.gsets("mouse")
kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
lapply(kegg.gs[1:3],head)
library(pathview)
gene.symbol.eg<- id2eg(ids=PR_diff$symbol, category='SYMBOL', org='Mm') 

PR_diff$eg<- gene.symbol.eg[,2]
keep1=is.na(PR_diff$eg)!=TRUE
PR_diff=PR_diff[keep1,]

data_for_gage=as.data.frame(PR_diff$PR_fold)
row.names(data_for_gage)=PR_diff$eg

fc.kegg.p<- gage(data_for_gage, gsets= kegg.gs, ref=NULL, samp=NULL)
write.table(fc.kegg.p,"b.txt",sep="\t")
rownames(fc.kegg.p)
###  up-regulate genes
sel<- fc.kegg.p$greater[,"p.val"] < 0.6 & !is.na(fc.kegg.p$greater[,"p.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
path.ids.great <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids.great, function(pid) pathview(gene.data = data_for_gage, pathway.id = pid,species = "mmu"))
###  down-regulate genes
sel.l<- fc.kegg.p$less[,"p.val"] < 0.6 & !is.na(fc.kegg.p$less[,"p.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids.l2 <- substr(path.ids.l, 1, 8)
path.ids.l2=path.ids.l2[-1]
pv.out.list <- sapply(path.ids.l2 [1:8], function(pid) pathview(gene.data = data_for_gage, pathway.id = pid,species = "mmu"))
###   Both
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
#path.ids3=path.ids2[-1]
pv.out.list <- sapply(path.ids2[1:10], function(pid) pathview(gene.data = data_for_gage, pathway.id = pid,species = "mmu"))



















#####   AMCre/WT
files=c("D6WT1.txt","D6WT2.txt","D6AM1.txt","D6AM2.txt")
counts=readDGE(files) $counts
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__alignment_not_unique","__too_low_aQual","__not_aligned")
cpms = cpm(counts)
keep=rowSums(cpms)>=3& !noint
counts=counts[keep,]
colnames(counts) = c("WT1","WT2","AMCre1","AMCre2")
d = DGEList(counts=counts)
d = calcNormFactors(d)
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
d$samples$group = rep(c("WT", "AM"), c(2, 2))
et <- exactTest(d, dispersion = "auto")
data_value=cbind(as.data.frame(d$counts),as.data.frame(d$pseudo.counts),as.data.frame(et$table))
write.table(data_value,"D6_WT_AMCre_all.txt",sep="\t")






### PRCre vs AMCre
files=c("D6PR1.txt","D6PR2.txt","D6AM1.txt","D6AM2.txt")
counts=readDGE(files) $counts
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__alignment_not_unique","__too_low_aQual","__not_aligned")
cpms = cpm(counts)
keep=rowSums(cpms)>=3& !noint
counts=counts[keep,]
colnames(counts) = c("PRCre1","PRCre2","AMCre1","AMCre2")
d = DGEList(counts=counts)
d = calcNormFactors(d)
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
d$samples$group = rep(c("PRCre", "AMCre"), c(2, 2))
et <- exactTest(d, dispersion = "auto")
data_value=cbind(as.data.frame(d$counts),as.data.frame(d$pseudo.counts),as.data.frame(et$table))
write.table(data_value,"D6_PR_AMCre_all.txt",sep="\t")






















####      Heatmap
b=as.matrix(data_value[,c(9,10,11,12)])
#b=as.matrix(b)
library(gplots); 
myheatcol <- bluered(100)

hr <- hclust(as.dist(1-cor(t(b), method="pearson")), method="complete"); 
hc <- hclust(as.dist(1-cor(b, method="spearman")), method="complete")
#mycl <- cutree(hr, h=max(hr$height)/1.2);
#mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=1); 
#mycolhc <- mycolhc[as.vector(mycl)]

heatmap.2(b, Rowv=as.dendrogram(hr), col=myheatcol,
          scale="row", density.info="none", trace="none",margins=c(8,8))







#a=read.table("all_heatmap.txt",sep="\t",header=T, row.names = 1)

b=as.matrix(data_value[,c(7,8,9,10,11,12)])
#b=as.matrix(as.data.frame(a))
#install.packages("Hmisc")
#c=t(b)
#d=c[,1:1000]
#library(Hmisc)
#c=a[1:10,]
#b=t(a)
b=a
mcor = cor(b,method = c("pearson"))
#mcor=rcorr(a)
library(gplots)
#col = bluered(100)
col=colorRampPalette(c("red", "white", "blue"))(200)

#col <- bluered(100)
#heatmap(x = mcor, col = col, symm = TRUE)

heatmap.2(mcor,col=col,density.info="none", trace="none",
          margins=c(8,8),symm = TRUE)


































###  sample to sample

dists = dist( t(b ) )
mat = as.matrix( dists )
library("RColorBrewer")
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol <- bluered(100)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13), density.info="none")
heatmap.2(mat, trace="none", col =bluered(75), margin=c(13, 13), density.info="none")








###    Scatter_plot
library(ggplot2)
#con=log2((data_value[,5]+data_value[,6])/2)
#Deci=log2((data_value[,7]+data_value[,8])/2) 
data_value=transform(data_value,
                     con_mean=(data_value[,5]+data_value[,6])/2,
                     Deci_mean=(data_value[,7]+data_value[,8])/2
)

data_value$flag=cut(data_value$logFC,breaks=c(-Inf,-1,1,Inf),labels=c("-1","0","1"))
#data_value$flag[data_value$logFC>=1]="1"
#data_value$flag[data_value$logFC<=-1]="-1"
#data_value$flag[data_value$logFC<1& data_value$logFC>-1]="0"

p <- ggplot(data_value, aes(log(con_mean,2), log(Deci_mean,2)))
pic=p+geom_point(aes(color=factor(flag)))+
  coord_cartesian(xlim=c(4,20),ylim=c(4,20))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)

pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid


p <- ggplot(data_value, aes(log(Deci1,2), log(Deci2,2)))

pic=p+geom_point(aes(color=factor(flag)))+
  coord_cartesian(xlim=c(0,20),ylim=c(0,20))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)

pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid

#####    MA plot  

p <- ggplot(data_value,aes(log2(sqrt(con_mean*Deci_mean)),log2(Deci_mean/con_mean)))
pic=p+geom_point(aes(color=factor(flag)))+coord_cartesian(xlim=c(5,20),ylim=c(-10,10))+ scale_color_discrete(h=c(100,350), c=100, l=60)
pic
write.table(data_value,"ESC_con_EP22.txt",sep="\t")




#################             11-2-2016
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
All_exp=read.table("D6_WT_PR_AM_all_rpkm_annoed.txt",header=TRUE,sep="\t")


All_exp[,21]=(rowMeans(All_exp[,10:11])+0.01)/(rowMeans(All_exp[,8:9])+0.01)
All_exp[,22]=(rowMeans(All_exp[,12:13])+0.01)/(rowMeans(All_exp[,8:9])+0.01)

keep=All_exp[,21]>=2 & All_exp[,22]>=2
All_up=All_exp[keep,]

keep_down=All_exp[,21]<=0.5 & All_exp[,22]<=0.5
All_down=All_exp[keep_down,]
head(All_down)
write.table(All_up,"All_up_anno.txt",sep="\t")
write.table(All_down,"All_down_anno.txt",sep="\t")








######    Common target of ER, PR and Ring1b 1-28-2018   #####  
setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")





















########### ###   old code ---------------------------------------------------------------------



#  source("http://bioconductor.org/biocLite.R")
#  biocLite("edgeR")


library(edgeR)

files=c("D6WT1.txt","D6WT2.txt","D6PR1.txt","D6PR2.txt","D6AM1.txt","D6AM2.txt")
counts=readDGE(files) $counts

noint = rownames(counts) %in% c("__no_feature","__ambiguous","__alignment_not_unique","__too_low_aQual","__not_aligned")
cpms = cpm(counts)
#   keep=rowSums(cpms)>=0& !noint
keep=rowSums(cpms)>=3 & !noint
counts=counts[keep,]
colnames(counts) = c("WT1","WT2","PRCre1","PRCre2","AMCre1","AMCre2")
d = DGEList(counts=counts)
d = calcNormFactors(d)

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

#d$samples$group = rep(c("Control", "other"), c(2, 2))

#et <- exactTest(d, dispersion = "auto")

data_value=cbind(as.data.frame(d$counts),as.data.frame(d$pseudo.counts))

write.table(data_value,"D6_WT_PR_AMCre_all.txt",sep="\t")








######################  Ring1b and ER,PR    1-28-2018   ########

setwd("F:/New_Genome_data/Fenghua/Ring1b_chip_seq")
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
library(ChIPseeker)
library(GenomicFeatures)


Ring1b_peaks=read.table("F:/New_Genome_data/Fenghua/Ring1b_chip_seq/Ring1b_control_peaks.txt",header=TRUE, sep="\t",stringsAsFactors = F)
ER_peaks=read.table("E:/Genomic_data/E2_ChIP_seq_Uterus/Cal_data/Macs2/ER_E2/WT_1h_E2_ERalpha_peaks.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
PR_peaks=read.table("E:/Genomic_data/E2_ChIP_seq_Uterus/Cal_data/Macs2/PR_P4/P4_PR_peaks.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

txdb=makeTxDbFromGFF("F:/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2013-03-06-15-01-24/Genes/genes.gtf", 
                     format=c("auto"), 
                     organism="Mus musculus",circ_seqs=character())


Ring1b_peak_anno=annotatePeak("F:/New_Genome_data/Fenghua/Ring1b_chip_seq/Ring1b_control_peaks.txt",tssRegion = c(-3000,3000),TxDb=txdb,
                          genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
write.table(Ring1b_peak_anno,file="F:/New_Genome_data/Fenghua/Ring1b_chip_seq/Ring1b_control_peaks_annoed.txt",sep="\t")


ER_peak_anno=annotatePeak("E:/Genomic_data/E2_ChIP_seq_Uterus/Cal_data/Macs2/ER_E2/WT_1h_E2_ERalpha_peaks.txt",tssRegion = c(-3000,3000),TxDb=txdb,
                       genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
write.table(ER_peak_anno,file="E:/Genomic_data/E2_ChIP_seq_Uterus/Cal_data/Macs2/ER_E2/WT_1h_E2_ERalpha_peaks_annoed.txt",sep="\t")


PR_peak_anno=annotatePeak("E:/Genomic_data/E2_ChIP_seq_Uterus/Cal_data/Macs2/PR_P4/P4_PR_peaks.txt",tssRegion = c(-3000,3000),TxDb=txdb,
                       genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"))
write.table(PR_peak_anno,file="E:/Genomic_data/E2_ChIP_seq_Uterus/Cal_data/Macs2/PR_P4/P4_PR_peaks_annoed.txt",sep="\t")



### find overlap between ER PR and Ring1b


library(GenomicRanges)
#Rin1b_gr=GRanges(seqnames=data.frame(Ring1b_peak_anno@anno)$seqnames,IRanges(start=as.numeric(as.character(Ring1b_peak_anno$start)),end=as.numeric(as.character(Ring1b_peak_anno$end))))
ER_gr=GRanges(seqnames=ER_peaks$chr,IRanges(start=as.numeric(as.character(ER_peaks$start)),end=as.numeric(as.character(ER_peaks$end))))
PR_gr=GRanges(seqnames=PR_peaks$chr,IRanges(start=as.numeric(as.character(PR_peaks$start)),end=as.numeric(as.character(PR_peaks$end))))

ER_gr=GRanges(ER_peaks)
PR_gr=GRanges(PR_peaks)

#m=data.frame(Ring1b_peak_anno@anno)
#m$peak=paste(m$seqnames,m$start,m$end)
#n=m[order(m$peak,-m$pileup),]
#n=n[!duplicated(n$peak),]


#Ring_db=GRanges(seqnames=Ring1b_site$chr,IRanges(start=as.numeric(as.character(Ring1b_site$start)),end=as.numeric(as.character(Ring1b_site$end))))

###  ER Ring1b
a=findOverlaps(Ring1b_peak_anno@anno, ER_gr)
Ring1b_comm=Ring1b_peak_anno@anno[queryHits(a),]
ER_comm=ER_gr[subjectHits(a),]


b=data.frame(Ring1b_comm)
#b$index_1=1:length(b$seqnames)
b$peak=paste(b$seqnames,b$start,b$end)
o=order(b$peak,-b$pileup)
c=b[o,]
Ring1b_comm_s=c[!duplicated(c$peak),]

b_ER=data.frame(ER_comm)
#b$index_1=1:length(b$seqnames)
b_ER$peak=paste(b_ER$seqnames,b_ER$start,b_ER$end)
o_ER=order(b_ER$peak,-b_ER$pileup)
c_ER=b[o_ER,]

ER_comm_s=c_ER[!duplicated(c_ER$peak),]



###   ER PR
a=findOverlaps(ER_peak_anno@anno, PR_peak_anno@anno)
ER_comm=ER_peak_anno@anno[queryHits(a),]
#PR_comm=PR_gr[subjectHits(a),]


b=data.frame(ER_comm)
#b$index_1=1:length(b$seqnames)
b$peak=paste(b$seqnames,b$start,b$end)
o=order(b$peak,-b$pileup)
c=b[o,]
ER_comm_s=c[!duplicated(c$peak),]


write.table(ER_comm_s,"ER_PR_Final_site.txt",sep="\t")



###  ER Ring1b
a=findOverlaps(Ring1b_peak_anno@anno, ER_gr)
Ring1b_comm=Ring1b_peak_anno@anno[queryHits(a),]
ER_comm=ER_gr[subjectHits(a),]


b=data.frame(Ring1b_comm)
#b$index_1=1:length(b$seqnames)
b$peak=paste(b$seqnames,b$start,b$end)
o=order(b$peak,-b$pileup)
c=b[o,]
Ring1b_comm_s=c[!duplicated(c$peak),]

b_ER=data.frame(ER_comm)
#b$index_1=1:length(b$seqnames)
b_ER$peak=paste(b_ER$seqnames,b_ER$start,b_ER$end)
o_ER=order(b_ER$peak,-b_ER$pileup)
c_ER=b[o_ER,]

ER_comm_s=c_ER[!duplicated(c_ER$peak),]














write.table(ER_comm_s,"ER_Final_site.txt",sep="\t")
write.table(Ring1b_comm_s,"Ring1b_Final_site.txt",sep="\t")


write.table(Final_site_ring,"Final_site_ring.txt",sep="\t")




diff_gene=read.table("diff_exp_2018.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)
a=data.frame(diff_gene$Symbol)
colnames(a)=c("gene_name")

annoed_gene=data.frame(Ring1b_peak_anno)
b=data.frame(annoed_gene[!duplicated(annoed_gene$geneId),])
b$gene_name=b$geneId

target_ab=merge(b,a,by="gene_name",all.a=TRUE)
write.table(target_ab,"target_ab.txt",sep="\t")




#######   Diff TFs   ####

setwd("F:/New_Genome_data/Fenghua/RNA-Seq/D6_WT_PGR-Cre_Amhr-Cre/RNA-seq")
All_exp=read.table("D6_WT_PR_AM_all_rpkm_annoed_with_chrom_v3.txt",header=TRUE,sep="\t")
#All_exp$Gene_name=All_exp$symbol

All_TF=read.table("F:/New_Genome_data/All_TFs/all_tfs.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)
All_DNA=read.table("F:/New_Genome_data/All_TFs/DNAbinding_plus_TF_db.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)


aaa=All_exp[All_exp$symbol %in% All_DNA$Gene_name,]

write.table(aaa,"temp.txt",sep="\t")
