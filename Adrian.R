#########  ########             edgeR  RPKM
setwd("/Volumes/EXT/Adrian")
source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")
biocLite("edgeR")
library(Rsubread)
library(edgeR)
################    10-23-2016

#temp<-list.files()
#for (i in 60:83){
#  list.myfiles=paste("/Volumes/EXT/Adrian/SRR61976",i,"accepted_hits.bam")
#}


files=dir(patter="^SRR")
myfiles=""
for (i in 1:24){
  myfiles[i]=paste(files[i],"/accepted_hits.bam",sep="")
  
}

Exp=featureCounts(myfiles, annot.ext = "/Volumes/DATA/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2013-03-06-15-01-24/Genes/genes.gtf",
                  isGTFAnnotationFile = TRUE,
                  nthreads = 8)
#write.table(Exp$counts,"CB1CB2KO.txt",sep="\t",quote = FALSE,col.names =FALSE)
#write.table(cbind(Exp$counts,Exp$annotation$Length),"CB1CB2KO_xc length.txt",sep="\t",quote = FALSE,col.names =FALSE)
D=DGEList(counts=Exp$counts)
D=calcNormFactors(D)
rpkm_exp=rpkm(D,gene.length = Exp$annotation$Length)
D$counts=rpkm_exp
colnames(D$counts)=c("E15-DSC-1","E7-DSC-1","E7-DSC-TI-1","E15-DSC-2","E7-DSC-2","E7-DSC-TI-2",
                     "E15-DSC-3","E7-DSC-3","E7-DSC-TI-3","E15-MSC-1","E7-MSC-1","E7-MSC-TI-1",
                     "E15-MSC-2","E7-MSC-2","E7-MSC-TI-2","E15-MSC-3","E7-MSC-3","E7-MSC-TI-3",
                     "USC-1","USC-TI-1","USC-2","USC-TI-2","USC-3","USC-TI-3"
)
#D=estimateCommonDisp(D)
#D=estimateTagwiseDisp(D)
#D$samples$group=c("WT1","WT2","Vangl1KO1","Vangl1KO2","Vangl2KO1","Vangl2KO2",
#                "R1R2KO1","R1R2KO2","Wnt5aGOF1","Wnt5aLOF1","Wnt5aLOF2")
#et=exactTest(D)
#p_val=cbind(et$table$logFC,et$table$PValue)
#colnames(p_val) = c("log_Fold","pValue")
#data_value=cbind(D$counts,D$pseudo.counts)

data_value=D$counts

keep1=!grepl("Mir",row.names(data_value))
all_exp_1=data_value[keep1,]

write.table(all_exp_1,"D8_D16_str_msc_Adrian.txt",sep="\t")
















#######  D8 immuno cells

#########  ########             edgeR  RPKM
setwd("/Volumes/EXT/Adrian")
source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")
biocLite("edgeR")
library(Rsubread)
library(edgeR)

Exp=featureCounts(files=c("SRR6197661/accepted_hits.bam",
                          "SRR6197664/accepted_hits.bam",
                          "SRR6197667/accepted_hits.bam",
                          "/Volumes/EXT/D4_6_8/D8.bam"
),

annot.ext = "/Volumes/DATA/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2013-03-06-15-01-24/Genes/genes.gtf",
isGTFAnnotationFile = TRUE,
nthreads = 8)
#write.table(Exp$counts,"CB1CB2KO.txt",sep="\t",quote = FALSE,col.names =FALSE)
#write.table(cbind(Exp$counts,Exp$annotation$Length),"CB1CB2KO_xc length.txt",sep="\t",quote = FALSE,col.names =FALSE)
D=DGEList(counts=Exp$counts)
D=calcNormFactors(D)
rpkm_exp=rpkm(D,gene.length = Exp$annotation$Length)
D$counts=rpkm_exp
colnames(D$counts)=c("E7-DSC-1","E7-DSC-2","E7-DSC-3","D8")


data_value=D$counts

keep1=!grepl("Mir",row.names(data_value))
all_exp_1=data_value[keep1,]

write.table(all_exp_1,"D8_immune_cells.txt",sep="\t")



###########   for Figures    #############

setwd("F:/New_Genome_data/Adrian")
Exp=read.table("D8_D16_str_msc_Adrian.txt",sep="\t", header = TRUE, stringsAsFactors = FALSE)
#target_gene=data.frame(paste("Tlr",c(1:9),sep=""))

Tlr_gene=subset(Exp,row.names(Exp) %in% c(paste("Tlr",c(1:9),sep="")))
Tlr_gene=Tlr_gene[,c(2,5,8,11,14,17,1,4,7,10,13,16)]
tlr_s1=Tlr_gene[,c(1,4,7,10)]
tlr_s1$symbol=row.names(tlr_s1)
colnames(tlr_s1)=c("D8_Dec","D8_Msc","D16_Dec","D16_Msc","symbol")

tlr_s2=Tlr_gene[,c(2,5,8,11)]
tlr_s2$symbol=row.names(tlr_s2)
colnames(tlr_s2)=c("D8_Dec","D8_Msc","D16_Dec","D16_Msc","symbol")

tlr_s3=Tlr_gene[,c(3,6,9,12)]
tlr_s3$symbol=row.names(tlr_s3)
colnames(tlr_s3)=c("D8_Dec","D8_Msc","D16_Dec","D16_Msc","symbol")

tlr_s11=rbind(tlr_s1,tlr_s2)
Tlr_s22=rbind(tlr_s11,tlr_s3)

D_8d=Tlr_s22[,c(1,5)]
D_8d$sample=c("D8_Deci")
colnames(D_8d)=c("Exp","symbol","sample")

D_8m=Tlr_s22[,c(2,5)]
D_8m$sample=c("D8_Msc")
colnames(D_8m)=c("Exp","symbol","sample")

D_16d=Tlr_s22[,c(3,5)]
D_16d$sample=c("D16_Deci")
colnames(D_16d)=c("Exp","symbol","sample")

D_16m=Tlr_s22[,c(4,5)]
D_16m$sample=c("D16_Msc")
colnames(D_16m)=c("Exp","symbol","sample")

TLR_E=rbind(D_8d,D_8m,D_16d,D_16m)

library(reshape2)
library(ggplot2)
#iPSC_small2exon_exp2_melt=melt(iPSC_small2exon_exp,id.vars = c("gene_name","gene_id", "readNumber", "flag"))

##  D8 deci
TLR_small=TLR_E[c(1:54),]
#TLR_small$Exp=log((TLR_small$Exp+1),10)
p <- ggplot(TLR_small[1:27,], aes(symbol, Exp,fill="1"))+stat_boxplot(geom ='errorbar')+ geom_boxplot() +xlab("") +ylab("TLRs expression in day 8 Deci (FPKM)")

##  D8 myo
p <- ggplot(TLR_small[28:54,], aes(symbol, Exp,fill="1"))+stat_boxplot(geom ='errorbar')+ geom_boxplot() +xlab("") +ylab("TLRs expression in day 8 Myomytrium (FPKM)")
p+theme(panel.grid.major=element_line(colour=NA)) +theme(panel.grid=element_blank())


TLR_small=TLR_E[c(55:108),]
#TLR_small$Exp=log((TLR_small$Exp+1),10)

## D16 Deci
p <- ggplot(TLR_small[1:27,], aes(symbol, Exp,fill="1"))+stat_boxplot(geom ='errorbar')+ geom_boxplot() +xlab("") +ylab("TLRs expression in day 16 Deci (FPKM)")
p+theme(panel.grid.major=element_line(colour=NA)) +theme(panel.grid=element_blank())
## D16 Myo

p <- ggplot(TLR_small[28:54,], aes(symbol, Exp,fill="1"))+stat_boxplot(geom ='errorbar') + geom_boxplot() +xlab("") +ylab("TLRs expression in day 16 Myomytrium (FPKM)")
p+theme(panel.grid.major=element_line(colour=NA)) +theme(panel.grid=element_blank())

#p+theme_set(theme_bw())   # no backgroud
#p+  theme(panel.grid=element_blank())









###############    TFs in Deci and MSC  #############
setwd("F:/New_Genome_data/Adrian")
Exp=read.table("D8_D16_str_msc_Adrian.txt",sep="\t", header = TRUE, stringsAsFactors = FALSE)

exp=exp[,c(2,5,8,11,14,17,1,4,7,10,13,16)]
exp$Gene_name=row.names(exp)
  
TFs=read.table("F:/New_Genome_data/All_TFs/all_tfs.txt",header=TRUE,sep="\t",stringsAsFactors = F)
DNA_pr=read.table("F:/New_Genome_data/All_TFs/DNAbinding_plus_TF_db.txt",header=TRUE,sep="\t",stringsAsFactors = F)


TF_merge=merge(TFs,exp,by="Gene_name",all.x=TRUE)
write.table(TF_merge,"TFs_deci_myo.txt",sep="\t",row.names=FALSE,quote=FALSE)



DNA_pr_merge=merge(DNA_pr,exp,by="Gene_name",all.x=TRUE)
write.table(DNA_pr_merge,"DNA_bingding_protein_deci_myo.txt",sep="\t",row.names=FALSE,quote=FALSE)




