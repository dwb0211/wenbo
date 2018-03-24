#########  ########             Mac code
setwd("/Volumes/EXT/D4_6_8")
source("https://bioconductor.org/biocLite.R")
#biocLite("Rsubread")
#biocLite("edgeR")
library(Rsubread)
library(edgeR)

Exp=featureCounts(files=c("D4.bam",
                          "D6.bam",
                          "D8.bam"),
                  annot.ext = "/Volumes/DATA/Mus_musculus_UCSC_mm9/Mus_musculus/UCSC/mm9/Annotation/Archives/archive-2013-03-06-15-01-24/Genes/genes.gtf",
                  isGTFAnnotationFile = TRUE,
                  nthreads = 8)

#write.table(cbind(Exp$counts,Exp$annotation$Length),"C2_CM_iPSC_mRNA_length.txt",sep="\t",quote = FALSE,col.names =FALSE)
rpkm_exp=rpkm(Exp$counts,gene.length = Exp$annotation$Length)
D=DGEList(counts=rpkm_exp)
d=calcNormFactors(D)
d=estimateCommonDisp(d)
d=estimateTagwiseDisp(d)

d$samples$group=c(1,2,1)
et=exactTest(d)
p_val=cbind(2^et$table$logFC,et$table$PValue)
colnames(p_val) = c("Fold","pValue")

data_value=cbind(d$counts,d$pseudo.counts,p_val)
colnames(data_value)=c("D4","D6","D8","D4_nor","D6_nor","D8_nor","Fold","p_Val" )
write.table(data_value,"D4_6_8_rpkm.txt",sep="\t")


##################################    windows code

setwd("F:/New_Genome_data/D4_6_8")
exp=read.table("D4_6_8_rpkm.txt",sep="\t",header=TRUE, stringsAsFactors = F)
colnames(exp)=c("gene_name","D4","D6","D8","b","D4_nor","D6_nor","D8_nor","a","Fold","p_Val" )

TF_db=read.table("F:/New_Genome_data/All_TFs/DNAbinding_plus_TF_db.txt",sep="\t", header=T,stringsAsFactors = FALSE)
colnames(TF_db)=c("gene_name","description")
TF_exp= merge(x =TF_db , y = exp, by = "gene_name", all.x=TRUE)
write.table(TF_exp,"D4_6_8_TFs_exp.txt",sep="\t", row.names = F,col.names=T,quote=F)



TF_db=read.table("F:/New_Genome_data/All_TFs/all_tfs.txt",sep="\t", header=T,stringsAsFactors = FALSE)
colnames(TF_db)=c("gene_name","description")
TF_exp= merge(x =TF_db , y = exp, by = "gene_name", all.x=TRUE)
write.table(TF_exp,"D4_6_8_TFs_only_exp.txt",sep="\t", row.names = F,col.names=T,quote=F)
