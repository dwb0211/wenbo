#########  ########             edgeR  RPKM
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
