#########        edgeR for RPKM calculation

setwd("F:/New_Genome_data/MEF_and_uterus")



#  source("http://bioconductor.org/biocLite.R")
#  biocLite("edgeR")


library(edgeR)

files=c("D4_nor.txt","p53WT_con.txt")
counts=readDGE(files) $counts

noint = rownames(counts) %in% c("__no_feature","__ambiguous","__alignment_not_unique","__too_low_aQual","__not_aligned")
cpms = cpm(counts)
keep=rowSums(cpms)>=0& !noint
keep=rowSums(cpms)>=3 & !noint
counts=counts[keep,]
colnames(counts) = c("D4","MEF")
d = DGEList(counts=counts)
d = calcNormFactors(d)

d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

d$samples$group = c("D4", "MEF")

et <- exactTest(d, dispersion = "auto")

data_value=cbind(as.data.frame(d$counts),as.data.frame(d$pseudo.counts),as.data.frame(et$table))

write.table(data_value,"D4_MEF.txt",sep="\t")



TF_db=read.table("All_transciption_factors.txt",sep="\t", quote = "")
exp_mat=read.table("D4_MEF.txt",sep="\t")
exp_mat$gene_id=row.names(exp_mat)

colnames(TF_db)=c("gene_id","gene_name")

TF_exp= merge(x =TF_db , y = exp_mat, by = "gene_id", all.x=TRUE)
write.table(TF_exp,"TF_exp.txt",sep="\t")




##############   Transription screen    1-16-2018  #####


setwd('F:/New_Genome_data/V2Ltf_LE_St')
v1v2_exp=read.table("LE_St_LtfVangl2.txt",sep="\t", header=TRUE,stringsAsFactors = FALSE)

Epi_st=v1v2_exp[,c(1:6)]
Epi_st$gene_name=row.names(Epi_st)


TF_db=read.table("F:/New_Genome_data/All_TFs/DNAbinding_plus_TF_db.txt",sep="\t", header=T,stringsAsFactors = FALSE)
colnames(TF_db)=c("gene_name","description")

TF_exp= merge(x =TF_db , y = Epi_st, by = "gene_name", all.x=TRUE)

write.table(TF_exp,"TF_exp_Jan_20th_2018.txt",sep="\t", row.names = F,col.names=T,quote=F)





###  for figure

exp=read.table("LE_Str_specific_TF_exp_for_R.txt",sep="\t", quote="", header=TRUE,stringsAsFactors = FALSE)

exp=exp[order(exp$Fold),]


exp$mean_LE=(rowSums(exp[,c(3:5)])+1)/3

exp$mean_ST=(rowSums(exp[,c(6:8)])+1)/3
row.names(exp)=exp$gene_name

exp_H=rbind(head(exp,20),tail(exp,20))


library(gplots); 
myheatcol <- bluered(10)

b=as.matrix(exp_H[,c(10,11)])
#b=b[order(b[,1],decreasing = TRUE),]
c=log2(b)


#c[is.infinite(c)]=0
c=c-max(c)/2
library("RColorBrewer")
hmcol=colorRampPalette(c("red", "yellow", "blue"))(200)
heatmap.2(c, trace="none", col = rev(hmcol), Colv=FALSE, margin=c(1, 13), density.info="none",keysize = 1.1)
heatmap.2(c, trace="none", col = rev(hmcol),Rowv=FALSE, Colv=FALSE, margin=c(5, 5), density.info="none",keysize = 1.1)




