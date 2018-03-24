








######################    win_code    ##############   
######################   Read ciRNA_data   3-17-2018    ################
setwd("F:/New_Genome_data/Hum_deci_ciRNA") 
filenames=list("con_circ.txt","deci_circ.txt")
datalist <- lapply(filenames, function(name) {
  read.table(name, sep = "\t",header=FALSE)
})
names(datalist)=c("con_ci","deci_ci")
colnames(datalist$con_ci)=c("chrom","start",	"end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	
                                  "itemRgb","exonCount","exonSizes","exonOffsets",	"readNumber",
                                  "circType",	"geneName","isoformName","flankIntron","IntronSite")
colnames(datalist$deci_ci)=c("chrom","start",	"end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	
                                  "itemRgb","exonCount","exonSizes","exonOffsets",	"readNumber",
                                  "circType",	"geneName","isoformName","flankIntron","IntronSite")
datalist$con_ci$rpb=datalist$con_ci$readNumber/85403972*10^9
datalist$deci_ci$rpb=datalist$deci_ci$readNumber/78786509*10^9

small_datalist2=lapply(datalist,function(x) {
  cbind(
    paste(x$chrom,x$start,x$end),data.frame(x$geneName),data.frame(x$isoformName),data.frame(x$rpb),data.frame(x$readNumber),data.frame(x$exonCount))
})
#names(small_datalist2)=filenames
for(i in 1:length(small_datalist2)){ 
  names(small_datalist2[[i]])=c('gene_id','gene_name','isoformName',paste(names(small_datalist2[i]),'rpb',sep="_"),
                                paste(names(small_datalist2[i]),'readNumber',sep="_"),paste(names(small_datalist2[i]),'exoncount',sep="_"))
  #small_datalist[[i]]=small_datalist[[i]][!duplicated(small_datalist[[i]]$gene_id),]
}
merged_data = Reduce(function(...) merge(..., all=T), small_datalist2)
merged_data[is.na(merged_data)]=0
write.table(merged_data,"hESC_ciRNA.txt",sep="\t")


############################  Exon_counts ##############################
################              3-17-2016  
###########################################==========>>>>>

setwd("F:/New_Genome_data/Hum_deci_ciRNA") 
all_ci=read.table("hESC_ciRNA.txt",head=TRUE,sep="\t")
library(reshape2)
library(ggplot2)
all_ci_melt=melt(all_ci,id.vars = c( "gene_id", "gene_name","isoformName","con_ci_rpb","con_ci_readNumber","deci_ci_rpb","deci_ci_readNumber"))
a=ggplot(all_ci_melt,aes(value,color=variable,fill=variable))+geom_bar()+xlim(0,15)+ylim(0,500)  +xlab("Number of exons") +ylab("The number of ciRNAs contain different exons")
a+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) +  theme(panel.grid=element_blank())      # no backgroud
ggsave("exon_number.png",dpi=300)





##############################  All_ciRNAs_length_distribution  6-10-2017  ############################
################                  3-17-2018
###########################################==========>>>>>
setwd("F:/New_Genome_data/Hum_deci_ciRNA") 
filenames=list("con_circ.txt","deci_circ.txt")
datalist <- lapply(filenames, function(name) {
  read.table(name, sep = "\t",header=FALSE)
})
names(datalist)=c("con_ci","deci_ci")
colnames(datalist$con_ci)=c("chrom","start",	"end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	
                            "itemRgb","exonCount","exonSizes","exonOffsets",	"readNumber",
                            "circType",	"geneName","isoformName","flankIntron","IntronSite")
colnames(datalist$deci_ci)=c("chrom","start",	"end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	
                             "itemRgb","exonCount","exonSizes","exonOffsets",	"readNumber",
                             "circType",	"geneName","isoformName","flankIntron","IntronSite")
datalist$con_ci$rpb=datalist$con_ci$readNumber/85403972*10^9
datalist$deci_ci$rpb=datalist$deci_ci$readNumber/78786509*10^9

small_datalist2=lapply(datalist,function(x) {
  cbind(
    paste(x$chrom,x$start,x$end),data.frame(x$geneName),data.frame(x$rpb),data.frame(x$readNumber),data.frame(x$exonCount),data.frame(x$exonSizes))
})

for(i in 1:length(small_datalist2)){ 
  names(small_datalist2[[i]])=c('gene_id','gene_name','rpb','readNumber','exonCount','exon_size')  
}
small_datalist2$con_ci$flag=c("con_ci")
small_datalist2$deci_ci$flag=c("deci_ci")
long_mtr=rbind(small_datalist2$con_ci,small_datalist2$deci_ci) 
long_mtr$exon_length=''
for(i in 1:nrow(long_mtr)){ 
  long_mtr$exon_length[i]=sum(as.numeric(unlist(strsplit(as.character(long_mtr[i,]$exon_size),"[,]"))))
}
long_mtr$exon_length=as.numeric(long_mtr$exon_length)
write.table(long_mtr,"ciRNA_for_exon_size_3_17_2018.txt",sep="\t")

library(ggplot2)
a=ggplot(long_mtr,aes(flag,exon_length,color=flag))+geom_boxplot()+geom_jitter(width=0.2)+
  stat_boxplot(geom ='errorbar', width = 0.4) + xlab("") +ylab("Length of ciRNAs ")+ylim(0,5000)
a+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) +  theme(panel.grid=element_blank())      # no backgroud
ggsave("Length_of_ciRNAs.png",dpi=300)



############################  heatmap for high expressed ciRNAs #####################################
################            (Top 100)         3-17-2018
###########################################       ==========>>>>>
setwd("F:/New_Genome_data/Hum_deci_ciRNA") 
all_ci=read.table("hESC_ciRNA.txt",head=TRUE,sep="\t")
library(reshape2)
library(ggplot2)
rownames(all_ci)=all_ci$gene_id
exp_mat=all_ci[,-c(1,2,3,5,6,8,9)]
library(gplots); 
myheatcol <- bluered(10)
b=as.matrix(exp_mat)
colnames(b)=c("con","deci")
b=b[order(b[,1],decreasing = TRUE),]
c=log2(b)
c[is.infinite(c)]=0
c=c-max(c)/2
library("RColorBrewer")
hmcol=colorRampPalette(c("red", "yellow", "blue"))(200)
heatmap.2(c, trace="none", col = rev(hmcol), Colv=FALSE, margin=c(5,13), density.info="none",keysize = 1.1)
#heatmap.2(c, trace="none", col = rev(hmcol),Rowv=FALSE, Colv=FALSE, density.info="none",keysize = 1.1)










########################  how many ciRNAs each gene contain ##################################
################          3-17-2018
###########################################==========>>>>>
setwd("F:/New_Genome_data/Hum_deci_ciRNA") 
filenames=list("con_circ.txt","deci_circ.txt")
datalist <- lapply(filenames, function(name) {
  read.table(name, sep = "\t",header=FALSE)
})
names(datalist)=c("con_ci","deci_ci")
colnames(datalist$con_ci)=c("chrom","start",	"end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	
                            "itemRgb","exonCount","exonSizes","exonOffsets",	"readNumber",
                            "circType",	"geneName","isoformName","flankIntron","IntronSite")
colnames(datalist$deci_ci)=c("chrom","start",	"end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	
                             "itemRgb","exonCount","exonSizes","exonOffsets",	"readNumber",
                             "circType",	"geneName","isoformName","flankIntron","IntronSite")
datalist$con_ci$rpb=datalist$con_ci$readNumber/85403972*10^9
datalist$deci_ci$rpb=datalist$deci_ci$readNumber/78786509*10^9

small_datalist2=lapply(datalist,function(x) {
  cbind(
    paste(x$chrom,x$start,x$end),data.frame(x$geneName),data.frame(x$isoformName),data.frame(x$rpb),data.frame(x$readNumber),data.frame(x$exonCount))
})

for(i in 1:length(small_datalist2)){ 
  names(small_datalist2[[i]])=c('gene_id','gene_name','isoformName',paste(names(small_datalist2[i]),'rpb',sep="_"),
                                paste(names(small_datalist2[i]),'readNumber',sep="_"),paste(names(small_datalist2[i]),'exoncount',sep="_"))
  #small_datalist[[i]]=small_datalist[[i]][!duplicated(small_datalist[[i]]$gene_id),]
}

con_ci_freq=data.frame(table(small_datalist2$con_ci[,c(1,2)]$gene_name))
con_ci_freq$flag=c("con_ci")
deci_ci_freq=data.frame(table(small_datalist2$deci_ci[,c(1,2)]$gene_name))
deci_ci_freq$flag=c("deci_ci")


long_mtr=rbind(con_ci_freq,deci_ci_freq) 
library(ggplot2)

a=ggplot(long_mtr,aes(Freq,color=flag,fill=flag))+geom_bar() +xlim(0,8)+xlab("Number of ciRNAs each gene contain") +ylab("Number of Genes containing different ciRNAs")
a+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA)) +  theme(panel.grid=element_blank())      # no backgroud
write.table(long_mtr,file="number_of_ciRNAs_each_gene_has.txt",sep="\t",row.names = FALSE,quote=FALSE)
ggsave("number_of_ciRNAs_each_gene_has.png")









#####################  compare ciRNAs with known ciRNAs  12-3-2016  ##############
################          
################            12-3-2016
###########################################==========>>>>>
setwd("F:/New_Genome_data/Hum_deci_ciRNA") 
#all_ci=read.table("hESC_ciRNA.txt",head=TRUE,sep="\t")
Esc_con=read.table("con_circ.txt",header=FALSE,sep="\t")
Esc_con_small=cbind(paste(Esc_con$V1,Esc_con$V2,Esc_con$V3),data.frame(Esc_con$V15))
colnames(Esc_con_small)=c("ciRNA_id","gene_symbol")



Esc_deci=read.table("deci_circ.txt",header=FALSE,sep="\t")
Esc_deci_small=cbind(paste(Esc_deci$V1,Esc_deci$V2,Esc_deci$V3),data.frame(Esc_deci$V15))
colnames(Esc_deci_small)=c("ciRNA_id","gene_symbol")



###  old verison   better
known_ciRNA=read.table("hsa_hg19_circRNA.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)

small_known_ciRNA=data.frame(
  cbind(paste(known_ciRNA$chrom,known_ciRNA$start,known_ciRNA$end,sep=" "),known_ciRNA$strand,known_ciRNA$circRNA.ID,known_ciRNA$gene.symbol)
)

names(small_known_ciRNA)=c("gene_id","strand","ciRNA_id","gene_symbol")
venn_in=  list( Esc_con_ID=Esc_con_small$ciRNA_id,  Esc_deci_ID=Esc_deci_small$ciRNA_id,
                 known_ciRNAs=small_known_ciRNA$gene_id)
library(VennDiagram)
venn.diagram(venn_in,"VennDiagram.tiff",col = "transparent",
             fill = c("skyblue", "pink1","darkorchid4" ),
             alpha = 0.50
             #label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
             # cat.col = c("skyblue", "pink1","mediumorchid", "orange")
             #c("darkblue", "darkgreen", "orange", "darkorchid4")
)

##### all_ciRNA  new verison   too much, not very good for figure
KnownciRNA=read.table("KnownciRNA.txt",header=FALSE, sep="|")
KnownciRNA$V2=sub("(chr[0-9]*)(:)([0-9]*)(-)([0-9]*)([-|+])", "\\1 \\3 \\5", KnownciRNA$V2)
names(KnownciRNA)=c("HSA_ID","ciRNA_id_2","NM_ID","gene_symbol")
venn_in=  list( Esc_con_ID=Esc_con_small$ciRNA_id,  Esc_deci_ID=Esc_deci_small$ciRNA_id,
                known_ciRNAs=KnownciRNA$ciRNA_id_2)
library(VennDiagram)
venn.diagram(venn_in,"VennDiagram.tiff",col = "transparent",
             fill = c("skyblue", "pink1","darkorchid4" ),
             alpha = 0.50
)



#####           ciRNA_above_readNumber>1

keep=small_datalist2$iPSC$readNumber>1
iPSC_filtered=small_datalist2$iPSC[keep,]

keep=small_datalist2$CM$readNumber>1
CM_filtered=small_datalist2$CM[keep,]

venn_in=  list( iPSC=iPSC_filtered$gene_id,  CM=CM_filtered$gene_id, known_ciRNAs=small_known_ciRNA$gene_id)


library(VennDiagram)
venn.diagram(venn_in,"VennDiagram.tiff",col = "transparent",
             fill = c("skyblue", "pink1","darkorchid4" ),
             alpha = 0.50
             #label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
             # cat.col = c("skyblue", "pink1","mediumorchid", "orange")
             #c("darkblue", "darkgreen", "orange", "darkorchid4")
)





#####    pick up group specific genes    1-2-2017

iPSC=data.frame(small_datalist2$iPSC)
CM=data.frame(small_datalist2$CM)
known_ciRNAs=data.frame(small_known_ciRNA$gene_id)
names(known_ciRNAs)=c("gene_id")


### iPSC_CM 270
iPSC_CM=merge(x=iPSC,y=CM,by="gene_id")
iPSC_CM=iPSC_CM[,-c(3,5,6,7,8,10,11)]
names(iPSC_CM)=c("gene_id","gene_name","iPSC_rpb","CM_rpb")
write.table(iPSC_CM,"iPSC_CM_common.txt",sep="\t",row.names = FALSE)


### iPSC_know 
iPSC_know=merge(x=iPSC,y=known_ciRNAs,by="gene_id")
iPSC_know=iPSC_know[,-c(3,5,6)]
names(iPSC_know)=c("gene_id","gene_name","iPSC_rpb")
write.table(iPSC_know,"iPSC_known_common.txt",sep="\t",row.names = FALSE)


### CM_know
CM_know=merge(x=CM,y=known_ciRNAs,by="gene_id")
CM_know=CM_know[,-c(3,5,6)]
names(CM_know)=c("gene_id","gene_name","CM_rpb")
write.table(CM_know,"CM_known_common.txt",sep="\t",row.names = FALSE)

### all_common  226

all_common=merge(x=iPSC_CM,y=iPSC_know,by="gene_id")
all_common=all_common[,-c(5,6)]
names(all_common)=c("gene_id","gene_name","iPSC_rpb","CM_rpb")
write.table(all_common,"all_common.txt",sep="\t",row.names = FALSE)



### iPSC_CM specific  44
keep =!iPSC_CM$gene_id %in%  all_common$gene_id
iPSC_CM_specific=iPSC_CM[keep,]
write.table(iPSC_CM_specific,"iPSC_CM_specific.txt",sep="\t",row.names = FALSE)



### iPSC_all specific  728
keep =!iPSC_know$gene_id %in%  all_common$gene_id
iPSC_all_specific=iPSC_know[keep,]
write.table(iPSC_all_specific,"iPSC_all_specific.txt",sep="\t",row.names = FALSE)


### CM_all specific  1072
keep =!CM_know$gene_id %in%  all_common$gene_id
CM_all_specific=CM_know[keep,]
write.table(CM_all_specific,"CM_all_specific.txt",sep="\t",row.names = FALSE)


### iPSC_novel 614
keep =!iPSC$gene_id %in%  iPSC_CM$gene_id
iPSC_novel_1=iPSC[keep,]

keep1 =!iPSC_novel_1$gene_id %in%  iPSC_all_specific$gene_id
iPSC_novel_1=iPSC_novel_1[keep1,]
iPSC_novel_1=iPSC_novel_1[,-c(3,5,6)]
write.table(iPSC_novel_1,"iPSC_novel.txt",sep="\t",row.names = FALSE)


### CM_novel 2918
keep =!CM$gene_id %in%  iPSC_CM$gene_id
CM_novel_1=CM[keep,]

keep1 =!CM_novel_1$gene_id %in%  CM_all_specific$gene_id
CM_novel_1=CM_novel_1[keep1,]
CM_novel_1=CM_novel_1[,-c(3,5,6)]
write.table(CM_novel_1,"CM_novel.txt",sep="\t",row.names = FALSE)








############################### correlationship of number of ciRNAs and length of mRNA ###################################################
################                   3-17-2017
######################################==========>>>>>
setwd("F:/New_Genome_data/Hu Shijun/Circular RNA/CircRNA_V2/Circ2_all") ## windows
filenames=list("F1-EiPSC_nor.csv","iPSC_CM_nor.csv")  ## two groups 12-12-2016

iPSC=read.table("F1-EiPSC_nor.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#keep=iPSC$readNumber>=2
#iPSC=iPSC[keep,]

CM=read.table("iPSC_CM_nor.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
#keep=CM$readNumber>=2
#CM=CM[keep,]

iPSC_freq=as.data.frame(table(iPSC$geneName))
names(iPSC_freq)=c("gene_name","Freq")
iPSC_freq$flag=c("iPSC")

CM_freq=as.data.frame(table(CM$geneName))
names(CM_freq)=c("gene_name","Freq")
CM_freq$flag=c("CM")

exon_read=read.table("CM_iPSC_mRNA_rpkm_length_exon.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
exon_read$gene_name=row.names(exon_read)

iPSC_freq_merged=merge(x=iPSC_freq,y=exon_read,by="gene_name",all.x=TRUE)
CM_freq_merged=merge(x=CM_freq,y=exon_read,by="gene_name",all.x=TRUE)

# remove NA row
iPSC_merged_clean=iPSC_freq_merged[!is.na(iPSC_freq_merged$length),]
iPSC_merged_clean=iPSC_merged_clean[,c(1,2,3,9)]

CM_merged_clean=CM_freq_merged[!is.na(CM_freq_merged$length),]
CM_merged_clean=CM_merged_clean[,c(1,2,3,9)]

ciRNA_isofrom_and_exon_count=rbind(iPSC_merged_clean,CM_merged_clean)

library(ggplot2)
a=ggplot(ciRNA_isofrom_and_exon_count,aes(log(exon,2),Freq,color=flag))+geom_point(alpha=0.5)+xlim(0,8)+ylim(0,10)+
  #geom_smooth(colour = "red", fill = "lightgreen", method = 'lm') +
  #geom_boxplot()+
  ylab("Number of ciRNA isoforms") +xlab("Number of Exons (Log2)")


#a=ggplot(ciRNA_isofrom_and_exon_count,aes(Freq,exon,color=flag))+geom_point(alpha=0.5)+xlim(0,20)+ylim(0,100)+
#xlab("Number of ciRNA isoforms") +ylab("Number of Exons")
a+theme_set(theme_bw())
a+theme(panel.grid.major=element_line(colour=NA)) +  theme(panel.grid=element_blank()) 




##############################    ciRNA and mRNA   9-7-2017  ############################
################          
################           
###########################################==========>>>>>

setwd("F:/New_Genome_data/Hu Shijun/Circular RNA/CircRNA_V2/Circ2_all") ## windows
all_mRNA=read.table("CM_iPSC_mRNA_rpkm_edger.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
mRNA_small_list=all_mRNA[,c(3:4)]
mRNA_small_list$gene_name=all_mRNA$genes

all_ciRNA=read.table("all_4_in_one_v1.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
all_ciRNA=all_ciRNA[,-c(4,6)]
all_ciRNA[is.na(all_ciRNA)]=0
keep=rowSums(all_ciRNA[,4:5])>0
all_ciRNA=all_ciRNA[keep,]


mRNA2ciRNA=merge(x = all_ciRNA, y = mRNA_small_list, by = "gene_name", all.x = TRUE)
## remove genes such 1-sep et al.
keep1=is.na(mRNA2ciRNA[,7])!=TRUE
mRNA2ciRNA2=mRNA2ciRNA[keep1,]
names(mRNA2ciRNA2)=c("gene_name","gene_id","isoformName", "ciPSC","ciCM","CM","iPSC")

library(ggplot2)

#mRNA2ciRNA2$C2_ratio=(mRNA2ciRNA2[,4]+0.1)/(mRNA2ciRNA2[,7]+0.1)
mRNA2ciRNA2$CM_ratio=(mRNA2ciRNA2$ciCM+0.1)/(mRNA2ciRNA2$CM+0.1)
mRNA2ciRNA2$iPSC_ratio=(mRNA2ciRNA2$ciPSC+0.1)/(mRNA2ciRNA2$iPSC+0.1)

mRNA2ciRNA2$ration_fold=(mRNA2ciRNA2$CM_ratio+0.01)/(mRNA2ciRNA2$iPSC_ratio+0.01)
mRNA2ciRNA2$flag=cut(mRNA2ciRNA2$ration_fold,breaks=c(0,0.2,5,Inf),labels=c("-1","0","1"))

write.table(mRNA2ciRNA2,"mRNA2ciRNA2_ratio.txt",sep="\t")

p <- ggplot(mRNA2ciRNA2, aes(log(iPSC_ratio,2), log(CM_ratio,2)))
pic=p+geom_point(aes(color=flag))+xlab("Ratio of ciRNA/mRNA in iPSC")+ylab("Ratio of ciRNA/mRNA in CM")
pic
pic+theme(panel.grid=element_blank())  












##            alternative
setwd("F:/New_Genome_data/Hu Shijun/Circular RNA/CircRNA_V2/Circ2_all") ## windows
all_mRNA=read.table("C2_CM_iPSC_mRNA.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
mRNA_small_list=data.frame(all_mRNA[,5])
mRNA_small_list$gene_name=row.names(all_mRNA)
names(mRNA_small_list)=c("mRNA","gene_name")


CM_ciRNA=read.table("iPSC_CM_nor.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
CM_ciRNA_small=data.frame(cbind(paste(CM_ciRNA$chrom,CM_ciRNA$start,CM_ciRNA$end),CM_ciRNA$readNumber,CM_ciRNA$rpb,CM_ciRNA$geneName))
names(CM_ciRNA_small)=c("gene_id","readNumber","rpb","gene_name")
mRNA2ciRNA=merge(x = CM_ciRNA_small, y = mRNA_small_list, by = "gene_name", all.x = TRUE)
keep=!is.na(mRNA2ciRNA$mRNA)
mRNA2ciRNA=mRNA2ciRNA[keep,]
mRNA2ciRNA$readNumber=as.numeric(mRNA2ciRNA$readNumber)
mRNA2ciRNA$rpb=as.numeric(mRNA2ciRNA$rpb)
mRNA2ciRNA$CM_ratio=(mRNA2ciRNA[,4]+0.1)/(mRNA2ciRNA[,5]+0.1)
library(ggplot2)
p <- ggplot(mRNA2ciRNA, aes(rpb, log(mRNA)))
pic=p+geom_point()
pic


#all_ciRNA=read.table("all_4_in_one_v1.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE)
#all_ciRNA=all_ciRNA[,-6]
#all_ciRNA[is.na(all_ciRNA)]=0
#mRNA2ciRNA=merge(x = all_ciRNA, y = mRNA_small_list, by = "gene_name", all.x = TRUE)
## remove genes such 1-sep et al.
#keep1=is.na(mRNA2ciRNA[,7])!=TRUE
#mRNA2ciRNA2=mRNA2ciRNA[keep1,]
#names(mRNA2ciRNA2)=c("gene_name","gene_id","isoformName","ciC2", "ciPSC","ciCM","C2","CM","iPSC")



#mRNA2ciRNA2$C2_ratio=(mRNA2ciRNA2[,4]+0.1)/(mRNA2ciRNA2[,7]+0.1)

#mRNA2ciRNA2$iPSC_ratio=(mRNA2ciRNA2[,5]+0.1)/(mRNA2ciRNA2[,9]+0.1)


head(mRNA2ciRNA2)

p <- ggplot(mRNA2ciRNA2, aes(log(C2_ratio,2), log(CM_ratio,2)))
p <- ggplot(mRNA2ciRNA2, aes(log(C2_ratio,2), log(iPSC_ratio,2)))
p <- ggplot(mRNA2ciRNA2, aes(log(CM_ratio,2), log(iPSC_ratio,2)))
pic=p+geom_point()
pic

pic+theme_set(theme_bw())   # no backgroud
pic+theme(panel.grid.major=element_line(colour=NA))    # no grid

###     fold_change
###     
mRNA2ciRNA2$ciC2toCM_fold=(mRNA2ciRNA2[,4]+0.1)/(mRNA2ciRNA2[,6]+0.1)
mRNA2ciRNA2$C2toCM_fold=(mRNA2ciRNA2[,7]+0.1)/(mRNA2ciRNA2[,8]+0.1)
p <- ggplot(mRNA2ciRNA2, aes(log(ciC2toCM_fold,2), log(C2toCM_fold,2)))
pic=p+geom_point()
pic


