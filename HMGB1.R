setwd("F:/New_Genome_data/D4_6_8")


Exp=read.table("D4_6_8_rpkm.txt", head=TRUE, sep="\t")


#HMG=subset(Exp,Exp$Gene %in% "Hmg")

HMG=subset(Exp,Exp$X %in% grep("^Hmg", Exp$X, value = TRUE))
D4_rep1=HMG[,c(1,2)]
colnames(D4_rep1)=c("Symbol","Exp")
D4_rep2=HMG[,c(1,6)]
colnames(D4_rep2)=c("Symbol","Exp")
D4_exp=rbind(D4_rep1,D4_rep2)
D4_exp$flag=c("D4")

D5_rep1=HMG[,c(1,3)]
colnames(D5_rep1)=c("Symbol","Exp")
D5_rep2=HMG[,c(1,7)]
colnames(D5_rep2)=c("Symbol","Exp")
D5_exp=rbind(D5_rep1,D5_rep2)
D5_exp$flag=c("D5")

D8_rep1=HMG[,c(1,4)]
colnames(D8_rep1)=c("Symbol","Exp")
D8_rep2=HMG[,c(1,8)]
colnames(D8_rep2)=c("Symbol","Exp")
D8_exp=rbind(D8_rep1,D8_rep2)
D8_exp$flag=c("D8")


Exp_all=rbind(D4_exp,D5_exp,D8_exp)

library(reshape2)
library(ggplot2)

p <- ggplot(Exp_all, aes(Symbol, Exp, fill=flag)) + 
  geom_boxplot() +stat_boxplot(geom ='errorbar')+xlab("") +ylab("HMG expression in day 4 Deci (FPKM)")


p+theme(panel.grid.major=element_line(colour=NA)) +theme(panel.grid=element_blank())










