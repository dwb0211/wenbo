############  All_DNA_binding protein from MGI    #############

setwd("F:/New_Genome_data/All_TFs")

#DNA_binding=read.table("GO_term_summary_20180121_010817.csv", sep=",")

##    读之前人工去掉文本中的  引号
DNA_binding=read.table("GO_0003677_from_MGI_small.txt",header=T, sep="\t",quote="",stringsAsFactors = F)

TF=DNA_binding[,c(2:3)]
colnames(TF)=c("Gene_name","description")

TF1=TF[!duplicated(TF$Gene_name),]
#len=nrow(TF1)
#for (i in 1:len){
#  if(grepl("\"",TF1[i,2])){
#    m=unlist(strsplit(as.character(TF1[i,2]),"[\"]"))
#   # print(m)
#    TF1[i,2] = m[2]
#  }
#}


#a=unlist(strsplit(as.character(TF1[5,2]), "[\"]"))

write.table(TF1,"All_DNA_binding_protein_from_MGI.txt",sep="\t",row.names=T,quote=F)




###########     All_transcription factors from Jaspar     ########

PFAM=read.table("JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt", header=FALSE, sep="\t",stringsAsFactors=FALSE,fill=TRUE)
a=grep("^>",PFAM$V1)
b=PFAM[a,]

colnames(b)=c("MA","Gene_name")
len=nrow(b)
n_b=b
for (i in 1:len){
  if(grepl("::",n_b[i,2])){
      m=unlist(strsplit(n_b[i,2],"::"))
      n_b[i,2] = m[1]
      n_b[nrow(n_b)+1,1]= n_b[i,1]
      n_b[nrow(n_b),2]= m[2]
  }
}


lenbn=nrow(n_b)
for (i in 1:lenbn){
  if(grepl("\\(",n_b[i,2])){
    #print(c("b",i))
    aa=unlist(strsplit(n_b[i,2],"\\("))
    n_b[i,2] = aa[1]
  }
}


library(Hmisc)
n_b$Gene_name=capitalize(tolower(n_b$Gene_name))
n_b=n_b[!duplicated(n_b$Gene_name),]

write.table(n_b,"jaspar_nr.txt",sep="\t")




##################   for meme file ###################

PFAM=read.table("JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt", header=FALSE, sep="\t",stringsAsFactors=FALSE,fill=TRUE)
a=grep("^MOTIF",PFAM$V1)
b=data.frame(PFAM[a,])

colnames(b)=c("Gene_name")
len=nrow(b)

for (i in 1:len){
    m=unlist(strsplit(as.character(b[i,1])," "))
   b[i,2]=m[1]
   b[i,3]=m[2]
   b[i,4]=m[3]
}



n_b=b[,c(3:4)]

for (i in 1:len){
  if(grepl("::",n_b[i,2])){
    m=unlist(strsplit(n_b[i,2],"::"))
    n_b[i,2] = m[1]
    n_b[nrow(n_b)+1,1]= n_b[i,1]
    n_b[nrow(n_b),2]= m[2]
  }
}


lenbn=nrow(n_b)
for (i in 1:lenbn){
  if(grepl("\\(",n_b[i,2])){
    #print(c("b",i))
    aa=unlist(strsplit(n_b[i,2],"\\("))
    n_b[i,2] = aa[1]
  }
}

colnames(n_b)=c("MA","Gene_name")

library(Hmisc)
n_b$Gene_name=capitalize(tolower(n_b$Gene_name))
n_b=n_b[!duplicated(n_b$Gene_name),]
write.table(n_b,"meme_nr.txt",sep="\t")




##################   for transfac file ###################

PFAM=read.table("JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.txt", header=FALSE, sep="\t",stringsAsFactors=FALSE,fill=TRUE)
a=grep("^DE",PFAM$V1)
b=data.frame(PFAM[a,])

#colnames(b)=c("Gene_name")
len=nrow(b)

for (i in 1:len){
  m=unlist(strsplit(as.character(b[i,1])," "))
  b[i,2]=m[2]
  b[i,3]=m[3]

}


n_b=b[,c(2:3)]

for (i in 1:len){
  if(grepl("::",n_b[i,2])){
    m=unlist(strsplit(n_b[i,2],"::"))
    n_b[i,2] = m[1]
    n_b[nrow(n_b)+1,1]= n_b[i,1]
    n_b[nrow(n_b),2]= m[2]
  }
}


lenbn=nrow(n_b)
for (i in 1:lenbn){
  if(grepl("\\(",n_b[i,2])){
    #print(c("b",i))
    aa=unlist(strsplit(n_b[i,2],"\\("))
    n_b[i,2] = aa[1]
  }
}

colnames(n_b)=c("MA","Gene_name")

library(Hmisc)
n_b$Gene_name=capitalize(tolower(n_b$Gene_name))
n_b=n_b[!duplicated(n_b$Gene_name),]
write.table(n_b,"transfac_nr.txt",sep="\t")






###############   merge all TF into a single file   ########

filenames=list("jaspar_nr.txt","jaspar_r.txt","meme_nr.txt","meme_r.txt","transfac_nr.txt","transfac_nr.txt")


mydatalist <- lapply(filenames, function(name) {
  read.table(name, sep = "\t",header=TRUE)
})



a=mydatalist[[1]]
for(i in 1:length(mydatalist)){
  a=merge(mydatalist[[i]],a,by="Gene_name",all =TRUE)
}

g=a[!duplicated(a$Gene_name),]

## equal to above
a=merge(mydatalist[[1]],mydatalist[[2]],by="Gene_name",all=TRUE)
b=merge(mydatalist[[3]],mydatalist[[4]],by="Gene_name",all=TRUE)
c=merge(mydatalist[[5]],mydatalist[[6]],by="Gene_name",all=TRUE)
d=merge(a,b,by="Gene_name",all=TRUE)
e=merge(d,c,by="Gene_name",all=TRUE)
f=e[!duplicated(e$Gene_name),]

write.table(g[,1:2],"all_tfs.txt",sep="\t",row.names = F,quote=F)








#########     to a single file      ######

all_DNAbp=read.table("All_DNA_binding_protein_from_MGI.txt", sep="\t",header=TRUE,quote="")
#All_DNA_binding_protein_from_MGI.txt
#write.table(all_DNAbp,"a.txt",sep="\t")

all_tf=read.table("all_tfs.txt", sep = "\t",header=TRUE)
colnames(all_tf)=c("Gene_name","description")

total_bp=merge(all_DNAbp,all_tf,by="Gene_name",all=TRUE)

#total_bp=total_bp[!duplicated(total_bp$Gene_name),]
a=total_bp
write.table(a,"DNAbinding_plus_TF_db.txt",sep="\t")



