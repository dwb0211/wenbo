setwd("F:/New_Genome_data/Single_cell/Mapping the Mouse Cell Atlas by Microwell-Seq")
library(Rtsne)


exp1=read.table("GSM2906478_Uterus1_dge.txt",header=TRUE,sep="")
colnames(exp1)=c(paste("cell",c(1:ncol(exp1)),sep=""))

epx2=read.table("GSM2906479_Uterus2_dge.txt",header=TRUE,sep="")
exp2=epx2
colnames(exp2)=c(paste("cell",c(2047:(2047+ncol(exp2)-1)),sep=""))

exp1$gene_name=row.names(exp1)
exp2$gene_name=row.names(exp2)
exp=merge(exp1,exp2,by="row.names",all=TRUE)
row.names(exp)=exp$Row.names
exp=exp[,-c(1,2048,3764)]
exp[is.na(exp)] = 0

#expa=exp[,1:500]
library(Matrix)

B <- Matrix(as.matrix(exp), sparse = TRUE)   

#B <- Matrix(as.matrix(epx2), sparse = TRUE)   
#expa=exp[,1:500]
#d <- dist(t(expa))



d[is.na(d)] = 0



setwd("F:/New_Genome_data/Single_cell/test/")

#install.packages('Seurat')

library(Seurat)
library(dplyr)
library(Matrix)

#pbmc.data <- Read10X(data.dir = "hg19/")


pbmc.data=B

dense.size <- object.size(x = as.matrix(x = pbmc.data))
sparse.size <- object.size(x = pbmc.data)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")


mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = pbmc, PCs = 1:12)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 1, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)

PrintFindClustersParams(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, max_iter = 2000, do.fast = FALSE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)
save(pbmc, file = "uterus_2.Robj")

cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), 
                                min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))



pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

write.table(pbmc.markers,"cell_types.txt",sep="\t")



save(pbmc, file = "uterus_2.Robj")





pbmc.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)



#cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(object = pbmc, features.plot = c("Pgr"))

VlnPlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Cxcl15","Pgr","Hoxa10","Ccl11","Sprr2f","Mif","Dio2","Dlx5"), use.raw = TRUE, y.log=TRUE)

FeaturePlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Cxcl15","Pgr","C3","Ccl11","Sprr2f","Mif","Dio2","Dlx5","Tlr4","Cxcl2","Ccl4","Ncr1","Adgre1","Pdgfra"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


#####  Macrophage
FeaturePlot(object = pbmc, features.plot = c("Adgre1","Csf1r","Cd74","C1qc","Cd68"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


### NK cells
FeaturePlot(object = pbmc, features.plot = c("Ncr1","Klra4","Xcl1"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


### Endothelium cells
FeaturePlot(object = pbmc, features.plot = c("Esam","Egfl7","Emcn","Flt1","Pecam1"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

### DC cells
FeaturePlot(object = pbmc, features.plot = c("Cd209a","Ccr7","Il1r2"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")



FeaturePlot(object = pbmc, features.plot = c("Hbegf"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")



 FeaturePlot(object = pbmc, features.plot = c("Tlr2"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
VlnPlot(object = pbmc, features.plot = c("Cxcl15"), use.raw = TRUE)







top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)

pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.8, print.output = FALSE)
## Warning in BuildSNN(object = object, genes.use = genes.use, reduction.type
## = reduction.type, : Build parameters exactly match those of already
## computed and stored SNN. To force recalculation, set force.recalc to TRUE.
# Demonstration of how to plot two tSNE plots side by side, and how to color
# points based on different criteria
plot1 <- TSNEPlot(object = pbmc, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
plot2 <- TSNEPlot(object = pbmc, do.return = TRUE, group.by = "ClusterNames_0.6", 
                  no.legend = TRUE, do.label = TRUE)
plot_grid(plot1, plot2)

tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
# can see that CCR7 is upregulated in C0, strongly indicating that we can
# differentiate memory from naive CD4 cells.  cols.use demarcates the color
# palette from low to high expression
FeaturePlot(object = pbmc, features.plot = c("S100A4", "CCR7"), cols.use = c("green", 
                                                                             "blue"))


pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
save(pbmc, file = "~/Projects/datasets/pbmc3k_final.Rda")






















###### ======================
tsne_out <- Rtsne(d, is_distance=TRUE, perplexity=10, verbose = TRUE) 

plot(tsne_out$Y, pch=16, main='tSNE')



a=data.frame(tsne_out$Y)
p <- ggplot(a,aes(a[,1],a[,2]))
pic=p+geom_point()







iris_unique <- unique(iris) # Remove duplicates
iris_matrix <- as.matrix(iris_unique[,1:4])
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(iris_matrix) # Run TSNE

# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col=iris_unique$Species)

# Using a dist object
tsne_out <- Rtsne(dist(iris_matrix))
plot(tsne_out$Y,col=iris_unique$Species)















d <- dist(t(exp))

tsne_out <- Rtsne(d, is_distance=TRUE, perplexity=10, verbose = TRUE) 

library(ggplot2)



a=data.frame(tsne_out$Y)
p <- ggplot(a,aes(a[,1],a[,2]))
pic=p+geom_point()
pic

+
  coord_cartesian(xlim=c(0,15),ylim=c(0,15))+ 
  scale_color_discrete(h=c(100,350), c=100, l=60)+
  xlab("WT_LE")+ylab("WT_ST") +
  geom_text_repel(
    data = SC_LE[SC_LE$Gene_name %in% c("Sox17","Klf5","Dlx5","Msx1","Dlx6","Pax2",  "Hoxa10","Hoxa11","Myc","WT1","Sox18"),],
    aes(label = Gene_name),
    size = 5,
    box.padding = unit(1, "lines"), 
    point.padding = unit(0.1, "lines"))



pic +theme(panel.grid.major=element_line(colour=NA)) +theme(panel.grid=element_blank())







plot(tsne_out$Y, pch=16, main='tSNE')


d <- dist(t(mat))
set.seed(0) # tsne has some stochastic steps (gradient descent) so need to set random 
tsne_out <- Rtsne(d, is_distance=TRUE, perplexity=10, verbose = TRUE) 


















######
library(ggplot2)
library(devtools)
library(scMCA)
??scMCA
scMCA(exp, plot_heat = TRUE, interactive_plot = TRUE, numbers_plot = 3)




















#######################   Epi    #########################

setwd("F:/New_Genome_data/Single_cell/test")

library(Rtsne)


exp=read.table("GSE98451_uterus_single_cell_RNA-Seq_counts.txt",header=TRUE,sep="\t")
exp=exp[!duplicated(exp$gene.name.sample.name),]
row.names(exp)=exp[,1]
exp=exp[,-1]



#exp[is.na(exp)] = 0

#expa=exp[,1:500]
library(Matrix)

B <- Matrix(as.matrix(exp), sparse = TRUE)   

library(Seurat)
library(dplyr)
library(Matrix)

#pbmc.data <- Read10X(data.dir = "hg19/")


pbmc.data=B

dense.size <- object.size(x = as.matrix(x = pbmc.data))
sparse.size <- object.size(x = pbmc.data)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")


mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = pbmc, PCs = 1:12)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)

save(pbmc, file = "uterus_2.Robj")

cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), 
                                min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))



pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

write.table(pbmc.markers,"Epi_cell_types.txt",sep="\t")



save(pbmc, file = "uterus_2.Robj")





pbmc.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)



#cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(object = pbmc, features.plot = c("Pgr"))

VlnPlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Cxcl15","Pgr","Hoxa10","Ccl11","Sprr2f","Mif","Dio2","Dlx5"), use.raw = TRUE, y.log=TRUE)

FeaturePlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Cxcl15","Pgr","C3","Ccl11","Sprr2f","Mif","Dio2","Dlx5","Tlr4","Cxcl2","Ccl4","Ncr1","Adgre1","Pdgfra"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


#####  Macrophage
FeaturePlot(object = pbmc, features.plot = c("Meg3"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


### NK cells
FeaturePlot(object = pbmc, features.plot = c("Ncr1","Klra4","Xcl1"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


### Endothelium cells
FeaturePlot(object = pbmc, features.plot = c("Esam","Egfl7","Emcn","Flt1","Pecam1"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

### DC cells
FeaturePlot(object = pbmc, features.plot = c("Cd209a","Ccr7","Il1r2"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")



FeaturePlot(object = pbmc, features.plot = c("Wnt7","Pcp4","Dcn","Vim"), cols.use = c("grey", "Blue"), 
            reduction.use = "tsne")



FeaturePlot(object = pbmc, features.plot = c("Tlr2"), cols.use = c("grey", "Red"), 
            reduction.use = "tsne")
VlnPlot(object = pbmc, features.plot = c("Cxcl15"), use.raw = TRUE)







top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)

pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.8, print.output = FALSE)
## Warning in BuildSNN(object = object, genes.use = genes.use, reduction.type
## = reduction.type, : Build parameters exactly match those of already
## computed and stored SNN. To force recalculation, set force.recalc to TRUE.
# Demonstration of how to plot two tSNE plots side by side, and how to color
# points based on different criteria
plot1 <- TSNEPlot(object = pbmc, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
plot2 <- TSNEPlot(object = pbmc, do.return = TRUE, group.by = "ClusterNames_0.6", 
                  no.legend = TRUE, do.label = TRUE)
plot_grid(plot1, plot2)

tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
# can see that CCR7 is upregulated in C0, strongly indicating that we can
# differentiate memory from naive CD4 cells.  cols.use demarcates the color
# palette from low to high expression
FeaturePlot(object = pbmc, features.plot = c("S100A4", "CCR7"), cols.use = c("green", 
                                                                             "blue"))


pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
save(pbmc, file = "~/Projects/datasets/pbmc3k_final.Rda")





####    placenta   3-5-2018   #####

setwd("F:/New_Genome_data/Single_cell/test/placenta/")
library(Rtsne)

exp1=read.table("GSM2906465_PlacentaE14.1_dge.txt",header=TRUE,sep="")
colnames(exp1)=c(paste("cell",c(1:ncol(exp1)),sep=""))
epx2=read.table("GSM2906466_PlacentaE14.2_dge.txt",header=TRUE,sep="")
exp2=epx2
colnames(exp2)=c(paste("cell",c((ncol(exp1)+1):(ncol(exp1)+ncol(exp2))),sep=""))

exp1$gene_name=row.names(exp1)
exp2$gene_name=row.names(exp2)
exp=merge(exp1,exp2,by="row.names",all=TRUE)
row.names(exp)=exp$Row.names
exp=exp[,-c(1,ncol(exp1)+1,ncol(exp))]
exp[is.na(exp)] = 0

#expa=exp[,1:500]
library(Matrix)
B <- Matrix(as.matrix(exp), sparse = TRUE)   
#setwd("F:/New_Genome_data/Single_cell/test/")
#install.packages('Seurat')

library(Seurat)
library(dplyr)
library(Matrix)

#pbmc.data <- Read10X(data.dir = "hg19/")


pbmc.data=B

dense.size <- object.size(x = as.matrix(x = pbmc.data))
sparse.size <- object.size(x = pbmc.data)

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "placenta")


#mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
#percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
#pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
#VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:30, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
#pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(object = pbmc, PCs = 1:12)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:20, 
                     resolution = 4, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:20, do.fast = TRUE, perplexity=30,max_iter=10000)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)
save(pbmc, file = "placenta_2.Robj")

cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), 
                                min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))



pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

write.table(pbmc.markers,"Placenta_cell_types.txt",sep="\t")



                                                                                                                                                                                                                                                                                    save(pbmc, file = "uterus_2.Robj")





pbmc.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)



#cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(object = pbmc, features.plot = c("Prl3b1"))

VlnPlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Arg1","Pgr","Hoxa10","Ccl11","Sprr2f","Mif","Ctss","Flt1"), use.raw = TRUE, y.log=TRUE)

FeaturePlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Cxcl15","Pgr","C3","Ccl11","Sprr2f","Mif","Dio2","Dlx5","Tlr4","Cxcl2","Ccl4","Ncr1","Adgre1","Pdgfra"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")



FeaturePlot(object = pbmc, features.plot = c("Prl3b1"), cols.use = c("grey", "red"), 
            reduction.use = "tsne")











#######################     batch removed  3-13-2018  ############# 





setwd("F:/New_Genome_data/Single_cell/test/placenta/")
library(Rtsne)

exp1=read.table("PlacentaE14.1_rm.batch_dge.txt",header=TRUE,sep="")
colnames(exp1)=c(paste("cell",c(1:ncol(exp1)),sep=""))
exp2=read.table("PlacentaE14.2_rm.batch_dge.txt",header=TRUE,sep="")
colnames(exp2)=c(paste("cell",c((ncol(exp1)+1):(ncol(exp1)+ncol(exp2))),sep=""))

exp1$gene_name=row.names(exp1)
exp2$gene_name=row.names(exp2)
exp=merge(exp1,exp2,by="row.names",all=TRUE)
row.names(exp)=exp$Row.names
exp=exp[,-c(1,ncol(exp1)+1,ncol(exp))]
exp[is.na(exp)] = 0
library(Matrix)
B <- Matrix(as.matrix(exp), sparse = TRUE)   
#setwd("F:/New_Genome_data/Single_cell/test/")
#install.packages('Seurat')

library(Seurat)
library(dplyr)
library(Matrix)
#pbmc.data <- Read10X(data.dir = "hg19/")
pbmc.data=B
dense.size <- object.size(x = as.matrix(x = pbmc.data))
sparse.size <- object.size(x = pbmc.data)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "placenta")


#mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
#percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
#pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
#VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
#PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
#pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(object = pbmc, PCs = 1:12)

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:20, 
                     resolution = 1, print.output = 0, save.SNN = TRUE,force.recalc =TRUE)

PrintFindClustersParams(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:20, do.fast = TRUE, perplexity=30,max_iter=1000)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)
save(pbmc, file = "placenta_2.Robj")

cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), 
                                min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))



pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

write.table(pbmc.markers,"Placenta_cell_types.txt",sep="\t")



save(pbmc, file = "uterus_2.Robj")





pbmc.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)



#cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(object = pbmc, features.plot = c("Prl3b1"))

VlnPlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Arg1","Pgr","Hoxa10","Ccl11","Sprr2f","Mif","Ctss","Flt1"), use.raw = TRUE, y.log=TRUE)

FeaturePlot(object = pbmc, features.plot = c("Msx1", "Foxa2","Col1a2","Cxcl14","Cxcl15","Pgr","C3","Ccl11","Sprr2f","Mif","Dio2","Dlx5","Tlr4","Cxcl2","Ccl4","Ncr1","Adgre1","Pdgfra"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")



FeaturePlot(object = pbmc, features.plot = c("Prl3b1"), cols.use = c("grey", "red"), 
            reduction.use = "tsne")


