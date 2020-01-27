#load packages
library("ggplot2")
library("DropletUtils") 
library("edgeR") 
library("dplyr")
library("Seurat") 
library("DoubletFinder")
library("sctransform")
library("data.table")
library("conos") 
library("pagoda2") 
library("MAST") 
library("graphics")

#load in source file with functions
source("/mnt/data/Maya/SHvsPH/Analysis/SHvsPH/single_cell_function.R")

#set working directory for plots to be saved
setwd("/mnt/data/Maya/SHvsPH/Analysis/SHvsPH/")

#define samples
sample_id_1="HS"
sample_id_2="PH"

#load in raw data matrices (outputs of CellRanger)
dataMatrix.1 = Load.Data.Matrices(list(SAMPLE=paste("/mnt/data/Maya/SHvsPH/Analysis/SHvsPH/",sample_id_1 ,"/outs/raw_feature_bc_matrix",sep="")))
dataMatrixDgC.1=dataMatrix.1[[1]]
dataMatrix.2 = Load.Data.Matrices(list(SAMPLE=paste("/mnt/data/Maya/SHvsPH/Analysis/SHvsPH/",sample_id_2 ,"/outs/raw_feature_bc_matrix",sep="")))
dataMatrixDgC.2=dataMatrix.2[[1]]

#Barcode ranks plot to find inflection and knee point
pdf("SH_barcodeRanks.pdf")
brRanks.1 <- barcodeRanks(dataMatrixDgC.1)
options(repr.plot.width=5, repr.plot.height=4)
plot(brRanks.1$rank, brRanks.1$total, log="xy", xlab="Barcodes Ranked", ylab="Total UMI per Barcode", pch=19, cex=0.3, main="SH")
ordered.1 <- order(brRanks.1$rank)
abline(h=metadata(brRanks.1)$knee, col="dodgerblue", lty=2)
abline(h=metadata(brRanks.1)$inflection, col="forestgreen", lty=2)
legend("topright", lty=2, col=c("dodgerblue", "forestgreen"),legend=c("knee", "inflection"), cex=1, pt.cex=0.1)
dev.off()

pdf("PH_barcodeRanks.pdf")
brRanks.2 <- barcodeRanks(dataMatrixDgC.2)
options(repr.plot.width=5, repr.plot.height=4)
plot(brRanks.2$rank, brRanks.2$total, log="xy", xlab="Barcodes Ranked", ylab="Total UMI per Barcode", pch=19, cex=0.3, main="PH")
ordered.2 <- order(brRanks.2$rank)
abline(h=metadata(brRanks.2)$knee, col="dodgerblue", lty=2)
abline(h=metadata(brRanks.2)$inflection, col="forestgreen", lty=2)
legend("topright", lty=2, col=c("dodgerblue", "forestgreen"),legend=c("knee", "inflection"), cex=1, pt.cex=0.1)
dev.off()

knee.1=brRanks.1@metadata$knee
inflection.1=brRanks.1@metadata$inflection
paste(inflection.1, knee.1)

knee.2=brRanks.2@metadata$knee
inflection.2=brRanks.2@metadata$inflection
paste(inflection.2, knee.2)

#emptydrops to filter out non cells
set.seed(1234)
emptyDropsOutput.1= emptyDrops(dataMatrixDgC.1, lower=inflection.1,retain=knee.1)
emptyDropsOutput.2 = emptyDrops(dataMatrixDgC.2, lower=inflection.2, retain=knee.2)

emptyDropsKeep.1 = emptyDropsOutput.1$FDR <= 0.01
emptyDropsKeep.na.1 = emptyDropsKeep.1
emptyDropsKeep.na.1[is.na(emptyDropsKeep.1)]= FALSE

emptyDropsKeep.2 = emptyDropsOutput.2$FDR <= 0.01
emptyDropsKeep.na.2 = emptyDropsKeep.2
emptyDropsKeep.na.2[is.na(emptyDropsKeep.2)]= FALSE

table(emptyDropsKeep.na.1)
table(emptyDropsKeep.na.2)

counts.1=dataMatrixDgC.1[,emptyDropsKeep.na.1]
counts.2=dataMatrixDgC.2[,emptyDropsKeep.na.2]

countsPerCell.1 = Matrix::colSums(counts.1)
countsGerGene.1 = Matrix::rowSums(counts.1)
genesPerCell.1 = Matrix::colSums(counts.1>0)
countsPerCell.2 = Matrix::colSums(counts.2)
countsGerGene.2 = Matrix::rowSums(counts.2)
genesPerCell.2 = Matrix::colSums(counts.2>0)

#histogram of counts per cell
pdf("SH_hist_countsPerCell.pdf")
hist(countsPerCell.1,breaks=42, col="powderblue", border="black", main="SH", xlab="Counts Per Cell", xlim=c(0, 25000), ylim=c(0, 17000), cex.lab=1, cex.axis=0.75, las=1)
dev.off()
paste("Median Counts per Cell = ",median(countsPerCell.1))
pdf("PH_hist_countsPerCell.pdf")
hist(countsPerCell.2,breaks=50, col="powderblue", border="black", main="PH", xlab="Counts Per Cell", xlim=c(0, 25000), ylim=c(0, 19000), cex.lab=1, cex.axis=0.75, las=1)
paste("Median Counts per Cell = ",median(countsPerCell.2))
dev.off()

#histogram of genes per cell
pdf("SH_hist_genesPerCell.pdf")
hist(genesPerCell.1,breaks=45,col="powderblue", border="black",main="SH", xlab="Genes Per Cell", xlim=c(0, 6000), ylim=c(0, 10000), cex.lab=1, cex.axis=0.75, las=1)
paste("Median Genes per Cell = ",median(genesPerCell.1))
dev.off()
pdf("PH_hist_genesPerCell.pdf")
hist(genesPerCell.2,breaks=50,col="powderblue", border="black",main="PH", xlab="Genes Per Cell", xlim=c(0, 6000), ylim=c(0, 13000), cex.lab=1, cex.axis=0.75, las=1)
paste("Median Genes per Cell = ",median(genesPerCell.2))
dev.off()

#genes vs counts
pdf("SH_genesVScells.pdf")
plot(countsPerCell.1, genesPerCell.1, main="SH", xlab="Counts Per Cell", ylab="Genes Per Cell",pch=19, cex=0.01,  cex.lab=1, cex.axis=0.75, las=1)
dev.off()
pdf("PH_genesVScells.pdf")
plot(countsPerCell.2, genesPerCell.2, main="PH", xlab="Counts Per Cell", ylab="Genes Per Cell",pch=19, cex=0.01,  cex.lab=1, cex.axis=0.75, las=1)
dev.off()

countsMx.1=as.matrix(counts.1)
countsMX.2=as.matrix(counts.2)

#create seurat object
SH.Sample= CreateSeuratObject(countsMx.1, min.cells=5, min.features=5, project = "Control")
PH.Sample= CreateSeuratObject(countsMX.2, min.cells=5, min.features=5, project = "Epilepsy")

#no of genes and cells per sample
dim(SH.Sample)
dim(PH.Sample)

#mitochondrial filtering
SH.Sample[["MT.PER"]] = PercentageFeatureSet(SH.Sample, pattern = "mt-")
PH.Sample[["MT.PER"]] = PercentageFeatureSet(PH.Sample, pattern = "mt-")

pdf("SH_hist_MTpercentage.pdf")
hist(SH.Sample@meta.data$MT.PER, breaks=45,col="powderblue", border="black", xlab="% Mitochondrial Genes", main="SH", cex.lab=1, cex.axis=0.75, las=1)
dev.off()
pdf("PH_hist_MTpercentage.pdf")
hist(PH.Sample@meta.data$MT.PER, breaks=45,col="powderblue", border="black", xlab="% Mitochondrial Genes", main="PH", cex.lab=1, cex.axis=0.75, las=1)
dev.off()

pdf("SH_MTvsUMI.pdf")
FeatureScatter(SH.Sample, feature1 = "nCount_RNA", feature2 = "MT.PER", cols="black", pt.size=0.1)
dev.off()
pdf("PH_MTvsUMI.pdf")
FeatureScatter(SH.Sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols="black", pt.size=0.1)
dev.off()

SH.MTfilt = subset(SH.Sample, subset = MT.PER < 10)
PH.MTfilt = subset(PH.Sample, subset = MT.PER < 10)

dim(SH.MTfilt)
dim(PH.MTfilt)

SH.MTfilt.counts=GetAssayData(object = SH.MTfilt, assay = "RNA", slot = "counts")
PH.MTfilt.counts=GetAssayData(object = PH.MTfilt, assay = "RNA", slot = "counts")

SH.filtSeuratOb= CreateSeuratObject(counts=as.matrix(SH.MTfilt.counts), min.cells=5, min.features=5, project = "Epilepsy")
SH.filtSeuratOb[["MT.PER"]] = PercentageFeatureSet(SH.filtSeuratOb, pattern = "mt-")
PH.filtSeuratOb= CreateSeuratObject(counts=as.matrix(PH.MTfilt.counts), min.cells=5, min.features=5, project = "Epilepsy")
PH.filtSeuratOb[["MT.PER"]] = PercentageFeatureSet(PH.filtSeuratOb, pattern = "mt-")

SH.FilteredCountsMatrix=GetAssayData(object = SH.filtSeuratOb, assay = "RNA", slot = "counts")
PH.FilteredCountsMatrix=GetAssayData(object = PH.filtSeuratOb, assay = "RNA", slot = "counts")

SH.MTfilt.data=GetAssayData(object = SH.filtSeuratOb, assay = "RNA", slot = "data")
PH.MTfilt.data=GetAssayData(object = SH.filtSeuratOb, assay = "RNA", slot = "data")

#histogram of expression before normalisation
pdf("SH_hist_totalExpr.pdf")
hist(colSums(SH.MTfilt.data),breaks = 100, main = "SH - Total expression before normalisation", xlab = "Counts Per Cell")
dev.off()
pdf("PH_hist_totalExpr.pdf")
hist(colSums(PH.MTfilt.data),breaks = 100, main = "PH - Total expression before normalisation", xlab = "Counts Per Cell")
dev.off()

#normalise
SH.filtNormData=NormalizeData(SH.MTfilt, normalization.method = "LogNormalize", scale.factor = 10000)
PH.filtNormData= NormalizeData(PH.MTfilt, normalization.method = "LogNormalize", scale.factor = 10000)

#histogram of expression before normalisation
SH.MTfilt.data=GetAssayData(object = SH.filtNormData, assay = "RNA", slot = "data")
PH.MTfilt.data=GetAssayData(object = SH.filtNormData, assay = "RNA", slot = "data")
pdf("SH_hist_normedExpr.pdf")
hist(colSums(SH.filtNormData),breaks = 100, main = "Total expression after normalisation", xlab = "Normalised Counts Per Cell")
dev.off()
pdf("PH_hist_normedExpr.pdf")
hist(colSums(PH.filtNormData),breaks = 100, main = "Total expression after normalisation", xlab = "Normalised Counts Per Cell")
dev.off()

#find variable features
SH.filtNormData = FindVariableFeatures(SH.filtNormData, selection.method = "mean.var.plot", verbose=TRUE)
PH.filtNormData = FindVariableFeatures(PH.filtNormData, selection.method = "mean.var.plot", verbose=TRUE)
#look into other parameters i.e. 'mean.function', 'dispersion.function', mean.cutoff, dispersion.cutoff

#scale data 
allGenes.SH = rownames(SH.filtNormData)
allGenes.PH = rownames(PH.filtNormData)
SH.filtNormScaled = ScaleData(SH.filtNormData, features = allGenes.SH)
PH.filtNormScaled = ScaleData(PH.filtNormData, features = allGenes.PH)
#can also scale based on just the variable features

#PCA
SH.filtNormScaledPCA = RunPCA(SH.filtNormScaled,npcs = 30, verbose=F)
PH.filtNormScaledPCA = RunPCA(PH.filtNormScaled,npcs = 30, verbose=F)
#can access PCA list using PrintPCAParams/PrintPCA

pdf("SH_PCAplot.pdf")
DimPlot(object = SH.filtNormScaledPCA, dim.1 = 1, dim.2 = 2)
dev.off()

pdf("PH_PCAplot.pdf")
DimPlot(object = PH.filtNormScaledPCA, dim.1 = 1, dim.2 = 2)
dev.off()

SH.filtNormScaledPCA = JackStraw(SH.filtNormScaledPCA, num.replicate = 200)
SH.filtNormScaledPCA = ScoreJackStraw(SH.filtNormScaledPCA, reduction = "pca", dims = 1:20)
pdf("SH_JackSrawPlot.pdf")
JackStrawPlot(object = SH.filtNormScaledPCA, dims= 1:20)
dev.off()

PH.filtNormScaledPCA = JackStraw(PH.filtNormScaledPCA, num.replicate = 200)
PH.filtNormScaledPCA = ScoreJackStraw(PH.filtNormScaledPCA, reduction = "pca", dims = 1:20)
pdf("PH_JackSrawPlot.pdf")
JackStrawPlot(object = PH.filtNormScaledPCA, dims= 1:20)
dev.off()

saveRDS(SH.filtNormScaledPCA, "SH_filtNormScaledPCA.rds")
saveRDS(PH.filtNormScaledPCA, "PH_filtNormScaledPCA.rds")

SH.filtNormScaledPCA=readRDS("/mnt/data/Maya/SHvsPH/Analysis/SHvsPH/SH_filtNormScaledPCA.rds")
PH.filtNormScaledPCA=readRDS("/mnt/data/Maya/SHvsPH/Analysis/SHvsPH/PH_filtNormScaledPCA.rds")

#clustering
SH.filtNormScaledPCA = FindNeighbors(SH.filtNormScaledPCA, dims = 1:20)
SH.filtNormScaledPCAclust = FindClusters(SH.filtNormScaledPCA, resolution = 1.2, pc.use=1:20)
SH.filtNormScaledPCAclust = RunUMAP(SH.filtNormScaledPCAclust, dims = 1:20)

PH.filtNormScaledPCA = FindNeighbors(PH.filtNormScaledPCA, dims = 1:20)
PH.filtNormScaledPCAclust = FindClusters(PH.filtNormScaledPCA, resolution = 1.2, pc.use=1:20)
PH.filtNormScaledPCAclust = RunUMAP(PH.filtNormScaledPCAclust, dims = 1:20)

pdf("SH_UMAP.pdf")
DimPlot(SH.filtNormScaledPCAclust, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()

pdf("PH_UMAP.pdf")
DimPlot(PH.filtNormScaledPCAclust, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()

#doubletFinder
Con.res.list.SH = paramSweep_v3(SH.filtNormScaledPCAclust, PCs = 1:20)
Con.res.list.stats.SH = summarizeSweep(Con.res.list.SH, GT = FALSE)
bcmvn.SH= find.pK(Con.res.list.stats.SH)
as.numeric(as.vector(bcmvn.SH$pK[which.max(bcmvn.SH$BCmetric)]))

Con.res.list.PH = paramSweep_v3(PH.filtNormScaledPCAclust, PCs = 1:20)
Con.res.list.stats.PH = summarizeSweep(Con.res.list.PH, GT = FALSE)
bcmvn.PH= find.pK(Con.res.list.stat.PH)
as.numeric(as.vector(bcmvn.PH$pK[which.max(bcmvn.PH$BCmetric)]))

annotations.SH = SH.filtNormScaledPCAclust@meta.data$RNA_snn_res.1.2
Prop.Homotypic.SH = modelHomotypic(annotations.SH)
nExp_poi.SH = round(0.07*length(SH.filtNormScaledPCAclust@active.ident))
nExp_poi.adj.SH= round(nExp_poi.SH*(1-Prop.Homotypic.SH))

annotations.PH = PH.filtNormScaledPCAclust@meta.data$RNA_snn_res.1.2
Prop.Homotypic.PH = modelHomotypic(annotations.PH)
nExp_poi.PH = round(0.07*length(PH.filtNormScaledPCAclust@active.ident))
nExp_poi.adj.PH = round(nExp_poi.PH*(1-Prop.Homotypic.PH))

SH.filtNormScaledPCAclust = doubletFinder_v3(SH.filtNormScaledPCAclust, PCs = 1:20, pN = 0.25, pK = pK.SH, nExp = nExp_poi.SH, reuse.pANN = FALSE,sct = FALSE)
SH.filtNormScaledPCAclust = doubletFinder_v3(SH.filtNormScaledPCAclust, PCs = 1:20, pN = 0.25, pK = pK.SH, nExp = nExp_poi.adj.SH, reuse.pANN = paste("pANN_",nExp_poi.SH,sep=""), sct = FALSE)
PH.filtNormScaledPCAclust = doubletFinder_v3(PH.filtNormScaledPCAclust, PCs = 1:20, pN = 0.25, pK = pK.PH, nExp = nExp_poi.PH, reuse.pANN = FALSE,sct = FALSE)
PH.filtNormScaledPCAclust = doubletFinder_v3(PH.filtNormScaledPCAclust, PCs = 1:20, pN = 0.25, pK = pK.PH, nExp = nExp_poi.adj.PH, reuse.pANN = paste("pANN_",nExp_poi.PH,sep=""), sct = FALSE)

SH.filtNormScaledPCAclust@meta.data[,"DF_Class"] = get('SH.filtNormScaledPCAclust')@meta.data[[paste("DF.classifications_",nExp_poi.SH,sep="")]]
#this next line doesn't make sense - telling it: which values in DF_Class are 'Doublets' and of those which are singlets, classify as Doublet_LOW
SH.filtNormScaledPCAclust@meta.data$DF_Class[which(SH.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet" & get('SH.filtNormScaledPCAclust')@meta.data[[paste("DF.classifications_0.25_0.09_",nExp_poi.SH,sep="")]] == "Singlet")] = "Doublet_LOW"
SH.filtNormScaledPCAclust@meta.data$DF_Class[which(SH.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet")] = "Doublet_HIGH"
SH.Doublets=as.character(colnames(SH.filtNormScaledPCAclust)[which(SH.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet_HIGH")])

PH.filtNormScaledPCAclust@meta.data[,"DF_Class"] = get('PH.filtNormScaledPCAclust')@meta.data[[paste("DF.classifications_",nExp_poi.PH,sep="")]]
PH.filtNormScaledPCAclust@meta.data$DF_Class[which(PH.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet" & get('PH.filtNormScaledPCAclust')@meta.data[[paste("DF.classifications_0.25_0.09_",nExp_poi.PH,sep="")]] == "Singlet")] = "Doublet_LOW"
PH.filtNormScaledPCAclust@meta.data$DF_Class[which(PH.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet")] = "Doublet_HIGH"
PH.Doublets=as.character(colnames(PH.filtNormScaledPCAclust)[which(PH.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet_HIGH")])

#plot to show doublets
pdf("SH_UMAP_doublets.pdf")
plot(DimPlot(SH.filtNormScaledPCAclust, reduction = "umap", group.by = colnames(as.matrix(SH.filtNormScaledPCAclust@meta.data[1,]))[9], label = TRUE))
dev.off()
pdf("PH_UMAP_doublets.pdf")
plot(DimPlot(PH.filtNormScaledPCAclust, reduction = "umap", group.by = colnames(as.matrix(PH.filtNormScaledPCAclust@meta.data[1,]))[9], label = TRUE))
dev.off()

#drop doublets
SH.filtNormScaledPCAclustDF=SH.filtNormScaledPCAclust[, !(colnames(SH.filtNormScaledPCAclust) %in% SH.Doublets), drop = FALSE]
PH.filtNormScaledPCAclustDF=PH.filtNormScaledPCAclust[, !(colnames(PH.filtNormScaledPCAclust) %in% PH.Doublets), drop = FALSE]

saveRDS(SH.filtNormScaledPCAclustDF, "SH_filtNormScaledPCAclustDF.rds")
saveRDS(SH.filtNormScaledPCAclustDF, "PH_filtNormScaledPCAclustDF.rds")

#find markers for each cluster
SH.dataMarkers= FindAllMarkers(SH.filtNormScaledPCAclustDF, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SH.dataMarkers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
SH.dataMarkers.LIST=unstack(SH.dataMarkers, SH.dataMarkers$gene ~ SH.dataMarkers$cluster)
saveRDS(SH.dataMarkers.LIST,"SH_markers.rds")

PH.dataMarkers= FindAllMarkers(PH.filtNormScaledPCAclustDF, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PH.dataMarkers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
PH.dataMarkers.LIST=unstack(PH.dataMarkers, PH.dataMarkers$gene ~ PH.dataMarkers$cluster)
saveRDS(PH.dataMarkers.LIST,"PH_markers.rds")

#load zeizal markers
zeizal.markers=readRDS("/mnt/data/Maya/SHvsPH/Analysis/SHvsPH/zeizal_markers.rds")

#functions for fishers exact test
SH.FisherTest.zeizal=function(x)
{
  TMP=matrix(ncol=9,nrow=1)
  Overlap.Genes=intersect(SH.Filtered.Genes,x)
  for(j in 1:9)
  {
    TMP.Genes=intersect(SH.Filtered.Genes,zeizal.markers[[j]])
    TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,Overlap.Genes))
    TMP.MAT[1,2]=length(setdiff(TMP.Genes,Overlap.Genes))
    TMP.MAT[2,1]=length(setdiff(Overlap.Genes,TMP.Genes))
    TMP.MAT[2,2]=length(SH.Filtered.Genes)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    TMP[1,j]=fisher.test(TMP.MAT,alternative="greater")$p.value
  }
  TMP
}

SH.CellTypeTest.zeizal=function(x)
{
  CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05 )]})
  CellType.Names=c(rep("Unclassified",length(SH.FisherTest.res.zeizal)))
  CellType.Names[which(CellType=="Mural")]="Mural"
  CellType.Names[which(CellType=="Endothelial")]="Endothelial"
  CellType.Names[which(CellType=="Ependymal")]="Ependymal"
  CellType.Names[which(CellType=="Microglia")]="Microglia"
  CellType.Names[which(CellType=="Oligodendrocyte")]="Oligodendrocyte"
  CellType.Names[which(CellType=="Interneuron")]="Interneurone"
  CellType.Names[which(CellType=="Pyramidal")]="CA1.Pyramidal"
  CellType.Names[which(CellType=="S1.Pyramidal")]="S1.Pyramidal"
  CellType.Names[which(CellType=="Astrocyte")]="Astrocyte"
  CellType.Names
}

PH.FisherTest.zeizal=function(x)
{
  TMP=matrix(ncol=9,nrow=1)
  Overlap.Genes=intersect(PH.Filtered.Genes,x)
  for(j in 1:9)
  {
    TMP.Genes=intersect(PH.Filtered.Genes,zeizal.markers[[j]])
    TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,Overlap.Genes))
    TMP.MAT[1,2]=length(setdiff(TMP.Genes,Overlap.Genes))
    TMP.MAT[2,1]=length(setdiff(Overlap.Genes,TMP.Genes))
    TMP.MAT[2,2]=length(PH.Filtered.Genes)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    TMP[1,j]=fisher.test(TMP.MAT,alternative="greater")$p.value
  }
  TMP
}

PH.CellTypeTest.zeizal=function(x)
{
  CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05 )]})
  CellType.Names=c(rep("Unclassified",length(PH.FisherTest.res.zeizal)))
  CellType.Names[which(CellType=="Mural")]="Mural"
  CellType.Names[which(CellType=="Endothelial")]="Endothelial"
  CellType.Names[which(CellType=="Ependymal")]="Ependymal"
  CellType.Names[which(CellType=="Microglia")]="Microglia"
  CellType.Names[which(CellType=="Oligodendrocyte")]="Oligodendrocyte"
  CellType.Names[which(CellType=="Interneuron")]="Interneurone"
  CellType.Names[which(CellType=="Pyramidal")]="CA1.Pyramidal"
  CellType.Names[which(CellType=="S1.Pyramidal")]="S1.Pyramidal"
  CellType.Names[which(CellType=="Astrocyte")]="Astrocyte"
  CellType.Names
}


#label clusters
SH.data.zeizal=SH.filtNormScaledPCAclustDF
SH.Filtered.Genes=unique(unlist(SH.dataMarkers.LIST))
SH.FisherTest.res.zeizal=lapply(SH.dataMarkers.LIST,SH.FisherTest.zeizal)
SH.TMP.zeizal = matrix(unlist(SH.FisherTest.res.zeizal), ncol = 9, byrow = TRUE)
colnames(SH.TMP.zeizal)=names(zeizal.markers)

SH.TMP.LABELS.zeizal=SH.CellTypeTest.zeizal(SH.TMP.zeizal)
names(SH.TMP.LABELS.zeizal)=names(SH.FisherTest.res.zeizal)
SH.zeizal.cluster.combined.ids = SH.TMP.LABELS.zeizal
names(SH.zeizal.cluster.combined.ids) = levels(SH.data.zeizal)
SH.data.zeizal = RenameIdents(SH.data.zeizal,SH.zeizal.cluster.combined.ids)
levels(SH.data.zeizal)

PH.data.zeizal=PH.filtNormScaledPCAclustDF
PH.Filtered.Genes=unique(unlist(PH.dataMarkers.LIST))
PH.FisherTest.res.zeizal=lapply(PH.dataMarkers.LIST,PH.FisherTest.zeizal)
PH.TMP.zeizal = matrix(unlist(PH.FisherTest.res.zeizal), ncol = 9, byrow = TRUE)
colnames(PH.TMP.zeizal)=names(zeizal.markers)

PH.TMP.LABELS.zeizal=PH.CellTypeTest.zeizal(PH.TMP.zeizal)
names(PH.TMP.LABELS.zeizal)=names(PH.FisherTest.res.zeizal)
PH.zeizal.cluster.combined.ids = PH.TMP.LABELS.zeizal
names(PH.zeizal.cluster.combined.ids) = levels(PH.data.zeizal)
PH.data.zeizal = RenameIdents(PH.data.zeizal,PH.zeizal.cluster.combined.ids)
levels(PH.data.zeizal)

pdf("SH_UMAP_labelled.pdf")
DimPlot(SH.data.zeizal, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()
table(Idents(SH.data.zeizal))
levels(SH.data.zeizal)

pdf("PH_UMAP_labelled.pdf")
DimPlot(PH.data.zeizal, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()
table(Idents(PH.data.zeizal))
levels(PH.data.zeizal)

SH.data.zeizal[["ClusterIdent"]] = Idents(object = SH.data.zeizal)
SH.ClusterIdent=as.character(SH.data.zeizal@meta.data$ClusterIdent)
names(SH.ClusterIdent)=rownames(SH.data.zeizal@meta.data)

PH.data.zeizal[["ClusterIdent"]] = Idents(object = PH.data.zeizal)
PH.ClusterIdent=as.character(PH.data.zeizal@meta.data$ClusterIdent)
names(PH.ClusterIdent)=rownames(PH.data.zeizal@meta.data)

saveRDS(as.factor(SH.ClusterIdent),"SH_CellID_Type.")
saveRDS(as.factor(PH.ClusterIdent),"PH_CellID_Type.")

#save matrix for conos
SH.data.finalMatrix=GetAssayData(object = SH.data.zeizal, assay = "RNA", slot = "data")
SH.data.finalMatrix.nonNorm=SH.MTfilt.data[,colnames(SH.data.finalMatrix)]
dim(SH.data.finalMatrix)
dim(SH.data.finalMatrix.nonNorm)
saveRDS(SH.data.finalMatrix.nonNorm,"SH_NonNormalised_filteredcounts_dgCMatrix.rds")

PH.data.finalMatrix=GetAssayData(object = PH.data.zeizal, assay = "RNA", slot = "data")
PH.data.finalMatrix.nonNorm=PH.MTfilt.data[,colnames(PH.data.finalMatrix)]
dim(PH.data.finalMatrix)
dim(PH.data.finalMatrix.nonNorm)
saveRDS(PH.data.finalMatrix.nonNorm,"PH_NonNormalised_filteredcounts_dgCMatrix.rds")

#read matrices for conos
tmp = list.files(pattern="*_NonNormalised_filteredcounts_dgCMatrix.rds")
Allfiles = lapply(tmp,readRDS)
filenames=list.files(pattern="*_NonNormalised_filteredcounts_dgCMatrix.rds", full.names=TRUE)
names(Allfiles) = gsub(".*/(.*)\\..*", "\\1", filenames)
names(Allfiles) = gsub("_NonNormalised_filteredcounts_dgCMatrix","", names(Allfiles))
Allfiles$PH@Dimnames[[2]]=gsub("SAMPLE", "PH", Allfiles$PH@Dimnames[[2]])
Allfiles$SH@Dimnames[[2]]=gsub("SAMPLE", "SH", Allfiles$SH@Dimnames[[2]])

#construct conos object
set.seed(1234)
Allfiles.Presprocessed = lapply(Allfiles, basicP2proc, n.cores=10, min.cells.per.gene=1, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
con=Conos$new(Allfiles.Presprocessed ,n.cores=4)
str(con$samples,1)

#conos clustering 
pdf("SH_PH_conos_clusters.pdf")
con$plotPanel(clustering="multilevel", use.local.clusters=T, title.size=6)
dev.off()

con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

#to clear cached PCAS
#con$pairs$PCA <- NULL
pdf("conos_PCS.pdf")
plotComponentVariance(con, space='PCA')
dev.off()

con$findCommunities(method=leiden.community, resolution=1)

pdf("SH_PH_same_conos_clusters.pdf")
con$plotPanel(font.size=4) 
dev.off()

pdf("combined_conos_clusters.pdf")
con$plotGraph(alpha=0.1)
dev.off()

pdf("combined_bysample_conos_clusters.pdf")
con$plotGraph(color.by='sample', mark.groups=F, alpha=0.1, show.legend=T)
dev.off()

CellAnnotate.SH =readRDS("SH_CellID_Type.rds")
names(CellAnnotate.SH)=gsub("SAMPLE", "SH", names(CellAnnotate.SH))

#CellAnnotate.SH= setNames(CellAnnotate.SH[,1], CellAnnotate.SH[,2])
#con$plotPanel(groups = CellAnnotate.SH)

LabelProb.info= con$propagateLabels(labels = CellAnnotate.SH, verbose=T)

pdf("conos_labelled.pdf")
con$plotPanel(colors=LabelProb.info$uncertainty, show.legend=T, legend.title="Uncertainty", legend.pos=c(1, 0))
dev.off()

pdf("conos_labelled_groups.pdf")
con$plotPanel(groups=LabelProb.info$labels, show.legend=F)
dev.off()


NewAnnotation.SH = LabelProb.info$labels
table(NewAnnotation.SH)

table(con$clusters[[1]][1])
sum(table(con$clusters[[1]][1]))
sum(table(NewAnnotation.SH))


write.table(NewAnnotation.SH,"SH_PH_annotations.csv")
saveRDS(con,"SH_PH_final_conos_object.rds")
saveRDS(PH.data.finalMatrix,"PH_final_matrix.rds")
saveRDS(PH.data.finalMatrix.nonNorm,"PH_final_matrix_nonNorm.rds")
saveRDS(SH.data.finalMatrix,"SH_final_matrix.rds")
saveRDS(SH.data.finalMatrix.nonNorm,"SH_final_matrix_nonNorm.rds")
saveRDS(SH.ClusterIdent,"SH_cluster_ident.rds")
saveRDS(SH.data.zeizal,"SH_data_zeizal.rds")
saveRDS(PH.data.zeizal,"PH_data_zeizal.rds")


NewAnnotation.SH=read.csv("SH_PH_annotations.csv")
con=readRDS("SH_PH_final_conos_object.rds")
PH.data.finalMatrix=readRDS("PH_final_matrix.rds")
PH.data.finalMatrix.nonNorm=readRDS("PH_final_matrix_nonNorm.rds")
SH.data.finalMatrix=readRDS("SH_final_matrix.rds")
SH.data.finalMatrix.nonNorm=readRDS("SH_final_matrix_nonNorm.rds")
SH.ClusterIdent=readRDS("SH_cluster_ident.rds")
SH.data.zeizal=readRDS("SH_data_zeizal.rds")
PH.data.zeizal=readRDS("PH_data_zeizal.rds")

GroupS= con$getDatasetPerCell() %>% substr(1,11) %>% setNames(names(con$getDatasetPerCell()))

MetaData=Reduce(cbind,list(as.data.frame(GroupS)))
write.csv(MetaData,"MetaData.csv")

#mast
#read in files
tmp = list.files(pattern="*_NonNormalised_filteredcounts_dgCMatrix.rds")
Allfiles = lapply(tmp,readRDS)
filenames=list.files(pattern="*_NonNormalised_filteredcounts_dgCMatrix.rds", full.names=TRUE)
names(Allfiles) = gsub(".*/", "", filenames)
names(Allfiles) = gsub("_NonNormalised_filteredcounts_dgCMatrix.rds", "", names(Allfiles))

Allfiles$PH@Dimnames[[2]]=gsub("SAMPLE", "PH", Allfiles$PH@Dimnames[[2]])
Allfiles$SH@Dimnames[[2]]=gsub("SAMPLE", "SH", Allfiles$SH@Dimnames[[2]])



freq_expressed <- 0.2
FcThreshold <- log2(1.5)

MetaData=read.csv("MetaData.csv",sep=",",row.names=1)
rownames(MetaData)=gsub("-1",".1",rownames(MetaData))
MetaData$Cell_IDs=rownames(MetaData)
Table=Reduce(merge, lapply(Allfiles, function(x) data.frame(x, rn = row.names(x))))
rownames(Table)=Table$rn
Table$rn=NULL
Table.MX.cpm=cpm(Table)
saveRDS(Table.MX.cpm,"Table.MX.cpm.rds")

NewAnnotation.SH2=as.data.frame(NewAnnotation.SH, row.names=names(NewAnnotation.SH))
#NewAnnotation.SH2=str_split_fixed(NewAnnotation.SH$x, " ", 2)
#NewAnnotation.SH.2=as.data.frame(NewAnnotation.SH2)
#head(colnames(Table.MX.cpm))
#NewAnnotation.SH.2$V1=gsub("-1",".1",NewAnnotation.SH.2$V1)
rownames(NewAnnotation.SH2)=gsub("-1",".1", rownames(NewAnnotation.SH2))
#rownames(NewAnnotation.SH.2)=NewAnnotation.SH.2$V1


#ASTROCYTE
NewAnnotation.SH.astro=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="Astrocyte")]
head(NewAnnotation.SH.astro)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.astro,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_astrocyte.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_astrocyte.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_astrocyte.csv")


#CA1.PYRAMIDAL
NewAnnotation.SH.CA1=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="CA1.Pyramidal")]
head(NewAnnotation.SH.CA1)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.CA1,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_CA1.Pyramidal.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_CA1.Pyramidal.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_CA1.Pyramidal.csv")

#MICROGLIA
NewAnnotation.SH.MIC=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="Microglia")]
head(NewAnnotation.SH.MIC)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.MIC,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_Microglia.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_Microglia.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_Microglia.csv")


#OLIGODENDROCYTE
NewAnnotation.SH.OLI=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="Oligodendrocyte")]
head(NewAnnotation.SH.OLI)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.OLI,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca.oligo = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca.oligo)>0)
colData(sca.oligo)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca.oligo)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca.oligo)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_Oligodendrocyte.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_Oligodendrocyte.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_Oligodendrocyte.csv")

length(NewAnnotation.SH.CA1)
#S1.PYRAMIDAL
NewAnnotation.SH.S1=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="S1.Pyramidal")]
head(NewAnnotation.SH.S1)
length(NewAnnotation.SH.S1)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.S1,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca.s1 = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_S1.Pyramidal.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_S1.Pyramidal.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_S1.Pyramidal.csv")

#INTERNEURON
NewAnnotation.SH.INT=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="Interneurone")]
head(NewAnnotation.SH.INT)
length(NewAnnotation.SH.INT)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.INT,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_Interneuron.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_Interneuron.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_Interneuron.csv")


#ENDOTHELIAL
NewAnnotation.SH.ENDO=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="Endothelial")]
head(NewAnnotation.SH.ENDO)
length(NewAnnotation.SH.ENDO)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.ENDO,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_Endothelial.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_Endothelial.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_Endothelial.csv")



#UNCLASSIFIED
NewAnnotation.SH.UNC=rownames(NewAnnotation.SH2)[which(NewAnnotation.SH2$NewAnnotation.SH=="Unclassified")]
head(NewAnnotation.SH.UNC)
length(NewAnnotation.SH.UNC)
Table.MX.cpm3=Table.MX.cpm[,intersect(NewAnnotation.SH.UNC,colnames(Table.MX.cpm))]

Ex.MetaData=MetaData[colnames(Table.MX.cpm3),] # MetaDATA Control or case Gender, 
fdata=data.frame(rownames(Table.MX.cpm3))
rownames(fdata)=rownames(Table.MX.cpm3)
saveRDS(fdata,"fdata_SH_PH.rds")
saveRDS(Table.MX.cpm3,"Table_MX_cpm_Ex_.rds")
saveRDS(Ex.MetaData,"Ex_MetaData.rds")

sca = FromMatrix(log2(as.matrix(Table.MX.cpm3)+1),Ex.MetaData,fdata)

Detected_Genes =colSums(assay(sca)>0)
colData(sca)$DetRate = scale(Detected_Genes)
cond=factor(colData(sca)$GroupS)
head(cond)
cond= relevel(cond,"SH")
colData(sca)$condition=cond

zlmCond=zlm(~condition + DetRate,sca)
summaryCond = summary(zlmCond, doLRT='conditionPH') 

saveRDS(summaryCond,"SummaryCond_SH_PH.rds")


summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='conditionPH' & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast=='conditionPH' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig = merge(fcHurdle[fdr <.05 & abs(coef)>FcThreshold], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
fcHurdleSig

fcHurdleSig2=fcHurdleSig
rownames(fcHurdleSig2)=fcHurdleSig2$primerid  
fcHurdleSig2$primerid=NULL
saveRDS(fcHurdle,"FCHURDLE.SH_PH_Unclassified.rds")
write.csv(fcHurdle,"FCHURDLE_SH_PH_Unclassified.csv")
write.csv(fcHurdleSig2,"Sig_Genes_SH_PH_Unclassified.csv")
