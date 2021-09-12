# Lung, Oesophagus, Spleen data downloaded from https://www.tissuestabilitycellatlas.org/
## Jing Jiao Li biorxiv 2020
dat = readRDS("C:/Users/ishaa/Desktop/lung_ts.rds")
jiao2020 = read.delim("C:/Users/ishaa/Downloads/COVID_Genes_Jingjiao_Li_biorxiv2020 - Sheet1.tsv",header=T)
j.split = split(jiao2020,jiao2020$Localization)
dat@meta.data$jiaoSet = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% jiao2020$Gene.Names,]))
dat@meta.data[,"test"] = colSums(as.matrix(dat@assays$RNA@counts[1:5,]))
for(x in names(j.split))
  dat@meta.data[,x] = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% jiao2020$Gene.Names[jiao2020$Localization==x],]))
Organelles = c("Cytoskeleton","Cytosol","Endoplasmic Reticulum","Golgi Apparatus","Mitochondrion","Nucleus","jiaoSet")
par(mar=c(15,10,2,2));
for(x in Organelles){
  boxplot( dat@meta.data[,x]/dat@meta.data$n_counts ~dat@meta.data$Celltypes,outline=F,main=paste("Lung",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
}
summary(aov(dat@meta.data[,"jiaoSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))

##################
dat = readRDS("C:/Users/ishaa/Desktop/oesophagus_ts.rds")
jiao2020 = read.delim("C:/Users/ishaa/Downloads/COVID_Genes_Jingjiao_Li_biorxiv2020 - Sheet1.tsv",header=T)
j.split = split(jiao2020,jiao2020$Localization)dat@meta.data$jiaoSet = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% jiao2020$Gene.Names,]))
dat@meta.data[,"test"] = colSums(as.matrix(dat@assays$RNA@counts[1:5,]))
for(x in names(j.split))
  dat@meta.data[,x] = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% jiao2020$Gene.Names[jiao2020$Localization==x],]))
Organelles = c("Cytoskeleton","Cytosol","Endoplasmic Reticulum","Golgi Apparatus","Mitochondrion","Nucleus","jiaoSet")

par(mar=c(15,10,2,2));
for(x in Organelles)
{
  boxplot( dat@meta.data[,x]/dat@meta.data$n_counts ~dat@meta.data$Celltypes,outline=F,main=paste("Oesophagus",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
}
summary(aov(dat@meta.data[,"jiaoSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))
###################
dat = readRDS("C:/Users/ishaa/Desktop/spleen_ts.rds")
jiao2020 = read.delim("C:/Users/ishaa/Downloads/COVID_Genes_Jingjiao_Li_biorxiv2020 - Sheet1.tsv",header=T)
j.split = split(jiao2020,jiao2020$Localization)dat@meta.data$jiaoSet = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% jiao2020$Gene.Names,]))
dat@meta.data[,"test"] = colSums(as.matrix(dat@assays$RNA@counts[1:5,]))
for(x in names(j.split))
  dat@meta.data[,x] = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% jiao2020$Gene.Names[jiao2020$Localization==x],]))Organelles = c("Cytoskeleton","Cytosol","Endoplasmic Reticulum","Golgi Apparatus","Mitochondrion","Nucleus","jiaoSet")
par(mar=c(15,10,2,2));
for(x in Organelles){
  boxplot( dat@meta.data[,x]/dat@meta.data$n_counts ~dat@meta.data$Celltypes,outline=F,main=paste("Spleen",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
}
summary(aov(dat@meta.data[,"jiaoSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))

## Krogan biorxiv 2020
dat = readRDS("C:/Users/ishaa/Desktop/lung_ts.rds")
dat = FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
krogan = read.delim("C:/Users/ishaa/Downloads/COVID_Genes_Krogan_biorxiv2020 - Sheet1.csv",sep=",",header=T,stringsAsFactors = F)
krogan.split = split(krogan, krogan$Viral.bait)
dat@meta.data$kroganSet = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene,]))
dat@meta.data$varGeneSum = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% dat@assays$RNA@var.features,]))
for(x in names(krogan.split))
  dat@meta.data[,x] = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene[krogan$Viral.bait==x],]))par(mar=c(15,10,2,2));
for(x in c(unique(krogan$Viral.bait),"kroganSet"))
{
  boxplot( dat@meta.data[,x]/dat@meta.data$n_counts ~dat@meta.data$Celltypes,
           outline=F,main=paste("Spleen",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)  
  #  boxplot( dat@meta.data[,x]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes,
  #            outline=F,main=paste("Lung",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
}
summary(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))
TukeyHSD(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))
summary(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes))
TukeyHSD(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes))
##################
dat = readRDS("C:/Users/ishaa/Desktop/oesophagus_ts.rds")
dat = FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
krogan = read.delim("C:/Users/ishaa/Downloads/COVID_Genes_Krogan_biorxiv2020 - Sheet1.csv",sep=",",header=T,stringsAsFactors = F)
krogan.split = split(krogan, krogan$Viral.bait)
dat@meta.data$kroganSet = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene,]))
dat@meta.data$varGeneSum = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% dat@assays$RNA@var.features,]))
for(x in names(krogan.split))
  dat@meta.data[,x] = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene[krogan$Viral.bait==x],]))par(mar=c(15,10,2,2));
for(x in c(unique(krogan$Viral.bait),"kroganSet"))
{
  boxplot( dat@meta.data[,x]/dat@meta.data$n_counts ~dat@meta.data$Celltypes,
           outline=F,main=paste("Spleen",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
  #  boxplot( dat@meta.data[,x]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes,
  #           outline=F,main=paste("Oesophagus",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
}

summary(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))
TukeyHSD(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))
summary(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes))
TukeyHSD(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes))

###################
dat = readRDS("C:/Users/ishaa/Desktop/spleen_ts.rds")
dat = FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
krogan = read.delim("C:/Users/ishaa/Downloads/COVID_Genes_Krogan_biorxiv2020 - Sheet1.csv",sep=",",header=T,stringsAsFactors = F)
krogan.split = split(krogan, krogan$Viral.bait)dat@meta.data$kroganSet = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene,]))
dat@meta.data$varGeneSum = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% dat@assays$RNA@var.features,]))
for(x in names(krogan.split))
  dat@meta.data[,x] = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene[krogan$Viral.bait==x],]))par(mar=c(15,10,2,2));
for(x in c(unique(krogan$Viral.bait),"kroganSet"))
{
  boxplot( dat@meta.data[,x]/dat@meta.data$n_counts ~dat@meta.data$Celltypes,
           outline=F,main=paste("Spleen",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
  #  boxplot( dat@meta.data[,x]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes,
  #           outline=F,main=paste("Spleen",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)}summary(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))
  TukeyHSD(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$n_counts ~dat@meta.data$Celltypes))
  summary(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes))
  TukeyHSD(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$varGeneSum ~dat@meta.data$Celltypes))
  ## PBMCs
  dat = readRDS("C:/Users/ishaa/Downloads/pbmc3k_final.rds")
  dat = FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
  krogan = read.delim("C:/Users/ishaa/Downloads/COVID_Genes_Krogan_biorxiv2020 - Sheet1.csv",sep=",",header=T,stringsAsFactors = F)
  krogan.split = split(krogan, krogan$Viral.bait)
  dat@meta.data$kroganSet = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene,]))
  for(x in names(krogan.split))
    dat@meta.data[,x] = colSums(as.matrix(dat@assays$RNA@counts[rownames(dat@assays$RNA@counts) %in% krogan$Host.gene[krogan$Viral.bait==x],]))
  dat@meta.data$Celltypes = Idents(dat)
  par(mar=c(15,10,2,2));
  for(x in c(unique(krogan$Viral.bait),"kroganSet")){
    boxplot( (dat@meta.data[,x]/dat@meta.data$nCount_RNA) ~dat@meta.data$Celltypes,
             outline=F,main=paste("PBMC-2.7k_",x),las=2 ,horizontal = T,xlab="",ylab="",cex.axis=0.75)
  }
  summary(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$nCount_RNA ~dat@meta.data$Celltypes))
  TukeyHSD(aov(dat@meta.data[,"kroganSet"]/dat@meta.data$nCount_RNA ~dat@meta.data$Celltypes))
}

## dECOVOLUTION using Bulk data
library(MuSiC)
library(Biobase)
library(xbioc)

dat = readRDS("C:/Users/ishaa/Desktop/lung_ts.rds")dat.m = as.matrix(dat@assays$RNA@scale.data)
pheno.m = dat@meta.datahlca.eset = ExpressionSet(assayData = as.matrix(dat.m), phenoData = new("AnnotatedDataFrame",data=pheno.m))bulk.data = read.delim("C:/Users/ishaa/Downloads/GSE147507_RawReadCounts_Human.tsv/GSE147507_RawReadCounts_Human.tsv")
bulk.mat = as.matrix(bulk.data[,grep("Lung",colnames(bulk.data))])
rownames(bulk.mat) = bulk.data[,1]bulk.mat.matched = (bulk.mat[which(rownames(bulk.mat) %in% rownames(dat.m)),])
bulk.eset = ExpressionSet(assayData = bulk.mat.matched )
Est.prop.bulk.annot = music_prop(bulk.eset = bulk.eset, 
                                 sc.eset = hlca.eset,
                                 clusters = 'Celltypes',
                                 samples="sample",
                                 select.ct = names(table(dat@meta.data$Celltypes)),
                                 verbose = T)
pheatmap::pheatmap(Est.prop.bulk.annot$Est.prop.weighted,cluster_rows = F,cluster_cols = F,display_numbers = T)
pheatmap::pheatmap(Est.prop.bulk.annot$Est.prop.allgene,cluster_rows = F,cluster_cols = F,display_numbers = T)

bulk.mat = as.matrix(bulk.data[,-1])
rownames(bulk.mat) = bulk.data[,1]
bulk.mat.matched = (bulk.mat[which(rownames(bulk.mat) %in% rownames(dat.m)),])
bulk.eset = ExpressionSet(assayData = bulk.mat.matched )
Est.prop.bulk.annot = music_prop(bulk.eset = bulk.eset,
                                 sc.eset = hlca.eset,
                                 clusters = 'Celltypes',
                                 samples="sample",
                                 select.ct = names(table(dat@meta.data$Celltypes)),
                                 verbose = T)
pheatmap::pheatmap(Est.prop.bulk.annot$Est.prop.weighted,cluster_rows = F,cluster_cols = F,display_numbers = T)
pheatmap::pheatmap(Est.prop.bulk.annot$Est.prop.allgene,cluster_rows = F,cluster_cols = F,display_numbers = T)