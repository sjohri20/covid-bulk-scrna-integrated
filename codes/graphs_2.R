
dat@meta.data[,"covid_up"] = colSums(dat[rownames(dat) %in% covid.up,])
dat@meta.data[,"covid_down"] = colSums(dat[rownames(dat) %in% covid.down,])
par(mar=c(5,7,2,2), mfrow=c(1,1))

boxplot(dat@meta.data$covid_up/dat@meta.data$n_counts ~dat@meta.data$Celltypes, horizontal=T, las=2, xlab="", ylab="", cex.axis=0.5, col="grey", )


dat.2 = readRDS("spleen_ts.rds")
dat.2@meta.data[,"covid_up"] = colSums(dat.2[rownames(dat.2) %in% covid.up,])
dat.2@meta.data[,"covid_down"] = colSums(dat.2[rownames(dat.2) %in% covid.down,])


par(mar=c(5,7,2,2), mfrow=c(1,1))

boxplot(dat@meta.data$covid_up/dat@meta.data$n_counts ~dat@meta.data$Celltypes, horizontal=T, las=2, xlab="", ylab="", cex.axis=0.5, col="grey", main="Spleen Covid Up")
boxplot(dat@meta.data$covid_down/dat@meta.data$n_counts ~dat@meta.data$Celltypes, horizontal=T, las=2, xlab="", ylab="", cex.axis=0.5, col="grey", main="Spleen Covid Down")

dat.3 = readRDS("oesophagus_ts.rds")
dat.3@meta.data[,"covid_up"] = colSums(dat.3[rownames(dat.3) %in% covid.up,])
dat.3@meta.data[,"covid_down"] = colSums(dat.3[rownames(dat.3) %in% covid.down,])


dat.4 = readRDS("pbmc3k_final.rds")
dat.4@meta.data[,"covid_up"] = colSums(dat.4[rownames(dat.4) %in% covid.up,])
dat.4@meta.data[,"covid_down"] = colSums(dat.4[rownames(dat.4) %in% covid.down,])


jiao2020 = read.delim("COVID_Genes_Jingjiao_Li_biorxiv2020 - Sheet1.tsv",header=T)
krogan2020 = read.delim("COVID_Genes_Krogan_biorxiv2020 - Sheet1.tsv",header=T)

j.split = split(jiao2020,jiao2020$Localization)
k.split = split(krogan2020,krogan2020$Annotation)



dat@meta.data[,"jiao2020"] = colSums(dat[rownames(dat) %in% jiao2020$Gene.Names,])
dat@meta.data[,"krogan2020"] = colSums(dat[rownames(dat) %in% krogan2020$Host.gene,])
dat@meta.data[,"both"] = dat@meta.data$jiao2020 + dat@meta.data$krogan2020

dat.2@meta.data[,"jiao2020"] = colSums(dat.2[rownames(dat.2) %in% jiao2020$Gene.Names,])
dat.2@meta.data[,"krogan2020"] = colSums(dat.2[rownames(dat.2) %in% krogan2020$Host.gene,])
dat.2@meta.data[,"both"] = dat.2@meta.data$jiao2020 + dat.2@meta.data$krogan2020

dat.3@meta.data[,"jiao2020"] = colSums(dat.3[rownames(dat.3) %in% jiao2020$Gene.Names,])
dat.3@meta.data[,"krogan2020"] = colSums(dat.3[rownames(dat.3) %in% krogan2020$Host.gene,])
dat.3@meta.data[,"both"] = dat.3@meta.data$jiao2020 + dat.3@meta.data$krogan2020

dat.4@meta.data[,"jiao2020"] = colSums(dat.4[rownames(dat.4) %in% jiao2020$Gene.Names,])
dat.4@meta.data[,"krogan2020"] = colSums(dat.4[rownames(dat.4) %in% krogan2020$Host.gene,])
dat.4@meta.data[,"both"] = dat.4@meta.data$jiao2020 + dat.4@meta.data$krogan2020



boxplot(dat@meta.data$covid_up/dat@meta.data$n_counts ~dat@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="Lung Covid Down")
boxplot(dat@meta.data$covid_down/dat@meta.data$n_counts ~dat@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="Lung Covid Up")

boxplot(dat.2@meta.data$covid_up/dat.2@meta.data$n_counts ~dat.2@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.5, col="grey", main="Spleen Covid Down")
boxplot(dat.2@meta.data$covid_down/dat.2@meta.data$n_counts ~dat.2@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.5, col="grey", main="Spleen Covid Up")

boxplot(dat.3@meta.data$covid_up/dat.3@meta.data$n_counts ~dat.3@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="Oesophagus Covid Down")
boxplot(dat.3@meta.data$covid_down/dat.3@meta.data$n_counts ~dat.3@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="Oesophagus Covid Up")

boxplot(dat.4@meta.data$covid_up/dat.4@meta.data$nCount_RNA ~dat.4@meta.data$seurat_clusters, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="PBMC Covid Down", names=c("Native CD4+ T", "Memory CD4+", "CD14+ Mono", "B", "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet"))
boxplot(dat.4@meta.data$covid_down/dat.4@meta.data$nCount_RNA ~dat.4@meta.data$seurat_clusters, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="PBMC Covid Up", names=c("Native CD4+ T", "Memory CD4+", "CD14+ Mono", "B", "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet"))

boxplot(dat.2@meta.data$both/dat.2@meta.data$n_counts ~dat.2@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.5, col="grey", main="Spleen Both")
boxplot(dat@meta.data$both/dat@meta.data$n_counts ~dat@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="Lung Both")
boxplot(dat.3@meta.data$both/dat.3@meta.data$n_counts ~dat.3@meta.data$Celltypes, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="Oesophagus Both")
boxplot(dat.4@meta.data$both/dat.4@meta.data$nCount_RNA ~dat.4@meta.data$seurat_clusters, horizontal=T, las=2, outline=F, xlab="", ylab="", cex.axis=0.6, col="grey", main="PBMC Both", names=c("Native CD4+ T", "Memory CD4+", "CD14+ Mono", "B", "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet"))

par(mfrow=c(2,2))
FeaturePlot(dat, features = "both")
FeaturePlot(dat.3, features = "both")
FeaturePlot(dat.4, features = "both")
FeaturePlot(dat.2, features = "both")


