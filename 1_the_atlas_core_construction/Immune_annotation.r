library(qs)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(patchwork)
library("scCustomize")
library(sceasy)
library(reticulate)
library(biomaRt)
library("DoubletFinder")

health_merge_4qc<-qread("/bmbl_data/qiguo/SCI/atlas/integration/processeddata/health_merge_qc_0227_allannotated_4qc.qs")
Idents(health_merge_4qc)<-health_merge_4qc$harmonized_majorct2
immunecells<-c("Microglia","Myeloid cells","Lymphocytes","Unknown")
health_merge_immune_4qc<-subset(health_merge_4qc,idents = immunecells)
embedding<-read.csv("/bmbl_data/qiguo/SCI/atlas/integration/output/scgen_healthy_immune_0227/healthy_corrected_immune_0227_latent.csv")
rownames(embedding)<-embedding$X
embedding<-embedding[,-1]
dim(embedding)
dim(health_merge_immune_4qc)
identical(colnames(health_merge_immune_4qc),rownames(embedding))#true, good
#
colnames(embedding)<-paste0("integratedscgen_",1:100)
embedding<-as.matrix(embedding)
health_merge_immune_4qc[["integrated.scgen"]]<-CreateDimReducObject(embeddings = embedding, key = 'scgenre_', assay = 'RNA')
health_merge_immune_4qc <- FindNeighbors(health_merge_immune_4qc,  reduction = "integrated.scgen", dims = 1:30)
health_merge_immune_4qc <- RunUMAP(health_merge_immune_4qc,reduction = "integrated.scgen", dims = 1:40)
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct2")

#need normalization for visualization
health_merge_immune_4qc[["RNA"]] <- split(health_merge_immune_4qc[["RNA"]], f = health_merge_immune_4qc$purestudy_firstname)
health_merge_immune_4qc <- NormalizeData(health_merge_immune_4qc)
health_merge_immune_4qc <- FindVariableFeatures(health_merge_immune_4qc)
health_merge_immune_4qc <- ScaleData(health_merge_immune_4qc)
health_merge_immune_4qc <- RunPCA(health_merge_immune_4qc)

#clustering
health_merge_immune_4qc <- FindClusters(health_merge_immune_4qc, resolution = 0.5, cluster.name = "cluster_scgen_r")
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "cluster_scgen_r",label = T)
DotPlot(health_merge_immune_4qc,features = c("Ccl5","Ms4a4b","Nkg7","Cd79a","Cd79b","P2ry12","Tmem119","Siglech","Mrc1","Gpnmb","Spp1","Fabp5",
                                       "H2-Aa","H2-Eb1","H2-Ab1","S100a9","S100a8","Retnlg","Atp1a2","Slc1a2","Atp1b2"))
#annotation

#B cells (7)
DotPlot(health_merge_immune_4qc,features = c("Cd79a","Cd79b","P2ry12"))
health_merge_immune_4qc$harmonized_majorct3<-as.character(health_merge_immune_4qc$cluster_scgen_r)
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r
health_merge_immune_4qc <- FindSubCluster(health_merge_immune_4qc, cluster = "7",resolution = 0.5, subcluster.name = "cluster_scgen_r_7",graph.name = "RNA_snn")
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r_7
DotPlot(health_merge_immune_4qc,features = c("Cd79a","Cd79b","Ccl5","Ms4a4b","Nkg7","P2ry12"))
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_7=="7_1")]<-"B cells"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_7=="7_3")]<-"B cells"

#T cells (4)
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_7=="7_2")]<-"T cells"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_7=="7_0")]<-"T cells"
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct3",label = T)

#Neutrophil
DotPlot(health_merge_immune_4qc,features = c("S100a9","S100a8","Retnlg","P2ry12"))
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="8")]<-"Neutrophils"

#microglia
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r
DotPlot(health_merge_immune_4qc,features = c("P2ry12","Tmem119","Ly6c1","Cldn5"))
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="6")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="4")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="3")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="0")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="1")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="11")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="10")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="12")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$harmonized_majorct3=="2")]<-"Microglia"

#endothelial
health_merge_immune_4qc <- FindSubCluster(health_merge_immune_4qc, cluster = "9",resolution = 0.5, subcluster.name = "cluster_scgen_r_9",graph.name = "RNA_snn")
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r_9
DotPlot(health_merge_immune_4qc,features = c("P2ry12","Tmem119","Ly6c1","Cldn5"))
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "cluster_scgen_r_9",label = T)
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_9=="9_2")]<-"Endothelial cells"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_9=="9_0")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_9=="9_1")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_9=="9_3")]<-"9_3"
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct3",label = T)


#BAM (3)
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$harmonized_majorct3
DotPlot(health_merge_immune_4qc,features = c("Lyve1","Mrc1", "Cd163","Ccr2","Gpnmb","Spp1","Fabp5","P2ry12","Tmem119"))
health_merge_immune_4qc <- FindSubCluster(health_merge_immune_4qc, cluster = "5",resolution = 1, subcluster.name = "cluster_scgen_r_5",graph.name = "RNA_snn")
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r_5
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "cluster_scgen_r_5",label = T)
DotPlot(health_merge_immune_4qc,features = c("Lyve1","Mrc1", "Cd163","Ccr2","Gpnmb","Spp1","Fabp5","P2ry12","Tmem119"))
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_10")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_6")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="9_3")]<-"Other PMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_9")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_8")]<-"Other PMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_7")]<-"BAMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_5")]<-"Other PMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_3")]<-"Other PMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_4")]<-"BAMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_2")]<-"5_2"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_1")]<-"5_1"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_5=="5_0")]<-"5_0"
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct3",label = T)
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct2",label = T)
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$harmonized_majorct3
DotPlot(health_merge_immune_4qc,features = c("Lyve1","Mrc1", "Cd163","Ccr2","Gpnmb","Spp1","Fabp5","P2ry12","Tmem119"))

#subcluster 5_0
health_merge_immune_4qc <- FindSubCluster(health_merge_immune_4qc, cluster = "5_0",resolution = 1, subcluster.name = "cluster_scgen_r_50",graph.name = "RNA_snn")
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r_50
DotPlot(health_merge_immune_4qc,features = c("Lyve1","Mrc1", "Cd163","Ccr2","Gpnmb","Spp1","Fabp5","P2ry12","Tmem119"))
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "cluster_scgen_r_50",label = T)
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_50=="5_0_2")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_50=="5_0_0")]<-"BAMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_50=="5_0_1")]<-"BAMs"

#subcluster 5_1
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$harmonized_majorct3
health_merge_immune_4qc <- FindSubCluster(health_merge_immune_4qc, cluster = "5_1",resolution = 1, subcluster.name = "cluster_scgen_r_51",graph.name = "RNA_snn")
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r_51
DotPlot(health_merge_immune_4qc,features = c("Lyve1","Mrc1", "Cd163","Ccr2","Gpnmb","Spp1","Fabp5","P2ry12","Tmem119"))
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct3",label = T)
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_51=="5_1_3")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_51=="5_1_2")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_51=="5_1_1")]<-"Microglia"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_51=="5_1_0")]<-"BAMs"

#subcluster 5_2
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$harmonized_majorct3
health_merge_immune_4qc <- FindSubCluster(health_merge_immune_4qc, cluster = "5_2",resolution = 1, subcluster.name = "cluster_scgen_r_52",graph.name = "RNA_snn")
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$cluster_scgen_r_52
DotPlot(health_merge_immune_4qc,features = c("Lyve1","Mrc1", "Cd163","Ccr2","Gpnmb","Spp1","Fabp5","P2ry12","Tmem119"))
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct3",label = T)
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_52=="5_2_0")]<-"BAMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_52=="5_2_1")]<-"BAMs"
health_merge_immune_4qc$harmonized_majorct3[which(health_merge_immune_4qc$cluster_scgen_r_52=="5_2_2")]<-"Microglia"
table(health_merge_immune_4qc$harmonized_majorct3)
DimPlot(health_merge_immune_4qc, reduction = "umap",group.by = "harmonized_majorct3",label = T)
Idents(health_merge_immune_4qc)<-health_merge_immune_4qc$harmonized_majorct3
DotPlot(health_merge_immune_4qc,features = c("Lyve1","Mrc1", "Cd163","Ccr2","Gpnmb","Spp1","Fabp5","P2ry12","Tmem119"))
qsave(health_merge_immune_4qc,"/bmbl_data/qiguo/SCI/atlas/integration/processeddata/health_immune_0303_allannotated_4qc.qs")

health_merge_immune_4qc<-qread("/bmbl_data/qiguo/SCI/atlas/integration/processeddata/health_immune_0303_allannotated_4qc.qs")

