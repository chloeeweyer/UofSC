## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------------
#installing necessary packages
#install.packages("Seurat")
#install.packages("tidyverse")
#install.packages("data.table")
#install.packages("LaF")


## ----------------------------------------------------------------------------------------------------------------------------------------
#loading libraries
library(Seurat) #"designed for QC, analysis, and exploration of single-cell RNA-seq data"
library(dplyr) #adding verbage for simplified data manipulation
library(patchwork) #for combining separate ggplots
library(readr) #for functions such as read.csv, read.delim
library("data.table") #for fread function, fastest way to read in large .csv files
library(knitr)


## ----------------------------------------------------------------------------------------------------------------------------------------
#reading in 25,000 columns (cells) of COPD csv file since raw file of ~95,000 cells is way too large. 
COPDhead<-fread("C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\W3_COPDdata.csv", header = FALSE, sep = ",", select = c(1:25000))


## ----------------------------------------------------------------------------------------------------------------------------------------
#now going to create a df we can use as row names
rnames<-COPDhead$V1 


## ----------------------------------------------------------------------------------------------------------------------------------------
#must set column v1 to row names to get gene names instead of numbers later
rownames(COPDhead) <- rnames


## ----------------------------------------------------------------------------------------------------------------------------------------
#now must set v1 to null 
COPDhead$V1 <- NULL


## ----------------------------------------------------------------------------------------------------------------------------------------
#creating seurat object. min.cells will include features that are in at least 3 cells. min features will include cells where over 200 features are detected.
COPD <- CreateSeuratObject(counts = COPDhead, project = "W3practical", min.cells = 3, min.features = 200)
#COPD


## ----------------------------------------------------------------------------------------------------------------------------------------
#next 25k cells (delete .csv dfs between these nexts 3 steps to save space)
COPDhead1<-fread("C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\W3_COPDdata.csv", header = FALSE, sep = ",", select = c(25001:50000))
rownames(COPDhead1) <- rnames
COPD1 <- CreateSeuratObject(counts = COPDhead1, project = "W3practical", min.cells = 3, min.features = 200)


## ----------------------------------------------------------------------------------------------------------------------------------------
remove(COPDhead1)
gc()
COPDhead2<-fread("C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\W3_COPDdata.csv", header = FALSE, sep = ",", select = c(50001:75000))
rownames(COPDhead2) <- rnames
COPD2 <- CreateSeuratObject(counts = COPDhead2, project = "W3practical", min.cells = 3, min.features = 200)


## ----------------------------------------------------------------------------------------------------------------------------------------
remove(COPDhead2)
gc()
COPDhead3<-fread("C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\W3_COPDdata.csv", header = FALSE, sep = ",", select = c(75001:94410))
rownames(COPDhead3) <- rnames
COPD3 <- CreateSeuratObject(counts = COPDhead3, project = "W3practical", min.cells = 3, min.features = 200)
remove(COPDhead3)
remove(COPDhead)
gc()


## ----------------------------------------------------------------------------------------------------------------------------------------
#merge 4 resulting seurat objects
COPDm <- merge(COPD, y = c(COPD1, COPD2, COPD3), add.cell.ids = c("1through25k", "25kthrough50k", "50kthrough75k" , "75kthrough95k"), project = "W3practical")
#remove old objects in terminal or using following chunk: 
rm(COPD)
rm(COPD1)
rm(COPD2)
rm(COPD3)
gc()


## ----------------------------------------------------------------------------------------------------------------------------------------
#rename new object to match what I previously used throughout code
COPD <- COPDm
#remove COPDm in terminal or using following chunk
rm(COPDm)


## ----------------------------------------------------------------------------------------------------------------------------------------
#save object without any further analysis
saveRDS(COPD, file = "C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\Weyer_W3_Practical_Rmark_files\\seuratobjraw.rds")


## ----------------------------------------------------------------------------------------------------------------------------------------
#running QC based on number of unique reads 
#first visualizing data as violin plot 
VlnPlot(COPD, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


## ----------------------------------------------------------------------------------------------------------------------------------------
#based on violin plots, will remove cells that have <200 or >10,000 unique feature counts. 
COPD <- subset(COPD, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)


## ----------------------------------------------------------------------------------------------------------------------------------------
#now to normalize data
COPD <- NormalizeData(COPD, normalization.method = "LogNormalize", scale.factor = 10000)


## ----------------------------------------------------------------------------------------------------------------------------------------
#identification of highly variable features. Seurat uses 2000 features by default
COPD <- FindVariableFeatures(COPD, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(COPD), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(COPD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


## ----------------------------------------------------------------------------------------------------------------------------------------
#from Seurat's website: "Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function: Shifts the expression of each gene, so that the mean expression across cells is 0, Scales the expression of each gene, so that the variance across cells is 1. This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate"
COPD <- ScaleData(COPD)


## ----------------------------------------------------------------------------------------------------------------------------------------
#running and visualizing PCA
COPD <- RunPCA(COPD, features = VariableFeatures(object = COPD))
print(COPD[["pca"]], dims = 1:5, nfeatures = 5)


## ----------------------------------------------------------------------------------------------------------------------------------------
VizDimLoadings(COPD, dims = 1:5, reduction = "pca")


## ----------------------------------------------------------------------------------------------------------------------------------------
#see how many PCs we have
DimPlot(COPD, reduction = "pca")


## ----------------------------------------------------------------------------------------------------------------------------------------
DimHeatmap(COPD, dims = 1:15, cells = 500, balanced = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------------
# We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
COPD <- JackStraw(COPD, dims = 50, num.replicate = 100)
COPD <- ScoreJackStraw(COPD, dims = 1:50)


## ----------------------------------------------------------------------------------------------------------------------------------------
#this takes a long time; elbow plot may be better for efficiency
JackStrawPlot(COPD, dims = 1:50)


## ----------------------------------------------------------------------------------------------------------------------------------------
ElbowPlot(COPD)
#elbow around PC 13 


## ----------------------------------------------------------------------------------------------------------------------------------------
COPD <- FindNeighbors(COPD, dims = 1:35)
COPD <- FindClusters(COPD, resolution = 0.5)


## ----------------------------------------------------------------------------------------------------------------------------------------
# Look at cluster IDs of the first 5 cells
head(Idents(COPD), 5)


## ----------------------------------------------------------------------------------------------------------------------------------------
# If you haven't installed UMAP, you can do so via:
#reticulate::py_install(packages = 'umap-learn')
COPD <- RunUMAP(COPD, dims = 1:35)


## ----------------------------------------------------------------------------------------------------------------------------------------
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(COPD, reduction = "umap")


## ----------------------------------------------------------------------------------------------------------------------------------------
#save your seurat object now
saveRDS(COPD, file = "C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\Weyer_W3_Practical_Rmark_files\\seuratobjprocessed.rds")


## ----------------------------------------------------------------------------------------------------------------------------------------
#finding markers for each cluster
COPD.markers <- FindAllMarkers(COPD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
COPD.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


## ----------------------------------------------------------------------------------------------------------------------------------------
write.csv(COPD.markers, "C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\Weyer_W3_Practical_Rmark_files\\COPD_markers.csv", row.names=FALSE)


## ----------------------------------------------------------------------------------------------------------------------------------------
#subsetting only CD45+ cells. decided this was the optimal way to go about this.
CD45 <- subset(COPD, PTPRC != 0)


## ----------------------------------------------------------------------------------------------------------------------------------------
#now to normalize data
CD45 <- NormalizeData(CD45, normalization.method = "LogNormalize", scale.factor = 10000)


## ----------------------------------------------------------------------------------------------------------------------------------------
#identification of highly variable features. Seurat uses 2000 features by default
CD45 <- FindVariableFeatures(CD45, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(CD45), 10)


## ----------------------------------------------------------------------------------------------------------------------------------------
#from Seurat's website: "Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function: Shifts the expression of each gene, so that the mean expression across cells is 0, Scales the expression of each gene, so that the variance across cells is 1. This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate"
CD45 <- ScaleData(CD45)


## ----------------------------------------------------------------------------------------------------------------------------------------
CD45 <- RunPCA(CD45, features = VariableFeatures(object = CD45))
DimPlot(CD45, reduction = "pca")


## ----------------------------------------------------------------------------------------------------------------------------------------
CD45 <- FindNeighbors(CD45, dims = 1:15)
CD45 <- FindClusters(CD45, resolution = 0.7)


## ----------------------------------------------------------------------------------------------------------------------------------------
CD45.markers <- FindAllMarkers(CD45, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CD45.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
write.csv(CD45.markers, "C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\Weyer_W3_Practical_Rmark_files\\CD45_markers.csv", row.names=FALSE)


## ----------------------------------------------------------------------------------------------------------------------------------------
CD45 <- RunUMAP(CD45, dims = 1:15)
png(filename = "CD45umap.png", width = 1200, height = 1200, units = "px", pointsize = 12, bg = "white")
DimPlot(CD45, reduction = "umap")
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------
#creating pngs of violin plots to determine cluster identity in more user-friendly space
#markers chosen based off what the human protein atlas uses
png(filename = "bcellVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("MS4A1",	"CD79A",	"BANK1",	"LTB",	"IGKC",	"MEF2C",	"BIRC3",	"STK17A",	"LIMD2",	"GPR183")) #B 
dev.off()

png(filename = "plasmaVlnPlotfeatures.png", width = 2000, height = 800, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("MZB1",	"JCHAIN",	"FKBP11",	"DERL3",	"ISG20",	"IGKC",	"ITM2C",	"CD79A",	"CD27",	"TNFRSF17")) #plasma
dev.off()

png(filename = "CD4TVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("IL32",	"CD3E",	"CD3D",	"CD2",	"LTB",	"ACAP1",	"TRAC",	"TRBC2")) #CD4T
dev.off()

png(filename = "CD8TVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("CCL5",	"IL32",	"CD3E",	"NKG7",	"CD2",	"GZMB",	"CTSW",	"CST7",	"GZMA",	"CD7")) #CD8T
dev.off()

png(filename = "pTVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("IL32",	"CD3D",	"TYMS",	"CD3E",	"KIAA0101",	"MKI67",	"RRM2",	"ANP32E",	"NUDT1",	"SMC4")) #pT
dev.off()

png(filename = "NKVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("NKG7",	"GNLY",	"GZMB",	"PRF1",	"CST7",	"KLRD1",	"CTSW",	"CD7",	"GZMA",	"CD247")) #NKs
dev.off()

png(filename = "DC1VlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("C1orf54",	"CPNE3", "LGALS2",	"CLEC9A",	"WDFY4",	"TAP1",	"NAAA",	"ICAM3",	"DNASE1L3",	"CADM1")) #DC1
dev.off()

png(filename = "DC2VlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("GPR183",	"CLEC10A",	"FGL2",	"RNASE6")) #Dc2
dev.off()

png(filename = "mDCVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("BIRC3",	"LAMP3",	"CCR7",	"CSF2RA",	"FSCN1",	"MARCKSL1",	"DAPP1",	"CCL22", "RASSF4",	"TMEM176A")) #Migratory DCs
dev.off()

png(filename = "pDCVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("GZMB",	"JCHAIN",	"ITM2C",	"IRF7",	"PPP1R14B",	"C12orf75",	"TSPAN13",	"UGCG",	"APP",	"TCF4")) #Plasmacytoid DCs
dev.off()

png(filename = "AMacroVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("FABP4",	"APOE",	"MSR1",	"SERPING1",	"MCEMP1",	"TREM1",	"GPNMB",	"CYP27A1",	"PPARG",	"NUPR1")) #Alveolar macrophages
dev.off()

png(filename = "AMaeCCL3VlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("MIR3945HG",	"CCL20",	"HLA-DQA2",	"PHLDA3",	"ITGB8",	"CCL4L2",	"NINJ1",	"SCGB3A1",	"TNIP3",	"CXCL5")) #Alveolar MœÜ CCL3+
dev.off()

png(filename = "AMaeMTVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("MT1X",	"MT1E",	"MT1G",	"SFTPB",	"KCNMA1",	"EMP1",	"MT1F",	"C15orf48",	"RP11-598F7.3")) #Alveolar MœÜ MT-positive
dev.off()

png(filename = "AMaeproflifVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("KIAA0101",	"CKS1B",	"TYMS",	"CENPW",	"NUSAP1",	"CDK1",	"UBE2C",	"PTTG1",	"MKI67", "SMC4")) #Alveolar MœÜ MT-positive
dev.off()

png(filename = "MDMaeVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("RNASE1",	"C15orf48",	"VCAN",	"CXCL8")) #Monocyte-derived MœÜ
dev.off()

png(filename = "perivascIMaeUVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("F13A1",	"RNASE1",	"STAB1",	"RNASE6",	"GPR183",	"FGL2",	"CSF1R",	"FCGR2B",	"IL13RA1", "METRNL")) #Interstitial MœÜ perivascular
dev.off()

png(filename = "cMonoVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("FCN1",	"VCAN",	"S100A12",	"CFP",	"EREG",	"CD300E",	"FGL2")) #Classical monocytes
dev.off()

png(filename = "ncMonoVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("FCN1",	"LILRB2",	"FGL2",	"CFP",	"LILRA5",	"CSF1R",	"LILRB1",	"WARS",	"FAM26F", "LIMD2")) #Interstitial MœÜ perivascular
dev.off()

png(filename = "MastVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("TPSAB1",	"TPSB2",	"CPA3",	"HPGDS",	"MS4A2",	"GATA2",	"KIT",	"CLU",	"VWA5A", "RGS13")) #Mast cells
dev.off()

png(filename = "mTVlnPlotfeatures.png", width = 2000, height = 1200, units = "px", pointsize = 12, bg = "white")
VlnPlot(CD45, features = c("CD3E",	"CCR7",	"CD27",	"SELL",	"CD69",	"ITGAE",	"FAS",	"IL2RB")) #memory t 
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------
png(filename = "toptenheatmap.png", width = 2000, height = 1200, units = "px", pointsize = 8, bg = "white")
CD45.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(CD45, features = top10$gene) + NoLegend()
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------
saveRDS(CD45, file = "C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\Weyer_W3_Practical_Rmark_files\\CD45subset.rds")


## ----------------------------------------------------------------------------------------------------------------------------------------
#annotate clusters
CD45 <- RenameIdents(object = CD45,
                               "0" = "CD8+ GZMB- NKG7- T cells 0", #CD3+CD8+CD161-
                               "1" = "CD4+ T cells 1", #CD3+CD4+CD161-
                               "2" = "DNT cells 2", #CD3+
                               "3" = "Monocyte-Derived Macrophages 3", #MSR1-
                               "4" = "MARCO+CD44+ unknowns 4",
                               "5" = "B cells 5",
                               "6" = "cD44+ Unknowns 6",
                               "7" = "CD4+ T cells 7", #CD3+CD4+CD161-
                               "8" = "Classical Monocytes", #MSR1+
                               "9" = "Effector Memory CD4+ T cells 9", #CD3+cD4+CD161+
                               "10" = "DNT cells 10",#CD3+
                               "11" = "DNT cells 11",#CD3+
                               "12" = "Unknowns 12",
                               "13" = "NK cells 13", 
                               "14" = "Alveolar Macrophage 14", #CD16+CD14+
                               "15" = "Effector Memory CD8+ T cells 15", #CD3+cD8+CD161+
			                         "16" = "Effector Memory CD8+ T cells 16", #CD3+cD8+CD161+ 
			                         "17" = "Alveolar Macrophages 17", #MSR1+
			                         "18" = "DNT cells 18" #CD3+
			                         )


## ----------------------------------------------------------------------------------------------------------------------------------------
plot <- DimPlot(CD45, reduction = "umap") + NoLegend()
LabelClusters(plot = plot, id = "ident")


## ----------------------------------------------------------------------------------------------------------------------------------------
png(filename = "umaplabelled.png", width = 2000, height = 1200, units = "px", pointsize = 18, bg = "white")
DimPlot(CD45, reduction = "umap",label = TRUE, label.size = 8.5, pt.size = 0.75)
dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------
saveRDS(CD45, file = "C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\Weyer_W3_Practical_Rmark_files\\CD45subset.rds")


## ----------------------------------------------------------------------------------------------------------------------------------------
purl("C:\\Users\\weyerce1\\Downloads\\CD45COPD_R_Analysis\\Weyer_W3_Practical_Rmark.Rmd")

