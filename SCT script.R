library(Seurat)
library(dplyr)
library(Matrix)
library(g.data)
library(DropletUtils)
library(patchwork)
library(harmony)
library(cowplot)
library(SoupX)
library(glmnet)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(metap)
library(Seurat)
library(SeuratData)
library(patchwork)

devtools::install_github(repo = "samuel-marsh/scCustomize")

remotes::install_github(repo = "samuel-marsh/scCustomize")
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(qs)

# Read in 
# Read in `matrix.mtx`
matrix<-Matrix::readMM('E:/EV project data CellrangerCount_withIntrons/CellRangerMatrices/ScRNAseq/Control/matrix.mtx.gz')

genes <- read_tsv("E:/EV project data CellrangerCount_withIntrons/CellRangerMatrices/ScRNAseq/Control/features.tsv.gz", col_names = FALSE)
gene_ids <- genes$X2

cell_ids <- read_tsv("E:/EV project data CellrangerCount_withIntrons/CellRangerMatrices/ScRNAseq/Control/barcodes.tsv.gz", col_names = FALSE)$X1

# Make the column names as the cell IDs and the row names as the gene IDs
rownames(matrix) <- make.unique(gene_ids)
colnames(matrix) <- cell_ids

Control.data<-matrix

CONTROL <- CreateSeuratObject(counts = Control.data, project = "CTRL", min.cells = 3, min.features = 200)

CONTROL

# Set up control object
CONTROL$EV <- "CTRL"
CONTROL[["percent.mt"]] <- PercentageFeatureSet(CONTROL, pattern = "^mt-")

CONTROL <- subset(CONTROL, subset = nFeature_RNA > 500& percent.mt < 5)
CONTROL <- NormalizeData(CONTROL, verbose = FALSE)
CONTROL <- FindVariableFeatures(CONTROL, selection.method = "vst", nfeatures = 2000)

# Read in `EV matrix.mtx`
mat<-Matrix::readMM('E:/EV project data CellrangerCount_withIntrons/CellRangerMatrices/ScRNAseq/EV/matrix.mtx.gz')

genesmat <- read_tsv("E:/EV project data CellrangerCount_withIntrons/CellRangerMatrices/ScRNAseq/EV/features.tsv.gz", col_names = FALSE)
genemat_ids <- genesmat$X2

cellmat_ids <- read_tsv("E:/EV project data CellrangerCount_withIntrons/CellRangerMatrices/ScRNAseq/EV/barcodes.tsv.gz", col_names = FALSE)$X1

# Make the column names as the cell IDs and the row names as the gene IDs
rownames(mat) <- make.unique(genemat_ids)
colnames(mat) <- cellmat_ids
EV.data<-mat

EV <- CreateSeuratObject(counts = EV.data, project = "EV", min.cells = 3, min.features = 200)

EV$EV <- "EV"
EV[["percent.mt"]] <- PercentageFeatureSet(EV, pattern = "^mt-")
EV <- subset(EV, subset = nFeature_RNA > 500& percent.mt < 5)
EV <- NormalizeData(EV, verbose = FALSE)
EV <- FindVariableFeatures(EV, selection.method = "vst", nfeatures = 2000)

# install glmGamPoi
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
# install sctransform from Github
devtools::install_github("satijalab/sctransform", ref = "develop")


CONTROL <- SCTransform(CONTROL, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.3, verbose = FALSE)

DimPlot(CONTROL, reduction = "umap")

EV<- SCTransform(EV, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(EV = 30, verbose = FALSE)
utero.list <- list(ctrl = CONTROL, stim = EV)
features <- SelectIntegrationFeatures(object.list = utero.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = utero.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
utero.combined <- FindClusters(immune.combined.sct, resolution = 0.3)

p1 <- DimPlot(utero.combined, reduction = "umap", group.by = "EV")
p2 <- DimPlot(utero.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)



DimPlot(utero.combined, reduction = "umap", split.by = "EV")

#find conserved markers of each cluster
DefaultAssay(utero.combined) <- "RNA"

markers.to.plot <- c("Ehd3", "Plat", "Kdr", "Nrp1", "S100a4", "Vim", "Nphs1", "Podxl",
                     "Wt1", "Synpo","Plac8","Tnc","Col1a1","S100a8","S100a9","Nkg7","Thsd4", "Havcr1","Scnn1g",
                     "Aqp2","Kit","Upk1b","Nphs2","Ncam1",
                     "Eln","Pdgfrb","Aqp1","C1qa","C1qb","Ccl21")

DotPlot(utero.combined, features = markers.to.plot, cols = c("blue", "Red"), dot.scale = 8, split.by = "EV") +
  RotatedAxis()

save(utero.combined, file ="utero.integrated.Rda")


#rename clusters
utero.combined <- RenameIdents(utero.combined,`0` = "PTS5",`1` = "GENC", `2` = "PTS1", `3` = "PTS2", `4` = "TAL", `5` = "Fib", `6` = "PTS6", `7` = "DCT2", `8` = "PTS3", `9` = "DCT1", `10` = "CNT",`11` = "Cycling PT", `14` = "Immunecells", `13` = "DL-TAL", `12` = "Podocytes",`15` = "IC1", `16` = "JGA",`17` = "IC2", `18` = "PTS4",`19` = "Endo", `20` = "CDC-PC")

DimPlot(utero.combined, label = TRUE)
# Visualization
p1 <- DimPlot(utero.combined, reduction = "umap", group.by = "EV")
p2 <- DimPlot(utero.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(utero.combined, reduction = "umap", split.by = "EV")


load(file = "immune.integrated.Rda")
view(immune.integrated.Rda)




immune.list <- subset(utero.combined,idents = "immunecells")


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
utero.combined$celltype.EV<- paste(Idents(utero.combined), utero.combined$EV, sep = "_")
utero.combined$celltype <- Idents(utero.combined)
Idents(utero.combined) <- "celltype.EV"

utero.combined <- PrepSCTFindMarkers(utero.combined)


#Immune cells cluster
Podocytes <- FindMarkers(utero.combined, assay= "SCT", ident.1 = "Podocytes_EV", ident.2 ="Podocytes_CTRL", verbose = FALSE)
head(Podocytes, n = 15)

rownames(Podocytes[Podocytes$avg_log2FC < -0.2 | Podocytes$avg_log2FC > 0.2  ,]) # deutliche Expression

View(Podocytes)
write.table(Podocytes, file="PodocyteslimmmaSCT1.txt", dec =",", sep = "\t")


#Immune cells cluster
CyclingPT <- FindMarkers(utero.combined, assay= "SCT", ident.1 = "CNT_EV", ident.2 ="CNT_CTRL", verbose = FALSE)
head(CyclingPT, n = 15)

rownames(CyclingPT[CyclingPT$avg_log2FC < -0.2 | CyclingPT$avg_log2FC > 0.2  ,]) # deutliche Expression

View(Podocytes)
write.table(CyclingPT, file="CyclingPTlimmmaSCT1.txt", dec =",", sep = "\t")

























immune.combined.sct.subset <- subset(immune.combined.sct, idents = c("B_STIM", "B_CTRL"))
immune.list <- subset(utero.combined,idents = c("Immunecells_EV", "Immunecells_CTRL"))
immune.list <- SplitObject(immune.list, split.by = "orig.ident")
immune.list <- lapply(X = immune.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = immune.list, nfeatures = 3000)
immune.list <- PrepSCTIntegration(object.list = immune.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = immune.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.integrated<-ScaleData(immune.combined.sct)
immune.integrated <- RunPCA(object = immune.integrated, npcs = 30)
immune.integrated <- RunUMAP(object = immune.integrated, dims = 1:10)
immune.integrated <- FindNeighbors(immune.integrated, dims = 1:10)
immune.integrated <- FindClusters(immune.integrated, resolution = 0.5)
# Visualization
p1immune <- DimPlot(immune.integrated, reduction = "umap", group.by = "EV")
p2immune <- DimPlot(immune.integrated, reduction = "umap", label = TRUE)
plot_grid(p1immune, p2immune)

DimPlot(immune.integrated, reduction = "umap", split.by = "EV")

save(immune.integrated, file = 'immune.renamedintegrated.Rda')

DefaultAssay(immune.integrated) <- "integrated"
#find conserved markers of each cluster
DefaultAssay(immune.integrated) <- "RNA"

markers.to.plot <- c("Ehd3", "Plat", "Kdr", "Nrp1", "S100a4", "Vim", "Nphs1", "Podxl",
                     "Wt1", "Synpo","Plac8","Tnc","Col1a1","S100a8","S100a9","Nkg7","Thsd4", "Nos1","Scnn1g",
                     "Aqp2","Kit","Upk1b","Nphs2","Ncam1",
                     "Eln","Pdgfrb","Aqp1","C1qa","C1qb")

DotPlot(immune.integrated, features = markers.to.plot, cols = c("blue", "Red"), dot.scale = 8, split.by = "EV") +
  RotatedAxis()


markers.to.plot <- c('Il1b','C1qb','Nkg7','Mertk','Tyrobp','Casp3','Clec10a','Hdc','Mmp9','Arg1','Adgre1','Mrc1','Cd68','Mrp8','Mpo',"S100a8",'Lcn2','Ly6c2','Ccr2','Cx3cr1','C1qa','Msr1','Ptprc','Cd3e','Cd4','Il7r','Cd8a','Mki67','Pcna','Itgax','Cd209a','Cd79a','Cd74')

DotPlot(immune.integrated, features = markers.to.plot, cols = c("blue", "Red"), dot.scale = 8, split.by = "EV") +
  RotatedAxis()


#rename clusters
immune.integrated <- RenameIdents(immune.integrated,`0` = "Mac_1",`1` = "Inf.Mac",`2` = "DC",`3` = "pDC",`4` = "T-cells")

DimPlot(immune.integrated, label = TRUE)

install.packages('BiocManager')
BiocManager::install('limma')


remove.packages("limma")



















library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
immune.integrated$celltype.EV <- paste(Idents(immune.integrated), immune.integrated$EV, sep = "_")
immune.integrated$celltype <- Idents(immune.integrated)
Idents(immune.integrated) <- "celltype.EV"


View(immune.integrated)

immune.integrated <- PrepSCTFindMarkers(immune.integrated)

#Macrophage cluster
InflmmatoryMac <- FindMarkers(immune.integrated, assay= "SCT", ident.1 = "Inf.Mac_EV", ident.2 ="Inf.Mac_CTRL", verbose = FALSE)
head(InflmmatoryMac, n = 15)

rownames(InflmmatoryMac[InflmmatoryMac$avg_log2FC < -0.2 | InflmmatoryMac$avg_log2FC > 0.2  ,]) # deutliche Expression

View(InflmmatoryMac)
write.table(InflmmatoryMac, file="InflmmatoryMac.txt", dec =",", sep = "\t")


#Mac_1 cluster
Mac_1 <- FindMarkers(immune.integrated, assay= "SCT", ident.1 = "Mac_1_EV", ident.2 ="Mac_1_CTRL", verbose = FALSE)
head(Mac_1 , n = 15)

rownames(Mac_1[Mac_1$avg_log2FC < -0.2 | Mac_1$avg_log2FC > 0.2  ,]) # deutliche Expression

View(Mac_1)
write.table(Mac_1, file="Mac_1.txt", dec =",", sep = "\t")

#AntiInfl.Mac cluster
AntiInfl.MAC <- FindMarkers(immune.integrated, ident.1 = "AntiInfl.Mac_EV", ident.2 ="AntiInfl.Mac_CTRL", verbose = FALSE)
head(AntiInfl.MAC , n = 15)

rownames(AntiInfl.MAC[AntiInfl.MAC$avg_log2FC < -0.4 | AntiInfl.MAC$avg_log2FC > 0.4  ,]) # deutliche Expression

View(AntiInfl.MAC)
write.table(AntiInfl.MAC, file="AntiInfl.MAC.txt", dec =",", sep = "\t")


#Dendritic cells cluster
pDC <- FindMarkers(immune.integrated, ident.1 = "pDC_EV", ident.2 ="pDC_CTRL", verbose = FALSE)
head(pDC, n = 15)

rownames(pDC[pDC$avg_log2FC < -0.4 | pDC$avg_log2FC > 0.4  ,]) # deutliche Expression

View(pDC)
write.table(pDC, file="pDC.txt", dec =",", sep = "\t")


#T-cells cluster
tcells <- FindMarkers(utero.combined, assay= "SCT",ident.1 = "CDC-PC_EV", ident.2 ="CDC-PC_CTRL", verbose = FALSE)
head(tcells, n = 15)

rownames(tcells[tcells$avg_log2FC < -0.2 | tcells$avg_log2FC > 0.2  ,]) # deutliche Expression

View(tcells)
write.table(tcells, file="t-cellsSCT.txt", dec =",", sep = "\t")



#B-cells cluster
bcells <- FindMarkers(immune.integrated, ident.1 = "B-cells_EV", ident.2 ="B-cells_CTRL", verbose = FALSE)
head(bcells, n = 15)

rownames(bcells[bcells$avg_log2FC < -0.4 | bcells$avg_log2FC > 0.4  ,]) # deutliche Expression

View(bcells)
write.table(bcells, file="b-cells.txt", dec =",", sep = "\t")


#Novel cluster
Novel <- FindMarkers(immune.integrated, ident.1 = "Novel_EV", ident.2 ="Novel_CTRL", verbose = FALSE)
head(Novel, n = 15)

rownames(Novel[Novel$avg_log2FC < -0.4 | Novel$avg_log2FC > 0.4  ,]) # deutliche Expression

View(Novel)


features <- c("Nlrp3", "Il1b", "Casp8", "Tlr4", "Myd88", "Nfkib","Casp1")
DoHeatmap(subset(utero.combined, downsample = 100), features = features, size = 3)

VlnPlot(utero.combined, features = features)



pal <- viridis(n = 10, option = "C", direction = -1)

gene_list_plot <- c("Nlrp3")
human_colors_list <- c("dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid",
                       "orange", "gold", "gray")

# Create Plots
Stacked_VlnPlot(seurat_object = utero.combined, features = gene_list_plot, x_lab_rotate = TRUE,
                )


FeaturePlot(object = utero.combined, features = "Jun")


FeaturePlot(object = utero.combined, features = c("P2ry12", "Fcrls", "Trem2", "Tmem119", "Cx3cr1", "Hexb", "Tgfbr1", "Sparc", "P2ry13",
                                                  "Olfml3", "Adgrg1", "C1qa", "C1qb", "C1qc", "Csf1r", "Fcgr3", "Ly86", "Laptm5"), cols = pal, order = T)

FeaturePlot_scCustom(seurat_object = utero.combined, features = features, split.by = "EV",
                     num_columns = 2)

Stacked_VlnPlot(seurat_object = utero.combined, features = c("P2ry12", "Fcrls", "Trem2", "Tmem119", "Cx3cr1", "Hexb", "Tgfbr1", "Sparc", "P2ry13",
                                                             "Olfml3", "Adgrg1", "C1qa", "C1qb", "C1qc", "Csf1r", "Fcgr3", "Ly86", "Laptm5"), x_lab_rotate = TRUE)


PC_Plotting(seurat_object = utero.combined, assay= "SCT", dim_number = 2)




#Immune cells cluster
Immunecells <- FindMarkers(utero.combined, ident.1 = "Immunecells_EV", ident.2 ="Immunecells_CTRL", verbose = FALSE)
head(Immunecells, n = 45)

rownames(Immunecells[Immunecells$avg_log2FC < -0.1 | Immunecells$avg_log2FC > 0.1  ,]) # deutliche Expression

View(Immunecells)
write.table(Immunecells, file="ImmunecellsL2F.txt", dec =",", sep = "\t")
