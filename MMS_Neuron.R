# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 11/16/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

# LOAD LIBRARIES ####
# Restart Rstudio or R

install.packages('ggplot2')
install.packages('cowplot')
install.packages('Matrix')
install.packages('ggridges')
install.packages('ggrepel')
install.packages('dplyr')
#install.packages('Seurat')
install.packages('plotly')
install.packages('clustree')
install.packages('patchwork')
install.packages('future')
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
BiocManager::install("EnhancedVolcano")
BiocManager::install("DoubletFinder")
BiocManager::install("glmGamPoi")
BiocManager::install("GOSemSim")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationHub")
BiocManager::install("GenomeInfoDb")
BiocManager::install("MeSHDbi")
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("DOSE")
BiocManager::install("dittoSeq")
BiocManager::install("escape")
BiocManager::install("ComplexHeatmap")
BiocManager::install(c("DropletUtils", "Nebulosa"))
BiocManager::install("hdf5r", force = TRUE)
BiocManager::install('multtest')
BiocManager::install("MAST")
BiocManager::install("enrichplot")

# install Seurat from Github (automatically updates sctransform)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

devtools::install_github("satijalab/seurat", ref = "develop")
devtools::install_github("satijalab/sctransform", ref = "develop", force = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("yanlinlin82/ggvenn")
devtools::install_github("gaospecial/ggVennDiagram")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages("harmony")
BiocManager::install("EnrichmentBrowser")
install.packages('SoupX')
install.packages('tidyverse')
install.packages("viridis")
install.packages("circlize")
install.packages("scCustomize")
install.packages("devtools")
install.packages("archive")
install.packages("R.utils")
install.packages("qs")
install.packages('metap')
install.packages('magick')
install.packages("corrplot")
install.packages("RColorBrewer")
install.packages("sunburstR")
install_github("didacs/ggsunburst")
install.packages("ggpubr")


# Run the following code once you have Seurat installed
suppressWarnings(
  {
    library(leiden)
    library(stringr)
    library(hdf5r)
    library(SoupX)
    library(Rcpp)
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(tidyverse)
    library(data.table)
    library(reticulate)
    library(Seurat)
    library(monocle3)
    library(harmony)
    library(Signac)
    library(EnsDb.Hsapiens.v86)
    library(GenomeInfoDb)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
    library(glmGamPoi)
    library(GOSemSim)
    library(org.Hs.eg.db)
    library(AnnotationHub)
    library(MeSHDbi)
    library(clusterProfiler)
    library(DOSE)
    library(dittoSeq)
    library(escape)
    library(EnrichmentBrowser)
    library(viridisLite)
    library(viridis)
    library(ComplexHeatmap)
    library(circlize)
    #library(scCustomize)
    library(Nebulosa)
    library(DropletUtils)
    library(ggvenn)
    library(ggVennDiagram)
    library(devtools)
    library(R.utils)
    library(qs)
    library(multtest)
    library(metap)
    library(MAST)
    library(magick)
    library(enrichplot)
    library(corrplot)
    library(DESeq2)
    library(RColorBrewer)
    library(sunburstR)
    library(d3r)
    library(ggpubr)
  }
)


# Set global environment parameter par-proc
#options(future.globals.maxSize = 8000 * 1024^2)
set.seed(1234)

# Python env
if(.Platform$OS.type == "windows") Sys.setenv(PATH= paste("C:/Users/mqadir/AppData/Local/r-miniconda/envs/r-reticulate",Sys.getenv()["PATH"],sep=";"))
py_config()

# WD
setwd(r"(C:\Users\mqadir\Box\Lab 2301\1. R_Coding Scripts\Sex Biology Study\WD)")
(WD <- getwd())
session_info()
sessionInfo()
# Check package versions
packageVersion("clusterProfiler")
packageVersion("dittoSeq")
packageVersion("escape")
packageVersion("Seurat")
packageVersion("signac")
packageVersion("EnrichmentBrowser")
packageVersion("org.Hs.eg.db")
packageVersion("DESeq2")

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
neuron.panc.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd shared with MMS\data\firstrun_raw)")

# STEP 2: Create Seurat objects ####
neuron.panc <- CreateSeuratObject(counts = neuron.panc.data, min.features = 500)

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
neuron.panc[["percent.mt"]] <- PercentageFeatureSet(object = neuron.panc, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(neuron.panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
neuron.panc <- subset(x = neuron.panc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 8)

# Normalise data
neuron.panc <- NormalizeData(neuron.panc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable genes
neuron.panc <- FindVariableFeatures(neuron.panc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(neuron.panc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(neuron.panc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling
all.genes <- rownames(neuron.panc)
neuron.panc <- ScaleData(neuron.panc, features = all.genes)

# Dimensions
neuron.panc <- RunPCA(neuron.panc, features = VariableFeatures(object = neuron.panc))
neuron.panc <- FindNeighbors(neuron.panc, dims = 1:10)
neuron.panc <- FindClusters(neuron.panc, resolution = 0.3)
ElbowPlot(neuron.panc)

# UMAP
neuron.panc <- RunUMAP(neuron.panc, dims = 1:5)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(neuron.panc, reduction = "umap")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
all.markers.neuron <- FindAllMarkers(neuron.panc, only.pos = TRUE)
all.markers.neuron %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Plot genes
FeaturePlot(neuron.panc, features = c("Nos1"))

# Save genes
write.csv(all.markers.neuron, file = r"(C:\Users\mqadir\Box\Fahd shared with MMS\data\all.markers.neuron.csv)")

# Subset cells
Idents(neuron.panc) <- "seurat_clusters"
nos1.cells <- subset(neuron.panc, idents = "3")

# re-analyze
# Dimensions
nos1.cells <- RunPCA(nos1.cells, features = VariableFeatures(object = nos1.cells))
nos1.cells <- FindNeighbors(nos1.cells, dims = 1:10)
nos1.cells <- FindClusters(nos1.cells, resolution = 0.5)
ElbowPlot(nos1.cells)

# UMAP
nos1.cells <- RunUMAP(nos1.cells, dims = 1:5)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(nos1.cells, reduction = "umap")

# Plot genes
FeaturePlot(nos1.cells, features = c("Nos1"))

# Selected genes
markers.to.plot <- c("Nos1", "Chat",  "Vip", "Gal",
                     "Chrm2", "Chrna3", "Chrnb4", "Htr3a", "Htr2b", "Htr2c", "Glp2r", "Glp1r", "Gipr")

# Dotplot
DotPlot(nos1.cells,  
        dot.scale = 8,
        col.min = -1, #minimum level
        col.max = 1,  #maximum level
        features = rev(markers.to.plot)) + 
  geom_point(aes(size=pct.exp), shape = 21, stroke=0.5) +
  theme_light() +
  #facet_wrap(~??? what metadata should be here??)
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.3, hjust=1, size =12, face = "bold", colour = "black")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title=element_text(size=12, face = "bold"), 
        legend.text=element_text(size=12, face = "bold")) +
  scale_colour_gradient2(low =c("dodgerblue"), mid = c("white"), high =c("red3")) +
  guides(color = guide_colorbar(title = 'Average Expression'))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
all.markers.nos1.cells <- FindAllMarkers(nos1.cells, only.pos = TRUE)
all.markers.nos1.cells %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Plot genes
FeaturePlot(nos1.cells, features = c("Nos1"))

# Save genes
write.csv(all.markers.nos1.cells, file = r"(C:\Users\mqadir\Box\Fahd shared with MMS\data\all.markers.nos1.cells.csv)")

