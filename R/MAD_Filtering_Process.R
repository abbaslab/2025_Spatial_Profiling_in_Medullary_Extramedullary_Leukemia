#####################################################################################
############################### MAD-Filtering #######################################
#####################################################################################

# Load necessary libraries
library(S4Vectors)
library(dplyr)
library(Seurat)
library(ggplot2)
library(SummarizedExperiment)
library(viridis)
library(ggExtra)
set.seed(1234)

# Load cell cycle gene data
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

##################################### srt_object #########################################
# Load spatial data
srt_object <- Load10X_Spatial("/path/to/10XGenomics/outs")

# Calculate mitochondrial and hemoglobin content
srt_object <- PercentageFeatureSet(srt_object, pattern = "^MT-", col.name = "percent_mito")
srt_object <- PercentageFeatureSet(srt_object, pattern = "^HB", col.name = "percent_hb")

# Generate violin plots before filtering
pdf("/path/to/output.pdf", width = 18, height = 14)
VlnPlot(srt_object, features = c("percent_mito"), pt.size = 0.1) + NoLegend() + ggtitle("Unfiltered Mito%")
VlnPlot(srt_object, features = c("percent_hb"), pt.size = 0.1) + NoLegend() + ggtitle("Unfiltered Hemoglobin %")
VlnPlot(srt_object, features = c("nFeature_Spatial"), pt.size = 0.1) + NoLegend() + ggtitle("Unfiltered nFeature")
VlnPlot(srt_object, features = c("nCount_Spatial"), pt.size = 0.1) + NoLegend() + ggtitle("Unfiltered nCount")

# Calculate log10 genes per UMI
srt_object$log10GenesPerUMI <- log10(srt_object$nFeature_Spatial) / log10(srt_object$nCount_Spatial)
metadata <- srt_object@meta.data
metadata$cells <- rownames(metadata)
metadata %>%
  ggplot(aes(x=nCount_Spatial, y=nFeature_Spatial, color=percent_mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  ggtitle("UMI")

################ Filtering by mitochondrial content ################
# Define MAD thresholds
max.mito.thr <- median(metadata$percent_mito) + 3*mad(metadata$percent_mito)
min.mito.thr <- median(metadata$percent_mito) - 3*mad(metadata$percent_mito)

# Plot with thresholds
p1 <- ggplot(metadata, aes(x=nFeature_Spatial, y=percent_mito)) +
  geom_point() +
  geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
  geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
  annotate(geom = "text", label = paste0(as.numeric(table(metadata$percent_mito > max.mito.thr | metadata$percent_mito < min.mito.thr)[2]), " cells removed\n",
                                         as.numeric(table(metadata$percent_mito > max.mito.thr | metadata$percent_mito < min.mito.thr)[1]), " cells remain"), x = 5000, y = 15)

ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100)

################ Filtering by gene and UMI counts ################
# Define MAD thresholds for genes and UMIs
min.Genes.thr <- median(log10(metadata$nFeature_Spatial)) - 3*mad(log10(metadata$nFeature_Spatial))
max.Genes.thr <- median(log10(metadata$nFeature_Spatial)) + 3*mad(log10(metadata$nFeature_Spatial))
max.nUMI.thr <- median(log10(metadata$nCount_Spatial)) + 3*mad(log10(metadata$nCount_Spatial))

# Gene/UMI scatter plot before filtering
p1 <- ggplot(metadata, aes(x=log10(nCount_Spatial), y=log10(nFeature_Spatial))) +
  geom_point() +
  geom_smooth(method="lm") +
  geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
  geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
  geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2)

ggMarginal(p1, type = "histogram", fill="lightgrey")

# Summary statistics
summary(srt_object$nFeature_Spatial)
summary(srt_object$nCount_Spatial)
summary(srt_object$percent_mito)
summary(srt_object$percent_hb)

# Apply filters
filtered_metadata <- metadata %>%
  filter(percent_mito <= max.mito.thr, percent_mito >= min.mito.thr) %>%
  filter(log10(nFeature_Spatial) <= max.Genes.thr, log10(nFeature_Spatial) >= min.Genes.thr) %>%
  filter(log10(nCount_Spatial) <= max.nUMI.thr)

# Update Seurat object with filtered data
srt_object_filtered <- subset(srt_object, cells = filtered_metadata$cells)

# Generate violin plots after filtering
VlnPlot(srt_object_filtered, "nCount_Spatial")
p1 <- VlnPlot(srt_object_filtered, "nFeature_Spatial") + ggtitle("MAD_Filtered nFeature")
p2 <- VlnPlot(srt_object, "nFeature_Spatial") + ggtitle("Unfiltered nFeature")

# Update coordinates in filtered object
srt_object_filtered@images$slice1@coordinates$tissue <- as.integer(srt_object_filtered@images$slice1@coordinates$tissue)
srt_object_filtered@images$slice1@coordinates$row <- as.integer(srt_object_filtered@images$slice1@coordinates$row)
srt_object_filtered@images$slice1@coordinates$col <- as.integer(srt_object_filtered@images$slice1@coordinates$col)
srt_object_filtered@images$slice1@coordinates$imagerow <- as.integer(srt_object_filtered@images$slice1@coordinates$imagerow)
srt_object_filtered@images$slice1@coordinates$imagecol <- as.integer(srt_object_filtered@images$slice1@coordinates$imagecol)

# Plot spatial features
SpatialFeaturePlot(srt_object_filtered, c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))

# SCTransform and Cell Cycle Scoring
srt_object_filtered <- SCTransform(srt_object_filtered, assay = "Spatial")
srt_object_filtered <- CellCycleScoring(srt_object_filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
SpatialDimPlot(srt_object_filtered, group.by = "Phase")
VlnPlot(srt_object_filtered, c("nCount_SCT", "nFeature_SCT"))
SpatialFeaturePlot(srt_object_filtered, c("nCount_SCT", "nFeature_SCT"))
srt_object_filtered <- FindVariableFeatures(srt_object_filtered, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(srt_object_filtered), 10)
plot1 <- VariableFeaturePlot(srt_object_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# PCA, Clustering, and UMAP
srt_object_filtered <- RunPCA(srt_object_filtered, assay = "SCT")
VizDimLoadings(srt_object_filtered, dims = 1:2, reduction = "pca")
DimPlot(srt_object_filtered, reduction = "pca")
DimHeatmap(srt_object_filtered, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(srt_object_filtered, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(srt_object_filtered, ndims = 50)

srt_object_filtered <- FindNeighbors(srt_object_filtered, reduction = "pca", dims = 1:25)
srt_object_filtered <- FindClusters(srt_object_filtered, resolution = seq(0, 1, 0.1), random.seed = 123)
srt_object_filtered <- FindClusters(srt_object_filtered, resolution = 0.14, random.seed = 123)
srt_object_filtered <- RunUMAP(srt_object_filtered, reduction = "pca", dims = 1:25)

# Feature and Spatial Plots
FeaturePlot(srt_object_filtered, c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
FeaturePlot(srt_object_filtered, c("nCount_SCT", "nFeature_SCT"))
SpatialDimPlot(srt_object_filtered)


# Read in pathology annotations
annotations <- read.csv("/path/to/pathology/annotation.csv", stringsAsFactors = FALSE)
srt_object_filtered[["Pathology"]] <- annotations$Leukemia[match(colnames(srt_object_filtered), annotations$Barcode)]
srt_object_filtered$Pathology[is.na(srt_object_filtered$Pathology)] <- "Marrow"

# Define colors for pathology annotation
my_cols <- c("Bone" = "#e3dac9", "Leuk-enriched" = "#ff681f", "Marrow" = "#91403d")
srt_object_filtered$Pathology <- factor(srt_object_filtered$Pathology, levels = names(my_cols))

# Plot pathology annotations
SpatialDimPlot(srt_object_filtered, group.by = "Pathology", crop = FALSE, alpha = 0.9, pt.size.factor = 1.6, cols = my_cols)
ggsave("/path/to/outputs/pathology_annotation.pdf")

# Close PDF device
dev.off()

# Save filtered Seurat object
saveRDS(srt_object_filtered, "/path/to/output/srt_object_filtered.RDS")
