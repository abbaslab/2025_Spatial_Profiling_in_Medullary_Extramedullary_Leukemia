################################################################################
# Project: Spatial Transcriptomics & scRNA-seq Integration in AML
# Script:  Analysis Code for AML Niche, Inflammation, Deconvolution, and Signaling
# Notes:   This script is organized into logical analysis sections (HLS/LLS,
#          SpaCET-based colocalization, inflammatory composite scores,
#          CellChat-based communication, AUCell-based pathway scoring, and
#          AML-only scRNA-seq analysis).
################################################################################

########################## Analysis Block: HLS vs LLS ##########################
############################ Data Import #######################################
# Load Libraries
required_packages <- c("Seurat", "ggplot2", "dplyr", "SCP", "SpaCET", "ggthemes", "RColorBrewer", "ComplexHeatmap", "circlize", "grid")
lapply(required_packages, require, character.only = TRUE)

# Load RDS files
scRNA_Ref <- readRDS("/path/to/scRNA/Reference.RDS")
BM1_filtered <- readRDS("/path/to/data/BM1_filtered.RDS")
BM2_filtered <- readRDS("/path/to/data/BM2_filtered.RDS")
EM1_filtered <- readRDS("/path/to/data/EM1_filtered.RDS")
EM2_filtered <- readRDS("/path/to/data/EM2_filtered.RDS")

################################################################################
####### Section: Define AML-high (HLS) vs AML-low (LLS) Spots & Spatial Map ####
################################################################################

# Define spots that have higher leukemic score than the median AML score.
# The Seurat_AML feature originates from the previously defined Seurat-based deconvolution.
Pos_AML_Spots <- WhichCells(BM1_filtered, expression = AML > median(BM1_filtered$Seurat_AML))
# Classify spots into HLS (High Leukemic Score) vs LLS (Low Leukemic Score)
BM1_filtered$AML_Spots <- ifelse(colnames(BM1_filtered) %in% Pos_AML_Spots, "HLS", "LLS")

# Spatial map of HLS vs LLS spots
AML_Spots_col = c("HLS" ="#D8511D", "LLS" = "#212E52")
SpatialDimPlot(BM1_filtered, group.by = "AML_Spots", crop = F, pt.size.factor = 1.25, cols = AML_Spots_col)
ggsave("/path/to/Figures/HLS_LLS_Spots.pdf")

################################################################################
######## Section: Cluster Composition of HLS vs LLS (Stacked Bar / Trend) #######
################################################################################

# Stacked bar / trend plot showing the distribution of unsupervised clusters
# within HLS vs LLS groups; seurat_cols is defined earlier as cluster colors.
StatPlot(BM1_filtered@meta.data, stat.by = "seurat_clusters", group.by = "AML_Spots", plot_type = "trend", label = T, palcolor = seurat_cols)
ggsave("/path/to/Figures/HLS_LLS_vs_Clusters_StackedBarPlot.pdf")

################################################################################
######## Section: Differential Deconvolution Between HLS and LLS ###############
################################################################################

# Define assay as deconvolution scores
DefaultAssay(BM1_filtered) <- "predictions"

# Differential “deconvolution score” test between HLS and LLS using the
# deconvolution assay. Each cell type’s deconvolution score is compared with
# a Wilcoxon rank-sum test to infer spatial co-localization patterns.
BM1_filtered <- RunDEtest(BM1_filtered, group_by = "AML_Spots", only.pos = F)

# Volcano plot summarizing differential co-localization of cell types between HLS and LLS
p <- VolcanoPlot(BM1_filtered, group_by = "AML_Spots", nlabel = 5, pt.size = 6, stroke.highlight = 4, label.size = 4, sizes.highlight = 3)
p
ggsave("/path/to/Figures/BM1_CellCoLocalization.pdf", height = 8, width = 12)

################################################################################
######## Section: Spatial & Cluster-Level Patterns of Monocytes and GMP ########
################################################################################

# Spatial map of Monocyte deconvolution scores
SpatialFeaturePlot(BM1_filtered, features = "Seurat_Monocytes", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_Monocytes.pdf")

# Box plot of Monocyte scores across unsupervised clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_Monocytes", plot_type = "box", palcolor = seurat_cols, y.max = 1, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/Monocytes_Cluster.pdf")

# Spatial map of GMP deconvolution scores
SpatialFeaturePlot(BM1_filtered, features = "Seurat_GMP", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_SpatialMap_GMP.pdf")

# Box plot of GMP scores across unsupervised clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_GMP", plot_type = "box", palcolor = seurat_cols, y.max = 1, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/BM1_GMP_Clusters.pdf")

################################################################################
######## Section: Spatial & Cluster-Level Patterns of Erythroid & CD8 ##########
################################################################################

# Spatial map for late erythroid-like deconvolution scores
SpatialFeaturePlot(BM1_filtered, features = "Seurat_Late Erythroid", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_SpatialMap_LateErythroid.pdf")

# Box plot of late erythroid-like scores across clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_Late Erythroid", plot_type = "box", palcolor = seurat_cols, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/BM1_LateEryth_Box.pdf")

# Spatial map for CD8 naive-like deconvolution scores
SpatialFeaturePlot(BM1_filtered, features = "Seurat_CD8 Naive", crop = F, pt.size.factor = 1.25, stroke = 0)
ggsave("/path/to/Figures/BM1_SpatialMap_CD8_Naive.pdf")

# Box plot of CD8 naive-like scores across clusters
FeatureStatPlot(BM1_filtered, group.by = "seurat_clusters", stat.by = "Seurat_CD8 Naive", plot_type = "box", palcolor = seurat_cols, y.max = 1, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
ggsave("/path/to/Figures/BM1_CD8_Naive_Box.pdf")

################################################################################
######## Analysis Block: SpaCET-Based Co-localization Matrix (EM1) #############
################################################################################

# Get the current SpaCET cell types present in the assay
current_features <- rownames(EM1_filtered[["propMatFromSpaCET"]])

# Select SpaCET cell types to retain for correlation analysis
features_to_keep <- setdiff(current_features, c("Malignant", "Malignant cell state A", "Malignant cell state B", "cDC1 CLEC9A", "Macrophage other", "Macrophage M2", "cDC3 LAMP3", "Macrophage M1", "B cell switched memory", "B cell exhausted", "B cell non-switched memory", "B cell naive", "Tfh", "Th17", "Th1", "Th2", "cDC2 CD1C", "Unidentifiable"))
features_to_keep

# Subset the SpaCET deconvolution assay to the selected features
data_matrix <- GetAssayData(EM1_filtered[["propMatFromSpaCET"]], slot = "data")
data_matrix_subset <- data_matrix[features_to_keep, ]

# Attach subsetted SpaCET assay back to the object
DefaultAssay(EM1_filtered) <- "SCT"
EM1_filtered[["SpaCET_Subset"]] <- CreateAssayObject(data = data_matrix_subset)
DefaultAssay(EM1_filtered) <- "SpaCET_Subset"

# Extract deconvolution scores
deconv_scores <- GetAssayData(EM1_filtered[["SpaCET_Subset"]], slot = "data")

# Correlation matrix between SpaCET cell types
correlation_matrix <- cor(t(deconv_scores))

# Mean abundance per cell type
cell_type_abundance <- rowMeans(deconv_scores)

# Helper for annotating correlation strengths as stars
get_correlation_stars <- function(correlation_value) {
        if (abs(correlation_value) > 0.7) {
                return("***")
        } else if (abs(correlation_value) > 0.5) {
                return("**")
        } else {
                return("")
        }
}

# Build a matrix of significance stars from the correlation matrix
correlation_stars_matrix <- matrix(nrow = nrow(correlation_matrix), ncol = ncol(correlation_matrix))
for (i in 1:nrow(correlation_matrix)) {
        for (j in 1:ncol(correlation_matrix)) {
                correlation_stars_matrix[i, j] <- get_correlation_stars(correlation_matrix[i, j])
        }
}

# Color function for correlation heatmap (blue–white–red)
custom_col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("#053061", "#2166AC", "white", "#B2182B", "#67001F"))

# Barplot annotation summarizing mean abundance per cell type
abundance_annotation <- rowAnnotation(
        Abundance = anno_barplot(
                cell_type_abundance, 
                gp = gpar(fill = "darkgray", col = "black"),
                border = TRUE, 
                axis_param = list(side = "top"),
                width = unit(2, "cm")
        )
)

# Correlation heatmap with overlaid stars and abundance annotation
heatmap <- Heatmap(
        correlation_matrix,
        name = "Correlation",
        col = custom_col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 10, fontface = "bold"),
        heatmap_legend_param = list(
                title = "Correlation",
                title_gp = gpar(fontsize = 12, fontface = "bold"),
                labels_gp = gpar(fontsize = 10),
                legend_height = unit(4, "cm")
        ),
        right_annotation = abundance_annotation,
        top_annotation = HeatmapAnnotation(
                text = anno_text(colnames(correlation_matrix), rot = 45, just = "right", gp = gpar(fontsize = 10))
        ),
        border = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA, lwd = 0.2))
                if (correlation_stars_matrix[i, j] != "") {
                        grid.text(correlation_stars_matrix[i, j], x = x, y = y, gp = gpar(fontsize = 15, col = "black"))
                }
        }
)

# Save co-localization (correlation) matrix as PDF
pdf("/path/to/Figures/EM1_CoLocalization_Matrix.pdf", width = 14, height = 14)
draw(heatmap)
dev.off()

################################################################################
######## Section: SpaCET Scores for Macrophage and cDC in EM1 ##################
################################################################################

# Add SpaCET deconvolution scores as metadata features (prefixed with Spacet_)
# The AddFeaturesFromAssay() helper is defined earlier in the workflow.
EM1_filtered <- AddFeaturesFromAssay(EM1_filtered, assay = "propMatFromSpaCET", prefix = "Spacet_")

# Spatial map for SpaCET macrophage scores
SpatialFeaturePlot(EM1_filtered, features = "Spacet_Macrophage", crop = F, pt.size.factor = 1.25, alpha = c(0,7), stroke = 0)
ggsave("/path/to/Figures/Macrophage_Alpha_Spatial.pdf")

# Violin plot of SpaCET macrophage scores across clusters
p <- FeatureStatPlot(EM1_filtered, group.by = "seurat_clusters", stat.by = "Spacet_Macrophage", add_point = T, palcolor = seurat_cols, add_box = T, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
panel_fix(p, height = 4, width = 3.5, save = "/path/to/Figures/EM1_Macrophages_Clusters_VlnPlot.pdf")

# Spatial map for SpaCET cDC scores
SpatialFeaturePlot(EM1_filtered, features = "Spacet_cDC", crop = F, pt.size.factor = 1.25, alpha = c(0,7), stroke = 0)
ggsave("/path/to/Figures/cDC_Alpha_Spatial.pdf")

# Violin plot of SpaCET cDC scores across clusters
p <- FeatureStatPlot(EM1_filtered, group.by = "seurat_clusters", stat.by = "Spacet_cDC", add_point = T, palcolor = seurat_cols, add_box = T, comparisons = list(c('2','1'), c('1', '0'), c('1', 0)))
panel_fix(p, height = 4, width = 3.5, save = "/path/to/Figures/EM1_cDCs_Clusters_VlnPlot.pdf")

########################## Analysis Block: Inflammation #########################
############################ Data Import #######################################
# Load Libraries
required_packages <- c("Seurat", "ggplot2", "dplyr", "SCP", "SpaCET", "ggthemes", "RColorBrewer", "ComplexHeatmap", "circlize", "grid")
lapply(required_packages, require, character.only = TRUE)

# Load RDS files
scRNA_Ref <- readRDS("/path/to/scRNA/Reference.RDS")
BM1_filtered <- readRDS("/path/to/data/BM1_filtered.RDS")
BM2_filtered <- readRDS("/path/to/data/BM2_filtered.RDS")
EM1_filtered <- readRDS("/path/to/data/EM1_filtered.RDS")
EM2_filtered <- readRDS("/path/to/data/EM2_filtered.RDS")

################################################################################
######## Section: Composite Inflammatory Pathway Score (Jenks Classes) #########
################################################################################

################################################################################
############################# Composite Score Function ##########################
################################################################################
# Function to add composite inflammatory pathway scores and Jenks-based
# categories to a Seurat object, and to return both plots and the underlying
# coordinate/score data frame.
addCompositeScores <- function(seurat_obj, pathways_list, image_name) {
        # Extract spatial coordinates from the specified image
        coords <- seurat_obj@images[[image_name]]@coordinates
        
        # Extract pathway scores from metadata
        pathways <- lapply(pathways_list, function(x) seurat_obj@meta.data[[x]])
        
        # Normalize and combine pathway scores
        normalized_pathways <- lapply(pathways, function(x) rescale(x, to = c(0, 1)))
        composite_score <- rowMeans(do.call(cbind, normalized_pathways))
        
        # Combine coordinates and composite scores into a data frame
        data <- data.frame(
                x = coords$imagerow,
                y = coords$imagecol,
                composite_score = composite_score
        )
        
        # Apply Jenks natural breaks to classify the composite score
        jenks_breaks <- classInt::classIntervals(data$composite_score, n = 4, style = "jenks")
        data$composite_score_category <- cut(data$composite_score, 
                                             breaks = jenks_breaks$brks, 
                                             labels = c("Low", "Medium-Low", "Medium-High", "High"))
        
        # Align with Seurat cell order
        data$cell <- rownames(coords)
        rownames(data) <- data$cell
        data <- data[rownames(seurat_obj@meta.data), ]
        
        # Store composite score and category in metadata
        seurat_obj@meta.data$composite_score <- data$composite_score
        seurat_obj@meta.data$composite_score_category <- data$composite_score_category
        
        # Spatial feature plots for continuous score and discrete category
        p1 <- SpatialFeaturePlot(seurat_obj, features = "composite_score") + 
                ggtitle("Spatial Distribution of Composite Scores")
        
        p2 <- SpatialDimPlot(seurat_obj, group.by = "composite_score_category") + 
                ggtitle("Spatial Distribution of Composite Score Categories")
        
        return(list(seurat_obj = seurat_obj, composite_plot = p1, category_plot = p2, data = data))
}

# Inflammatory pathway list (Hallmark signatures)
pathways_list <- c("HALLMARK_INFLAMMATORY_RESPONSE", 
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                   "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                   "HALLMARK_COMPLEMENT", 
                   "HALLMARK_IL2_STAT5_SIGNALING")

# Run composite score function on BM1 spatial sample
result <- addCompositeScores(BM1_filtered, pathways_list, image_name = "slice1")
# result <- addCompositeScores(EM1_filtered, pathways_list, image_name = "slice1")
BM1_filtered <- result$seurat_obj

# Rotated spatial plotting (e.g., for orientation / layout adjustment)
data_rotated <- data.frame(
        x = result$data$y,
        y = -result$data$x,
        composite_score = data$composite_score,
        composite_score_category = data$composite_score_category
)

# Define region corresponding to high activity in rotated coordinates
top_right_rotated <- data_rotated[data_rotated$x > quantile(data_rotated$x, 0.75) & data_rotated$y > quantile(data_rotated$y, 0.75), ]

# Spatial heatmap of composite inflammatory pathway activity with Jenks classes
spatial_plot <- ggplot(data_rotated, aes(x = x, y = y, fill = composite_score_category)) +
        geom_tile(alpha = 0.8, width = 15, height = 15) +
        scale_fill_manual(values = c("Low" = "#2166AC", "Medium-Low" = "#67A9CF", "Medium-High" = "#D1E5F0", "High" = "#B2182B")) +
        labs(title = "Composite Inflammatory Pathway Activity",
             fill = "Composite Score") +
        theme_minimal() +
        coord_fixed() +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12)) +
        annotate("text", x = mean(top_right_rotated$x), y = mean(top_right_rotated$y), label = "High Pathway Activity", color = "black", size = 5, fontface = "bold")
# Print
spatial_plot
# Save the plot
ggsave("/path/to/Figures/SpatialInflammationClass.pdf")

################################################################################
######## Analysis Block: Cell-Cell Communication with CellChat #################
################################################################################
library(CellChat)

# Helper to assign deconvolution-based labels per spot using prediction assay
assignLabels <- function(object, prediction = "predictions") {
        pred <- object[[prediction]]@data
        pred <- pred[1:(nrow(pred)-1), ]
        # Label each spot based on the maximum prediction probability
        labels = rownames(pred)[apply(pred, 2, which.max)]
        names(labels) <- colnames(pred)
        object$labels <- factor(labels)
        Idents(object) <- "labels"
        return(object)
}

# Assign labels to BM1 based on Seurat prediction assay
BM1_filtered <- assignLabels(BM1_filtered, prediction = "predictions")

# Create normalized data input for CellChat from SCT assay
data.input = Seurat::GetAssayData(BM18_filtered, layer = "data", assay = "SCT")

# Metadata for CellChat (cell labels and sample identity)
meta = data.frame(labels = Seurat::Idents(BM18_filtered), samples = "labels", row.names = names(Seurat::Idents(BM18_filtered)))
meta$samples <- factor(meta$samples)
unique(meta$labels)
unique(meta$samples)

# Spatial coordinates for CellChat (10x Visium / spatial transcriptomics)
spatial.locs = Seurat::GetTissueCoordinates(BM18_filtered, scale = NULL, cols = c("imagerow", "imagecol")) 

# Spatial factors: convert pixels to micrometers using Visium spot calibration
scalefactors = jsonlite::fromJSON(txt = file.path("./data/BM/outs/spatial/", 'scalefactors_json.json'))
spot.size = 65 # theoretical spot diameter (µm) for Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0]) # should be ~100µm for Visium data

# Create CellChat object using spatial mode
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

# Set the ligand–receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# Select secreted signaling pathways for communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use

# Subset expression data to signaling genes and run CellChat pre-processing
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

# Compute communication probabilities and construct spatial communication network
ptm = Sys.time()

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                              contact.dependent = TRUE, contact.range = 100)

# Filter out low-cell-number groups and compute pathway-level communication
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
pathway_analysis <- netAnalysis_computeCentrality(cellchat)
pathway_analysis <- netAnalysis_signalingRole_network(pathway_analysis)

# Visualize signaling roles of each cell group
netVisual_signalingRole_scatter(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Aggregate communication network across all pathways
cellchat <- aggregateNet(cellchat)

# Heatmaps of interaction count and strength
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

# Bubble plots summarizing interactions from AML, Monocytes, GMP, etc.
# Sources:
## 1 = AML
## 11 = GMP
## 14 = Monocytes
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:25), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1,14), targets.use = c(1,11), remove.isolate = FALSE) # AML↔Monocyte
netVisual_bubble(cellchat, sources.use = c(1,11,14), targets.use = c(1,11,14), remove.isolate = FALSE) # AML–Monocyte–GMP triad

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Individual pathway visualization (e.g., CXCL)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Outgoing vs incoming signaling role heatmaps
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Example: interactions received by a specific cell group (e.g., Inflam.DC)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4,5,6,7), targets.use = 1, legend.pos.x = 15)
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")
plotGeneExpression(cellchat, signaling = "CD99", enriched.only = TRUE, type = "violin")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
gg2

# Re-plot signaling role heatmaps (incoming/outgoing)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Check available signaling pathways
cellchat@netP$pathways
pathways.show <- c("CXCL") 

par(mfrow=c(1,1), xpd = TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord-style aggregate visualization
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Spatial layout visualization of ligand–receptor network
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 1, vertex.size.max = 0.5, alpha.image = 0.2, vertex.label.cex = 3.5, point.size = 1.5, weight.scale = T)

# Save example spatial CXCL signaling plot
ggsave("/Path/to/Figures/BM18_CXCLPathSpatial.pdf")

# Define cell types of interest
cell_types_of_interest <- c("AML", "Monocytes", "GMP")

# Subset communication network to cells of interest (skeleton placeholder)
cellchat_subset@net <- cellchat@net[cells_of_interest, cells_of_interest, ]

#...

################################################################################
######## Analysis Block: AUCell & CXCL12–CXCR4 / PI3K / Inflammation ###########
################################################################################
library(AUCell)

# Helper to read GMT gene sets (e.g. from MSigDB)
load_genesets <- function(path, ignore_cols = 2){
        x <- scan(path, what="", sep="\n")
        y <- strsplit(x, "\t")
        # y=lapply(y,str_remove_all,'\"')
        names(y) <- sapply(y, `[[`, 1)
        for(i in 1:ignore_cols){
                y <- lapply(y, `[`, -1)
        }
        return(y)
}

# Function to run AUCell on Visium counts and add scores to metadata
run_aucell<-function(sobject,path_list){
        exprMatrix=sobject@assays$Spatial$counts
        cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=2, plotStats=TRUE,splitByBlocks=TRUE)
        for (i in 1:length(path_list)){
                geneset=load_genesets(path_list[[i]],ignore_cols=2)
                cells_AUC <- AUCell_calcAUC(geneset, cells_rankings)
                score_mat <- t(SummarizedExperiment::assay(cells_AUC, 'AUC'))
                colnames(score_mat)<-make.names(colnames(score_mat))
                sobject=AddMetaData(sobject, as.data.frame(score_mat))
        }
        return(sobject)
}

# Run AUCell for Hallmark gene sets on BM18 Visium sample
BM18_filtered <- run_aucell(BM18_filtered, "/path/to/Data/h.all.v2023.1.Hs.symbols.gmt.txt")

# Define CXCL12–CXCR4 gene pair pathway
CXCL_Path <- c("CXCL12", "CXCR4")

# Gene set list used for AUCell scoring; can be extended with additional signatures
cxcl_list <- list(
        # ChrisTest = ChrisTest,
        CXCL_Path = CXCL_Path
        # HIF1A_Path = HIF1A_Path,
        # Adipo_MSC = Adipo_MSC
)

# Run AUCell scoring for each gene set in cxcl_list
for (gene_set_name in names(cxcl_list)) {
        BM18_filtered <- calculate_AUCell(BM18_filtered, cxcl_list, gene_set_name)
        cat("Processed AUCell for", gene_set_name, "\n")
}

# Fetch Hallmark and AUCell scores for CXCL12–CXCR4, PI3K/AKT/mTOR, IFNγ, EMT, and composite inflammation
cxcl12_cxcr4_scores <- FetchData(BM18_filtered, vars = "AUCell_CXCL_Path")
pi3k_data <- FetchData(BM18_filtered, vars = "HALLMARK_PI3K_AKT_MTOR_SIGNALING")
ifn_gamma_data <- FetchData(BM18_filtered, vars = "HALLMARK_INTERFERON_GAMMA_RESPONSE")
emt_data <- FetchData(BM18_filtered, vars = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
inflammation <- FetchData(BM18_filtered, vars = "composite_score")

# Combine into a single data frame for correlation and visualization
combined_data <- data.frame(
        CXCL12_CXCR4 = cxcl12_cxcr4_scores,
        PI3K_AKT_mTOR = pi3k_data,
        IFN_Gamma = ifn_gamma_data,
        EMT = emt_data,
        inflammation_Score <- inflammation
)

# Ensure column names are correct and human-readable
colnames(combined_data) <- c("CXCL12_CXCR4", "PI3K_AKT_mTOR", "IFN_Gamma", "EMT", "inflammation_Score")

# Sanity checks for combined_data
print(head(combined_data))
print(nrow(combined_data))

# Correlations between CXCL12–CXCR4 score and selected Hallmark signatures
cor_pi3k <- cor.test(combined_data$CXCL12_CXCR4, combined_data$PI3K_AKT_mTOR, method = "pearson")
cor_ifn <- cor.test(combined_data$CXCL12_CXCR4, combined_data$IFN_Gamma, method = "pearson")
cor_emt <- cor.test(combined_data$CXCL12_CXCR4, combined_data$EMT, method = "pearson")

# Additional correlation: inflammation composite vs PI3K/AKT/mTOR
cor_pi3k <- cor.test(combined_data$inflammation_Score, combined_data$PI3K_AKT_mTOR, method = "pearson")

# Print correlation results
cat("Correlation between CXCL12-CXCR4 and PI3K/AKT/mTOR:\n")
print(cor_pi3k)
cat("\nCorrelation between CXCL12-CXCR4 and IFN Gamma:\n")
print(cor_ifn)
cat("\nCorrelation between CXCL12-CXCR4 and EMT:\n")
print(cor_emt)

# Scatter plot: CXCL12–CXCR4 vs PI3K/AKT/mTOR colored by inflammation class
p1 <- ggplot(combined_data, aes(x = CXCL12_CXCR4, y = PI3K_AKT_mTOR, color = InflammationClass)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", size = 1) +
        scale_color_manual(values = inflammation_colors, name = "Inflammation Class") +
        xlab("CXCL12-CXCR4 Co-expression Score") +
        ylab("PI3K/AKT/mTOR Signaling") +
        ggtitle("CXCL12-CXCR4 vs PI3K/AKT/mTOR Modulated by Inflammation") +
        theme_bw() +
        theme(
                plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                axis.title = element_text(size = 14, face = "bold"),
                axis.text = element_text(size = 12),
                legend.position = "right",
                legend.title = element_text(size = 12, face = "bold"),
                legend.text = element_text(size = 10)
        )

print(p1)
ggsave("/path/to/Figures/CXCL12_PI3K_Inflammation.pdf", height = 6, width = 7)

################################################################################
######## Section: Spatial CXCR4 & Pathway Maps #################################
################################################################################
# Correlation between CXCR4 and inflammation score can be inspected spatially.

# Spatial map for CXCR4 expression on EM1
SpatialFeaturePlot(EM1_filtered, features = "CXCR4", crop = F, pt.size.factor = 1.25, stroke = 0)

# Spatial map for PI3K/AKT/mTOR Hallmark pathway
SpatialFeaturePlot(BM1_filtered, features = "HALLMARK_PI3K_AKT_MTOR_SIGNALING", crop = F, pt.size.factor = 1.25, stroke = 0)
# Spatial map for EMT/trans-differentiation Hallmark pathway
SpatialFeaturePlot(BM1_filtered, features = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", crop = F, pt.size.factor = 1.25, stroke = 0)

################################################################################
######## Analysis Block: T cell Exhaustion / Senescence Gene Sets ##############
################################################################################

# T cell dysfunction markers (e.g., Desai et al., 2023)
cd8_dysfunction=c("RPL41", "TXNIP", "CD2", "PRF1", "PSMB9", "LIMD2", "LTB", "GIMAP1", "GIMAP4", "CX3CR1", "TRAF3IP3", "GIMAP7", "GIMAP5", "DENND2D", "CLIC1", "MYL12A", "S100A11", "RAC2", "PSMB10", "ARPC1B", "CORO1A", "PCBP1", "UCP2", "PLEK", "CD8A")
senescence <- c("CD27", "CD28", "IL2","GZMB", "B3GAT1", "KLRG1", "KLRC1", "KLRC2", "IL6", "TNF", "IFNG", "TIGIT", "HAVCR2", "CCL5", "CCL16", "CCL23")
ex_papers <- c("HNF1A", "GZMB", "CD28", "CD69", "TOX", "CD38", "CTLA4", "ENTPD1", "HAVCR2", "PDCD1", "HLA-DRA", "MKI67", "IL2RB")
ex_markers <- c("TOX", "CD38", "LAG3", "CTLA4", "TIGIT", "HAVCR2", "PDCD1", "CXCL13", "TNFRSF9", "TOX2", "ENTPD2", "LAYN","EOMES","HLA-DRA", "GZMB")
exhaustion_markers <- unique(c(ex_markers,ex_papers))

# Treg markers (e.g., Braun et al., 2021)
Treg <- c("FOXP3", "TNFRSF4", "TNFRSF18", "TBC1D4", "IL2RA", "CTLA4", "NGFRAP1", "CORO1B", "RTKN2", "LAYN", "AC017002.1", "CISH", "SLAMF1", "NCF4", "DNPH1", "F5", "LTB", "CTSC", "BATF", "ICA1", "TIGIT", "UGP2", "PIM2", "FBLN7", "CD4", "IL32", "IKZF2", "MIR4435-1HG", "GBP2", "CARD16", "PHTF2", "GPX1", "IL1R2", "GBP5", "S100A4", "PBXIP1", "GLRX", "CLEC7A", "TBC1D8", "SPOCK2", "RPS12", "RPS27", "RPL23A", "ZFP36L2", "RARRES3", "CD99", "GIMAP1", "APMAP", "LITAF", "CLEC2B", "APOBEC3G", "GSTP1", "GIMAP5", "ANKRD32", "LYST", "MCTP2", "AC069363.1", "THEMIS", "NELL2", "CRTAM", "TC2N", "F2R", "EOMES", "PRF1", "PLEK", "GIMAP7", "AOAH", "KLRG1", "GIMAP4", "ITM2C", "CMC1", "KLRD1", "GPR56", "IFNG", "LAG3", "HOPX", "ANXA1", "GZMB", "CST7", "CCL3", "GZMH", "CTSW", "CD8B", "GZMA", "GZMK", "CD8A", "CCL4", "CCL5", "NKG7", "RPS15A", "FASLG", "PTGDR", "RPS27A", "RPS4Y1", "RPS27L", "PARP8", "MGAT4A", "ACTB", "LYAR")

# Organized T cell gene sets used for AUCell scoring
T_groups <- list(
        Poonam_exhaustion = exhaustion_markers,
        CD8_Dysfunction = cd8_dysfunction,
        T_Senescence =senescence
)

# Run AUCell scoring on each T-cell gene set for BM1_filtered
for (gene_set_name in names(T_groups)) {
        BM1_filtered <- calculate_AUCell(BM1_filtered, T_groups, gene_set_name)
        cat("Processed AUCell for", gene_set_name, "\n")
}

# Heatmap of exhaustion/senescence/Treg scores across inflammatory categories
SCP::GroupHeatmap(S18_NONa, group.by = "composite_score_category_v2", features =  c("AUCell_Poonam_exhaustion", "AUCell_CD8_Dysfunction", "AUCell_T_Senescence", "AUCell_Treg"), exp_method = "zscore", row_names_side = "left", show_row_names = T, add_dot = T, add_bg = T, group_palcolor = list(c("Low" = "#2166AC", "Medium-Low" = "#67A9CF", "Medium-High" = "#F9C6AB", "High" = "#B2182B", "Unkown" = "#B2182B")), limits = c(-1.5,1.5))

########################## Analysis Block: AML scRNA ###########################
############################ Data Import #######################################
# Load Libraries
required_packages <- c("Seurat", "ggplot2", "dplyr", "SCP", "SpaCET", "ggthemes", "RColorBrewer", "ComplexHeatmap", "circlize", "grid")
lapply(required_packages, require, character.only = TRUE)

# Load RDS files
scRNA_Ref <- readRDS("/path/to/scRNA/Reference.RDS")
BM1_filtered <- readRDS("/path/to/data/BM1_filtered.RDS")
BM2_filtered <- readRDS("/path/to/data/BM2_filtered.RDS")
EM1_filtered <- readRDS("/path/to/data/EM1_filtered.RDS")
EM2_filtered <- readRDS("/path/to/data/EM2_filtered.RDS")

# Subset only AML cells from scRNA reference
aml <- subset(scRNA_Ref, subset = class2 == "AML")

# Function implementing the Seurat + Harmony pipeline for AML sub-clustering
analyze_seurat <- function(seurat_object, number_features=2500, res=0.8, md=0.01, k=20){
        seurat_object <- FindVariableFeatures(seurat_object,nfeatures = number_features)
        seurat_object <- ScaleData(seurat_object)
        seurat_object <- RunPCA(object = seurat_object, verbose = FALSE,features = VariableFeatures(seurat_object),npcs=50)
        seurat_object <- RunHarmony(seurat_object, group.by.vars="orig.ident",assay.use ="RNA", max.iter.cluster=500)
        seurat_object <- FindNeighbors(seurat_object, dims=1:30,reduction = "harmony",k.param = k)
        seurat_object <- FindClusters(seurat_object, resolution = seq(0,1,0.1), random.seed=123)
        seurat_object <- FindClusters(object = seurat_object,resolution=res,random.seed=123,graph.name = 'SCT_snn')
        seurat_object <- RunUMAP(seurat_object,reduction = "harmony",seed.use = 123,dims=1:50,
                                 a=0.5,b=1.5, verbose=FALSE)
        return(seurat_object)
}

# Run AML-specific clustering and embedding
aml <- analyze_seurat(aml)

# Fine-tune AML UMAP parameters
aml <- RunUMAP(aml,reduction = "harmony",seed.use = 123,dims=1:50,
               a=3,b=1.3, verbose=FALSE)

# Subset deconvolved Visium sample to AML-high spots only
BM1_AML <- subset(BM1_filtered, subset = AML > median(BM1_filtered$Seurat_AML))

################################################################################
######## Section: AML Region Deconvolution in BM18 Visium ######################
################################################################################
DefaultAssay(BM18_AML) <- "SCT"

anchors <- FindTransferAnchors(reference = aml, query = BM18_AML, normalization.method = "SCT", recompute.residuals = F)
predictions.assay <- TransferData(anchorset = anchors, refdata = aml$class4, prediction.assay = TRUE,
                                  weight.reduction = BM18_AML[["pca"]], dims = 1:25)
BM18_AML[["predictions4"]] <- predictions.assay

DefaultAssay(BM18_AML) <- "predictions4"

SpatialFeaturePlot(BM18_AML, features = rownames(BM18_AML))
ggsave("./figs/Playground/Biological/PT2/BM18/Predictions4_SpatialPlot.pdf")

################################################################################
######## Section: Spatial Distribution of AML States vs Inflammation ###########
################################################################################

# Example violin plot: spatial AML “Committed-like” states vs composite inflammation categories
p1 <- FeatureStatPlot(S18_Leukemia, group.by = "composite_score_category", stat.by = "Seurat_Committed-like", palcolor = list(c("Low" = "#2166AC", "Medium-Low" = "#67A9CF", "Medium-High" = "#F9C6AB", "High" = "#B2182B")), add_box = T, plot_type = "violin")

# Save the violin plot of AML Committed-like states across inflammatory categories
panel_fix(p1, height = 4, width = 3, dpi = 300, save = "/path/to/Figures/EM1_CommittedViolin.pdf")

################################################################################
################# Section: Proteomics Data Integration #########################
################################################################################
# Read CSV file for Opal intensities
BM1_OPAL <- read.csv('Opal_Intensities/BM1.csv')
#Assign barcodes to rownames
rownames(BM1_OPAL) <- BM1_OPAL$Barcode
# Define spot IDs
spot_ids <- BM1_OPAL$Barcode
# Subset seurat object for only Opal matched spots
BM1_OPAL_Subsetted <- subset(BM1, cells = spot_ids)
# Add Opal intensities into seurat object
BM1_OPAL_Subsetted <- AddMetaData(BM1_OPAL_Subsetted, metadata = BM1_OPAL)

SpatialFeaturePlot(BM1_OPAL_Subsetted, features = 'Mean.Arcsinh.CD68', crop = F, pt.size.factor = 2.3)

