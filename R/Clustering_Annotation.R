################################################################################
#                Integrated scRNA-seq & Spatial Analysis                       #
#                                                                              #
#                                                                              #
# This script performs:                                                        #
#  - Data import for scRNA-seq and Visium-ST samples                           #
#  - UMAP visualization of scRNA-seq reference                                 #
#  - Spatial clustering and pathology annotation                               #
#  - Supervised label transfer (Seurat) deconvolution                          #
#  - Heatmap analyses using SCP::GroupHeatmap                                  #
#  - Cross-sample comparison of deconvolution outputs                          #
#  - SpaCET-based malignant cell state analysis                                #
################################################################################

############################# 1. DATA IMPORT ###################################

required_packages <- c(
        "Seurat", "ggplot2", "dplyr", "SCP", "SpaCET",
        "ggthemes", "RColorBrewer"
)
lapply(required_packages, require, character.only = TRUE)

scRNA_Ref    <- readRDS("/path/to/scRNA/Reference.RDS")
BM1_filtered <- readRDS("/path/to/data/BM1_filtered.RDS")
BM2_filtered <- readRDS("/path/to/data/BM2_filtered.RDS")
EM1_filtered <- readRDS("/path/to/data/EM1_filtered.RDS")
EM2_filtered <- readRDS("/path/to/data/EM2_filtered.RDS")

################################################################################
############################# 2. scRNA UMAP ####################################

ColAssign <- function(Var, palettes = "Classic 20") {
        require(ggthemes); require(RColorBrewer)
        pal <- tableau_color_pal(palette = palettes, direction = 1, type = "regular")
        
        if (length(Var) > 20) {
                palOut <- colorRampPalette(pal(20))(length(Var))
        } else if (length(Var) == 20) {
                palOut <- pal(20)
        } else {
                palOut <- pal(20)
                palOut <- setdiff(palOut, c("#7F7F7F", "#C7C7C7"))
                palOut <- c(palOut, "#7F7F7F", "#C7C7C7")
                palOut <- palOut[1:length(Var)]
        }
        names(palOut) <- Var
        return(palOut)
}

DimPlot(
        scRNA_Ref,
        group.by = "class2",
        cols = ColAssign(unique(scRNA_Ref$class2)),
        raster = FALSE
) +
        NoAxes() + ggtitle("") +
        theme(
                legend.title = element_text(size = 0),
                legend.text  = element_text(size = 8),
                legend.key.height = unit(0.25, "cm"),
                legend.key.width  = unit(0.25, "cm")
        ) +
        NoAxes()

ggsave("/path/to/Figures/AllCells_DimPlot.pdf", height = 12, width = 18, unit = "cm")

################################################################################
########################### 3. SPATIAL CLUSTERS ################################

seurat_cols <- c(`0` = "#FFC312", `1` = "#C4E538", `2` = "#12CBC4")
BM1_filtered$seurat_clusters <- factor(BM1_filtered$seurat_clusters, levels = names(seurat_cols))

SpatialDimPlot(
        BM1_filtered,
        group.by = "seurat_clusters",
        cols = seurat_cols,
        crop = FALSE,
        pt.size.factor = 1.25,
        stroke = 0
)
ggsave("/path/to/Figures/SpatialMap_UnsupervisedClusters.pdf")

# Pathology annotation
annotations <- read.csv("/path/to/Pathology/BM1_Pathology.csv", stringsAsFactors = FALSE)
BM1_filtered[["Pathology"]] <- annotations$Regions[match(colnames(BM1_filtered), annotations$Barcode)]

region_cols <- c(
        `Region1` = "#ff7f0e",
        `Region2` = "#2ca02c",
        `Region3` = "#1f77b4"
)

BM1_filtered$Regions <- factor(BM1_filtered$Regions, levels = names(region_cols))

SpatialDimPlot(
        BM1_filtered,
        group.by = "Regions",
        crop = FALSE,
        pt.size.factor = 1.25,
        cols = region_cols,
        stroke = 0
)
ggsave("/path/to/Figures/SpatialMap_PathologyAnnotation.pdf")

################################################################################
######################## 4. SUPERVISED DECONVOLUTION ###########################

anchors <- FindTransferAnchors(
        reference = scRNA_Ref,
        query = BM1_filtered,
        normalization.method = "SCT",
        recompute.residuals = FALSE
)

predictions.assay <- TransferData(
        anchorset = anchors,
        refdata = scRNA_Ref$class2,
        prediction.assay = TRUE,
        weight.reduction = BM1_filtered[["pca"]],
        dims = 1:25
)

BM1_filtered[["predictions"]] <- predictions.assay

SpatialFeaturePlot(
        BM1_filtered,
        features = "AML",
        crop = FALSE,
        pt.size.factor = 1.25,
        stroke = 0
)
ggsave("/path/to/Figures/SpatialMap_AML_Deconvolution.pdf")

SpatialFeaturePlot(
        BM1_filtered,
        features = "Late Erythroid",
        crop = FALSE,
        pt.size.factor = 1.25,
        stroke = 0
)
ggsave("/path/to/Figures/SpatialMap_Erythroid_Deconvolution.pdf")

################################################################################
############################ 5. HEATMAP ANALYSIS ################################

ht <- GroupHeatmap(
        BM1_filtered,
        group.by = "Regions",
        cell_annotation = "seurat_clusters",
        features = c(
                "HBB","HBD","HBA1","HBA2","GATA1","GATA2",
                "S100A12","FCGR3A","CD14","CD33","MS4A7"
        ),
        group_palcolor = list(region_cols),
        add_dot = TRUE,
        add_bg = TRUE,
        flip = TRUE,
        cell_annotation_palcolor = list(c("#FFC312", "#C4E538", "#12CBC4")),
        show_column_names = TRUE,
        cell_annotation_params = list(width = unit(12, "mm"))
)

ht$plot
panel_fix(
        ht$plot,
        height = 4.5,
        width = 7,
        dpi = 500,
        raster = FALSE,
        save = "/path/to/Figures/PutativeMarkers_HeatPlot.pdf"
)

################################################################################
############################ 6. SAMPLE MERGING #################################

BM1_filtered$Sample  <- "BM1"
BM2_filtered$Sample  <- "BM2"
BM1_filtered$Patient <- "PT1"
BM2_filtered$Patient <- "PT2"

AddFeaturesFromAssay <- function(srt, assay, prefix) {
        metadata <- srt[[assay]]@data
        transposed_data <- as.data.frame(t(metadata))
        rownames(transposed_data) <- gsub("\\.", "-", rownames(transposed_data))
        colnames(transposed_data) <- paste0(prefix, colnames(transposed_data))
        srt <- AddMetaData(srt, metadata = transposed_data)
        return(srt)
}

BM1_filtered <- AddFeaturesFromAssay(BM1_filtered, assay = "predictions", prefix = "Seurat_")
BM2_filtered <- AddFeaturesFromAssay(BM2_filtered, assay = "predictions", prefix = "Seurat_")

BMs <- merge(
        BM1_filtered,
        y = BM2_filtered,
        add.cell.ids = c("BM1", "BM2"),
        project = "BoneMarrows"
)

cell_types_ordered <- c(
        "Seurat_HSC","Seurat_AML","Seurat_GMP","Seurat_Monocytes",
        "Seurat_cDC","Seurat_pDC","Seurat_Platelet.Mega","Seurat_CLP",
        "Seurat_B","Seurat_Plasma","Seurat_NK","Seurat_UnconvT",
        "Seurat_CD8.Naive","Seurat_CD4.Naive","Seurat_CD4.Treg",
        "Seurat_CD4.Mem","Seurat_CD4.Eff","Seurat_CD8.Eff",
        "Seurat_CD8.Mem","Seurat_Early.Erythoid","Seurat_Late.Erythroid"
)

ht <- GroupHeatmap(
        BMs,
        features = cell_types_ordered,
        group.by = "Sample",
        add_dot = TRUE,
        flip = TRUE,
        show_column_names = TRUE,
        column_names_side = "top",
        show_row_names = TRUE,
        group_palcolor = list(c("#56CCF2", "#6FCF97")),
        dot_size = unit(12, "mm")
)

ht
ggsave("/path/to/Figures/BoneMarrows_DotPlot_CellTypes.pdf", height = 5, width = 10)

################################################################################
############################### 7. EXTRAMEDULLARY ###############################

EM1_filtered$seurat_clusters <- factor(EM1_filtered$seurat_clusters, levels = names(seurat_cols))

SpatialDimPlot(
        EM1_filtered,
        group.by = "seurat_clusters",
        cols = seurat_cols,
        crop = FALSE,
        pt.size.factor = 1.25,
        stroke = 0
)
ggsave("/path/to/Figures/SpatialMap_UnsupervisedClusters_EM1.pdf")

annotations <- read.csv("/path/to/Pathology/EM1_Pathology.csv", stringsAsFactors = FALSE)
EM1_filtered[["Pathology"]] <- annotations$Pathology[match(colnames(EM1_filtered), annotations$Barcode)]

em1_cols <- c(
        "Sarcoma"  = "#ff793f",
        "Dermis"   = "#218c74",
        "Epidermis"= "#B53471",
        "Gland"    = "#1f77b4"
)

EM1_filtered$Pathology <- factor(EM1_filtered$Pathology, levels = names(em1_cols))

SpatialDimPlot(
        EM1_filtered,
        group.by = "Pathology",
        crop = FALSE,
        pt.size.factor = 1.25,
        cols = region_cols,
        stroke = 0
)
ggsave("/path/to/Figures/SpatialMap_PathologyAnnotation_EM1.pdf")

################################################################################
############################### 8. SpaCET ANALYSIS #############################

library(SpaCET)
source("/path/to/Customization/utilities.R")

EM1_spacet <- convert.Seurat(EM1_filtered)
EM1_spacet <- SpaCET.quality.control(EM1_spacet)

SpaCET.visualize.spatialFeature(
        EM1_spacet,
        spatialType = "QualityControl",
        spatialFeatures = c("UMI", "Gene"),
        imageBg = FALSE
)

EM1_spacet <- SpaCET.deconvolution(EM1_spacet, cancerType = "PANCAN", coreNo = 4)

SpaCET.visualize.spatialFeature(
        EM1_spacet,
        spatialType = "CellFraction",
        spatialFeatures = c("Malignant","Macrophage"),
        imageBg = FALSE
)

SpaCET.visualize.spatialFeature(
        EM1_spacet,
        spatialType = "CellFraction",
        spatialFeatures = "All",
        sameScaleForFraction = TRUE,
        pointSize = 0.1,
        nrow = 5,
        imageBg = TRUE
)

EM1_spacet <- SpaCET.deconvolution.malignant(EM1_spacet, coreNo = 4)
EM1_filtered <- addTo.Seurat(EM1_spacet, EM1_filtered)

DefaultAssay(EM1_filtered) <- "propMatFromSpaCET"

SpatialFeaturePlot(
        EM1_filtered,
        features = "Malignant",
        crop = FALSE,
        pt.size.factor = 1.25,
        cols = region_cols,
        stroke = 0
)
ggsave("/path/to/Figures/EM1_MalignantCellDeconvolution.pdf")

################################################################################
############################### 9. TISSUE HEATMAP ###############################

DefaultAssay(EM1_filtered) <- "SCT"

ht <- GroupHeatmap(
        EM1_filtered,
        group.by = "Pathology",
        cell_annotation = "seurat_clusters",
        features = c(
                "CD33","CD34","FCGR1A",        # Leukemia
                "COL1A1","FBN1","VWF",         # Dermis
                "KRT14","FLG","LOR",           # Epidermis
                "KRT79","KRT7","ATP1A1"        # Gland
        ),
        group_palcolor = list(em1_cols),
        add_dot = TRUE,
        add_bg = TRUE,
        flip = TRUE,
        cell_annotation_palcolor = list(c("#FFC312", "#C4E538", "#12CBC4")),
        show_column_names = TRUE,
        cell_annotation_params = list(width = unit(12, "mm"))
)

ht$plot
panel_fix(
        ht$plot,
        height = 4.15,
        width  = 8,
        dpi = 500,
        raster = FALSE,
        save = "/path/to/Figures/EM1_PutativeMarkersHeatMap.pdf"
)
