########################## Supplementary QC Script ##############################
############################ Data Import ########################################

# Load 10X Genomics SpaceRanger output directories as Seurat spatial objects
BM1 <- Load10X_Spatial("/path/to/BM1/outs/")
BM2 <- Load10X_Spatial("/path/to/BM2/outs/")
EM1 <- Load10X_Spatial("/path/to/EM1/outs/")
EM2 <- Load10X_Spatial("/path/to/EM2/outs/")
EM1_v1 <- Load10X_Spatial("/path/to/EM1_v1/outs/")
EM2_v1 <- Load10X_Spatial("/path/to/EM2_v1/outs/")

#################################################################################
########################### Section: Median nCount vs nFeature ##################
#################################################################################

# Load required libraries for correlation and visualization
library(ggplot2)
library(ggpubr)  

# Construct per-sample data frames for sequencing depth (nCount_Spatial)
# and feature complexity (nFeature_Spatial)
sample_data <- list(
        BM1 = data.frame(x_var = BM1@meta.data[,"nCount_Spatial"], y_var = BM1@meta.data[,"nFeature_Spatial"]),
        BM2 = data.frame(x_var = BM2@meta.data[,"nCount_Spatial"], y_var = BM2@meta.data[,"nFeature_Spatial"]),
        EM2 = data.frame(x_var = EM2@meta.data[,"nCount_Spatial"], y_var = EM2@meta.data[,"nFeature_Spatial"]),
        EM1 = data.frame(x_var = EM1@meta.data[,"nCount_Spatial"], y_var = EM1@meta.data[,"nFeature_Spatial"]),
        EM2_v1 = data.frame(x_var = EM2_v1@meta.data[,"nCount_Spatial"], y_var = EM2_v1@meta.data[,"nFeature_Spatial"]),
        EM1_v1 = data.frame(x_var = EM1_v1@meta.data[,"nCount_Spatial"], y_var = EM1_v1@meta.data[,"nFeature_Spatial"])
)

# Prepare containers for per-sample summary points (medians) and colors
representative_points <- list()
colors <- rainbow(length(sample_data))

# Compute per-sample median depth and complexity as representative summary points
for (i in 1:length(sample_data)) {
        df <- sample_data[[i]]
        representative_point <- data.frame(
                Median_x = median(df$x_var, na.rm = TRUE),
                Median_y = median(df$y_var, na.rm = TRUE),
                Sample = names(sample_data)[i]
        )
        representative_points[[i]] <- representative_point
}

# Combine summary statistics into one data frame for plotting
combined_data <- do.call(rbind, representative_points)

# Correlation plot of median nCount_Spatial vs median nFeature_Spatial across samples
gg <- ggplot(combined_data, aes(x = Median_x, y = Median_y, color = Sample)) +
        geom_point(size = 4) +
        geom_smooth(method = "lm", se = FALSE, color = "darkblue", size = 1) +
        xlab("Median nCount_Spatial") +
        ylab("Median nFeature_Spatial") +
        ggtitle("Correlation between Median nCount_Spatial and Median nFeature_Spatial") +
        theme_minimal() +
        theme(text = element_text(size = 14))

# Add Pearson correlation p-value annotation
gg <- gg + stat_cor(aes(label = paste("p =", ..p..)), method = "pearson", size = 4, color = "black")

gg

# Save high-resolution correlation plot
ggsave("./figs/Playground/All_correlation_plot.pdf", gg, width = 8, height = 6, dpi = 700)

# Display the figure (useful in interactive sessions)
print(gg)

#################################################################################
########################### Section: DV200 vs Mito Radar Plot ###################
#################################################################################

# Compute mean mitochondrial percentage per sample to compare with DV200
mito_percentages <- c(
        mean(EM2@meta.data$percent_mito),
        mean(EM1@meta.data$percent_mito),
        mean(BM2@meta.data$percent_mito),
        mean(BM1@meta.data$percent_mito)
)

# DV200 values for EM2, EM1, BM2, BM1 (RNA quality metrics)
dv200_values <- c(24, 53, 62, 39)

# Assemble data frame for radar plot (DV200 and mitochondrial content)
data_for_radar <- data.frame(
        DV200 = dv200_values,
        Mito_Percentage = mito_percentages
)

rownames(data_for_radar) <- c("EM2", "EM1", "BM2", "BM1")
data_for_radar_t <- t(data_for_radar)
colnames(data_for_radar_t) <- c("EM2", "EM1", "BM2", "BM1") # ensure consistency with original sample names
data_for_radar_t <- as.data.frame(data_for_radar_t)

# Add max/min rows required by radarchart API
gm1 <- rbind(
        max = rep(100, ncol(data_for_radar_t)), # upper bound for percentage scale
        min = rep(0, ncol(data_for_radar_t)),   # lower bound
        data_for_radar_t                        # observed data
)

# Radar chart visualizing DV200 and mitochondrial percentage across samples
radarchart(gm1, axistype = 1,
           pcol = c("#4198FF", "#D02D25"), plwd = 3.5, plty = 1,
           caxislabels = seq(0, 100, by = 20),
           cglcol = "grey", cglty = 1, cglwd = 1.5, calcex = 1.8,
           axislabcol = "grey", vlcex = 2
)

# Legend describing samples (rows beyond max/min in gm1)
legend("topright", legend = rownames(gm1)[-c(1,2)], horiz = TRUE,
       bty = "n", pch = 20, col = c("#4198FF", "#D02D25"),
       text.col = "black", cex = 1.5, pt.cex = 2
)

#################################################################################
########################### Section: Spatial QC Maps (BM2) ######################
#################################################################################

# Spatial QC map for nCount_Spatial (BM2)
SpatialFeaturePlot(BM2, features = "nCount_Spatial", crop = F, pt.size.factor = 1.25)
ggsave("/path/to/figures/FigS1M_nCount.pdf")

# Spatial QC map for nFeature_Spatial (BM2)
SpatialFeaturePlot(BM2, features = "nFeature_Spatial", crop = F, pt.size.factor = 1.25)
ggsave("/path/to/figures/FigS1M_nFeature.pdf")
