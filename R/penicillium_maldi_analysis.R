# INTRO----
#SET YOUR WORKING DIRECTORY
setwd("YourWorkingDirectory")

# R PACKAGES
library(readBrukerFlexData)
library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(Rtsne)
library(randomForest)

## NOT AVAILABLE DUE TO DATA ANONYMIZATION (remove "#-" to use)##
#-# DATA IMPORT----
#-spectra <- importBrukerFlex("YourWorkingDirectory", verbose = FALSE)
#-
#-plot(spectra[[1]])
#-
#-spectra <- lapply(spectra, function(x) trim(x, range = c(2500, 4700)))
#-
#-# DATA PROCESSING----
#-# NORMALIZATION
#-spectra_norm <- calibrateIntensity(spectra, method = "TIC")
#-
#-plot(spectra_norm[[1]])
#-
#-# NAME EXTRACTION
#-sample_names <- sapply(spectra_norm, function(x) metaData(x)$sampleName)
#-
#-# AVERAGING
#-spectra_avg <- averageMassSpectra(spectra_norm, labels = sample_names, method = "mean")
#-
#-plot(spectra_avg[[1]], 
#-     main = paste("Average Spectrum:", names(spectra_avg)[1]),
#-     col = "blue")
#-
#-# PEAK DETECTION (SNR = 3)
#-peaks_avg <- detectPeaks(spectra_avg, method = "MAD", SNR = 3, halfWindowSize = 20)
#-
#-plot(spectra_avg[[1]], main = "Peak Detection (SNR=3)", col = "black")
#-points(peaks_avg[[1]], col = "red", pch = 4, lwd = 1.5)
#-
#-peak_masses <- mass(peaks_avg[[1]])
#-peak_masses
#-
#-# ALIGNMENT AND FILTRATION
#-peaks_aligned <- binPeaks(peaks_avg, tolerance = 0.0009)
#-
#-peaks_filtered <- filterPeaks(peaks_aligned, minFrequency = 0.1)
#-
# MATRIX CREATION
#-data_matrix <- intensityMatrix(peaks_filtered)
#-
#-data_matrix[is.na(data_matrix)] <- 0
#-
#-rownames(data_matrix) <- names(spectra_avg)

# DUE TO DATA ANONYMIZATION THE SCRIPT STARTS HERE (after setting work directory and loading the R packages)----
data_matrix <- readRDS("data_matrix.rds")

# CLUSTERING----
distances <- dist(data_matrix, method = "euclidean")

clusters <- hclust(distances, method = "ward.D2")

# PLOT (DENDOGRAM)----
plot(clusters, 
     cex = 0.8, 
     main = "Dendrogram", 
     xlab = "Strains", 
     ylab = "Distance", 
     sub = "")

# OBTAINING SPECIES NAMES
strains <- rownames(data_matrix)
species <- sapply(strsplit(strains, "_"), function(x) paste(x[1], x[2], sep=" "))
species <- as.factor(species)
print(table(species))

# STATISTICAL TESTING OF ALL PEAKS
p_values <- apply(data_matrix, 2, function(column) {
  test <- kruskal.test(column ~ species)
  return(test$p.value)
})

# FALSE POSITIVE CORRECTION
p_adjusted <- p.adjust(p_values, method = "BH")

# SELECTING THE BEST BIOMARKERS----
significant_peaks <- sort(p_adjusted[p_adjusted < 1])
#print(head(significant_peaks, 10))

# BEST BIOMARKER VISUALIZATION
best_peak <- names(significant_peaks)[1]
par(mar = c(12, 4, 4, 2) + 0.1)
boxplot(data_matrix[, best_peak] ~ species, 
        main = paste("Distribution of the best biomarker (m/z:", best_peak, ")"),
        ylab = "Peak intensity", 
        col = "lightblue",
        las = 2,
        cex.axis = 0.8)

# BIOMARKER TABLE
biomarker_table <- data.frame(
  Species = character(),
  Biomarker_1 = character(),
  Biomarker_2 = character(),
  Biomarker_3 = character(),
  stringsAsFactors = FALSE
)
for (current_species in levels(species)) {
  is_current <- species == current_species
  p_values <- apply(data_matrix, 2, function(column) {
    tryCatch({
      wilcox.test(column[is_current], column[!is_current], exact = FALSE)$p.value
    }, error = function(e) 1) 
  })
  mean_current <- colMeans(data_matrix[is_current, , drop = FALSE])
  mean_other <- colMeans(data_matrix[!is_current, , drop = FALSE])
  candidates <- names(p_values)[p_values < 0.05 & mean_current > mean_other]
  sorted_candidates <- candidates[order(p_values[candidates])]
  top_3 <- head(sorted_candidates, 3)
  length(top_3) <- 3 
  biomarker_table[nrow(biomarker_table) + 1, ] <- c(current_species, top_3)
}
print(biomarker_table)

# PCA----
# PC1 X PC2
pca_model <- prcomp(data_matrix, center = TRUE, scale. = FALSE)
pca_data <- data.frame(
  PC1 = pca_model$x[, 1],  
  PC2 = pca_model$x[, 2], 
  Species = species
)

colors <- c("red", "blue", "green3", "orange", "purple", "cyan", "brown", "magenta")

# Calculation of explained variance
variance_explained <- (pca_model$sdev^2) / sum(pca_model$sdev^2) * 100
pc1_var <- round(variance_explained[1], 1)
pc2_var <- round(variance_explained[2], 1)
pc3_var <- round(variance_explained[3], 1)

par(mar = c(5, 4, 4, 12), xpd = TRUE)

plot(pca_data$PC1, pca_data$PC2, 
     col = colors[as.numeric(pca_data$Species)], 
     pch = 19,
     cex = 1.5,
     main = "PC1 x PC2", 
     xlab = paste0("PC1 (", pc1_var, " %)"),    
     ylab = paste0("PC2 (", pc2_var, " %)"))   

legend("topleft", 
       inset = c(1.02, 0), 
       legend = levels(pca_data$Species), 
       col = colors, 
       pch = 19, 
       title = "Species",
       cex = 0.8, 
       text.font = 3, 
       bty = "n")

par(mar = c(5, 4, 4, 2), xpd = FALSE)

# HIGHLIGHTING STRAINS IN PCA PLOT (THIS WAS CHANGED DUE TO DATA ANONYMIZATION)
target_codes <- c("sample1") #Change "sample1" to what you want to highlight

all_names <- rownames(data_matrix)

is_highlighted <- sapply(all_names, function(name) {
  any(sapply(target_codes, function(code) grepl(code, name)))
})

pca_data$Strain <- all_names
pca_data$Highlighted <- is_highlighted

p_pca_highlight <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(data = subset(pca_data, !Highlighted), 
             aes(color = Species), alpha = 0.3, size = 3, shape = 19) +
  geom_point(data = subset(pca_data, Highlighted), 
             aes(color = Species), alpha = 1, size = 3, shape = 19) +
  theme_minimal() +
  scale_color_manual(values = colors) + 
  labs(title = "PCA with highlighted purified/subcultured strains",
       x = "PC1", 
       y = "PC2",
       color = "Species") +
  theme(legend.position = "right")

print(p_pca_highlight)

# PCA DOT IDENTIFICATION (if you want to use this remove "#*")----
#*plot(pca_data$PC1, pca_data$PC2, main="PC1 x PC2", col=colors[as.numeric(pca_data$Species)], pch=19)
#*identify(pca_data$PC1, pca_data$PC2, labels = rownames(pca_data))

# PC3
summary(pca_model)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# PC1 X PC3
plot(pca_model$x[, 1], pca_model$x[, 3], 
     col = colors[as.numeric(species)], 
     pch = 19, cex = 1.5,
     main = "PC1 x PC3",
     xlab = paste0("PC1 (", pc1_var, " %)"), 
     ylab = paste0("PC3 (", pc3_var, " %)"))

# PC2 X PC3
plot(pca_model$x[, 2], pca_model$x[, 3], 
     col = colors[as.numeric(species)], 
     pch = 19, cex = 1.5, 
     main = "PC2 x PC3",
     xlab = paste0("PC2 (", pc2_var, " %)"), 
     ylab = paste0("PC3 (", pc3_var, " %)"))

par(mfrow = c(1, 1))

# BIOMARKER VALIDATION (SENSITIVITY AND SPECIFICITY)----
validation_results <- data.frame(
  Species = character(),
  Biomarker_mz = character(),
  Sensitivity_pct = numeric(),
  Specificity_pct = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(biomarker_table)) {
  examined_species <- biomarker_table$Species[i]
  
  peaks_to_check <- c(biomarker_table$Biomarker_1[i], 
                      biomarker_table$Biomarker_2[i], 
                      biomarker_table$Biomarker_3[i])
  
  for (peak in peaks_to_check) {
    if (!is.na(peak) && peak %in% colnames(data_matrix)) {
      strains_with_peak <- data_matrix[, peak] > 0
      is_our_species <- species == examined_species
      sensitivity <- (sum(strains_with_peak[is_our_species]) / sum(is_our_species)) * 100
      specificity <- (sum(!strains_with_peak[!is_our_species]) / sum(!is_our_species)) * 100
      validation_results[nrow(validation_results) + 1, ] <- c(
        examined_species, 
        peak, 
        round(sensitivity, 1), 
        round(specificity, 1)
      )
    }
  }
}

print(validation_results)


# HEATMAP----
all_biomarkers <- na.omit(unique(c(biomarker_table$Biomarker_1, 
                                   biomarker_table$Biomarker_2, 
                                   biomarker_table$Biomarker_3)))

biomarker_matrix <- data_matrix[, all_biomarkers]

matrix_log <- log2(biomarker_matrix + 1)

row_annotation <- data.frame(Species = species) 
rownames(row_annotation) <- rownames(matrix_log)

pheatmap(matrix_log,
         annotation_row = row_annotation, 
         cluster_cols = FALSE,            
         cluster_rows = TRUE,             
         show_rownames = FALSE, #set to "TRUE" if you want to see the strain names          
         color = colorRampPalette(c("white", "navyblue"))(10), 
         main = "Biomarkers across strains")


# t-SNE----
#set.seed(42) #t-SNE is random, if you set a seed, the results will always be the same

tsne_result <- Rtsne(as.matrix(data_matrix), 
                     dims = 2, 
                     perplexity = 2, #for the original data set this was set to "10"
                     check_duplicates = FALSE)

tsne_data <- data.frame(
  X = tsne_result$Y[, 1], 
  Y = tsne_result$Y[, 2], 
  Species = species
)

plot(tsne_data$X, tsne_data$Y, 
     col = colors[as.numeric(tsne_data$Species)], 
     pch = 19, 
     cex = 1.5, 
     main = "t-SNE Analysis (perplexity = 10)", 
     xlab = "t-SNE 1", 
     ylab = "t-SNE 2")


# DENDOGRAM (2-5 kDa)
mz_values <- as.numeric(colnames(data_matrix))
min_mz <- 2000
max_mz <- 5000

columns_in_range <- mz_values >= min_mz & mz_values <= max_mz

matrix_cropped <- data_matrix[, columns_in_range]

print(paste("Original peak count:", ncol(data_matrix)))
print(paste("Peak count (2-5 kDa):", ncol(matrix_cropped)))

distances_cropped <- dist(matrix_cropped, method = "euclidean")
clusters_cropped <- hclust(distances_cropped, method = "ward.D2")

plot(clusters_cropped, 
     cex = 0.8, 
     main = "Dendrogram (Range 2000 - 5000 Da only)", 
     xlab = "Strains", 
     ylab = "Distance", 
     sub = "")


# MACHINE LEARNING----
#set.seed(42)  #randomForest is random, if you set a seed, the results will always be the same

rf_model <- randomForest(x = data_matrix, 
                         y = species, 
                         ntree = 500, 
                         importance = TRUE)

print(rf_model)


# IDENTIFICATION OF MOST IMPORTANT PEAKS
peak_importance <- importance(rf_model)
sorted_peaks <- peak_importance[order(peak_importance[, "MeanDecreaseAccuracy"], decreasing = TRUE), ]

rounded_names <- sprintf("%.2f", as.numeric(rownames(rf_model$importance)))
rownames(rf_model$importance) <- rounded_names
if(!is.null(rf_model$importanceSD)) {
  rownames(rf_model$importanceSD) <- rounded_names
}

varImpPlot(rf_model, 
           n.var = 15, # Shows best 15
           main = "Peak Importance (Random Forest)",
           color = "darkblue",
           pch = 19)


# NOISE HUNTING (SNR = 2) (NOT AVAILABLE DUE TO DATA ANONYMIZATION - if you want to use this remove "#^")
#^peaks_noise <- detectPeaks(spectra_avg, method = "MAD", SNR = 2, halfWindowSize = 20)
#^peaks_aligned_noise <- binPeaks(peaks_noise, tolerance = 0.0009)
#^peaks_filtered_noise <- filterPeaks(peaks_aligned_noise, minFrequency = 0.03)
#^matrix_noise <- intensityMatrix(peaks_filtered_noise)
#^matrix_noise[is.na(matrix_noise)] <- 0
#^rownames(matrix_noise) <- names(spectra_avg)
#^distances_noise <- dist(matrix_noise, method = "euclidean")
#^clusters_noise <- hclust(distances_noise, method = "ward.D2")
#^
#^par(mar = c(5, 4, 4, 2))
#^plot(clusters_noise, 
#^     cex = 0.6,
#^     main = "Dendrogram (Noise Hunting: SNR = 2, MinFreq = 3%)", 
#^     xlab = "Strains", 
#^     ylab = "Distance (Euclidean)", 
#^     sub = "")

# COMPARISON OF SHARED PEAK INTENSITY (if you want to use this remove "#;")
#;colnames(data_matrix)
#;selected_peak <- "2509.15480534363" #run "colnames(data_matrix)" and choose a peak
#;
#;df_intensity <- data.frame(
#;  Species = species,
#;  Intensity = data_matrix[, selected_peak]
#;)
#;
#;target_group <- c("Penicillium polonicum", "Penicillium cyclopium", "Penicillium melanoconidium", "Penicillium aurantiogriseum", "Penicillium tricolor", "Penicillium viridicatum")
#;df_target <- df_intensity[df_intensity$Species %in% target_group, ]
#;
#;p_box <- ggplot(df_target, aes(x = Species, y = Intensity, fill = Species)) +
#;  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#;  geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black", shape = 21) +
#;  
#;  theme_minimal() +
#;  scale_fill_manual(values = colors) +
#;  labs(title = paste("Intensity variability for the biomarker m/z", round(as.numeric(selected_peak), 2)),
#;       subtitle = "Do the intensities overlap across the target species?",
#;       x = "Species", 
#;       y = "Relative intensity (TIC)") +
#;  theme(legend.position = "none",
#;        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
#;
#;print(p_box)