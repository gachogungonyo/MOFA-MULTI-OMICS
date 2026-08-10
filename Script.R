library(reticulate)
library(ggplot2)
py_require("cptac")
cptac <- import("cptac")
luad <- cptac$Luad()
luad$list_data_sources()

#Clinical
clin_df <- py_to_r(luad$get_clinical(source = "mssm"))
str(clin_df)
colnames(clin_df)
stage_df <- data.frame(sample = rownames(clin_df), stage = clin_df$tumor_stage_pathological)

#CNV
cnv_df <- py_to_r(luad$get_CNV(source = "washu"))
cnv_df <- as.data.frame(apply(cnv_df, 2, function(x) { x[is.nan(x)] <- mean(x, na.rm = TRUE); x }))
cnv_df_v <- cnv_df[, apply(cnv_df, 2, var, na.rm = TRUE) > 0]
cnv_pca <- prcomp(cnv_df_v, scale. = TRUE, center = TRUE)
cnv_pca_df <- data.frame(PC1 = cnv_pca$x[, 1], PC2 = cnv_pca$x[, 2], sample = rownames(cnv_pca$x))
cnv_pca_df <- merge(cnv_pca_df, stage_df, by = "sample")

#circular_RNA
circrna_df <- py_to_r(luad$get_circular_RNA(source = "bcm"))
circrna_df <- as.data.frame(apply(circrna_df, 2, function(x) { x[is.nan(x)] <- mean(x, na.rm = TRUE); x }))
circrna_df_v <- circrna_df[, apply(circrna_df, 2, var, na.rm = TRUE) > 0]
circrna_pca <- prcomp(circrna_df_v, scale. = TRUE, center = TRUE)
circrna_pca_df <- data.frame(PC1 = circrna_pca$x[, 1], PC2 = circrna_pca$x[, 2], sample = rownames(circrna_pca$x))
circrna_pca_df <- merge(circrna_pca_df, stage_df, by = "sample")

#miRNA
mirna_df <- py_to_r(luad$get_miRNA(source = "washu"))
mirna_df <- as.data.frame(apply(mirna_df, 2, function(x) { x[is.nan(x)] <- mean(x, na.rm = TRUE); x }))
mirna_df_v <- mirna_df[, apply(mirna_df, 2, var, na.rm = TRUE) > 0]
mirna_pca <- prcomp(mirna_df_v, scale. = TRUE, center = TRUE)
mirna_pca_df <- data.frame(PC1 = mirna_pca$x[, 1], PC2 = mirna_pca$x[, 2], sample = rownames(mirna_pca$x))
mirna_pca_df <- merge(mirna_pca_df, stage_df, by = "sample")

#Proteomics
prot_df <- py_to_r(luad$get_proteomics(source = "umich"))
prot_df <- as.data.frame(apply(prot_df, 2, function(x) { x[is.nan(x)] <- mean(x, na.rm = TRUE); x }))
prot_df_v <- prot_df[, apply(prot_df, 2, var, na.rm = TRUE) > 0]
prot_pca <- prcomp(prot_df_v, scale. = TRUE, center = TRUE)
prot_pca_df <- data.frame(PC1 = prot_pca$x[, 1], PC2 = prot_pca$x[, 2], sample = rownames(prot_pca$x))
prot_pca_df <- merge(prot_pca_df, stage_df, by = "sample")

#Phosphoproteomics
phospho_df <- py_to_r(luad$get_phosphoproteomics(source = "umich"))
phospho_df <- as.data.frame(apply(phospho_df, 2, function(x) { x[is.nan(x)] <- mean(x, na.rm = TRUE); x }))
phospho_df_v <- phospho_df[, apply(phospho_df, 2, var, na.rm = TRUE) > 0]
phospho_pca <- prcomp(phospho_df_v, scale. = TRUE, center = TRUE)
phospho_pca_df <- data.frame(PC1 = phospho_pca$x[, 1], PC2 = phospho_pca$x[, 2], sample = rownames(phospho_pca$x))
phospho_pca_df <- merge(phospho_pca_df, stage_df, by = "sample")

#Acetylproteomics
acetylprot_df <- py_to_r(luad$get_acetylproteomics(source = "umich"))
acetylprot_df <- as.data.frame(apply(acetylprot_df, 2, function(x) { x[is.nan(x)] <- mean(x, na.rm = TRUE); x }))
acetylprot_df_v <- acetylprot_df[, apply(acetylprot_df, 2, var, na.rm = TRUE) > 0]
acetylprot_pca <- prcomp(acetylprot_df_v, scale. = TRUE, center = TRUE)
acetylprot_pca_df <- data.frame(PC1 = acetylprot_pca$x[, 1], PC2 = acetylprot_pca$x[, 2], sample = rownames(acetylprot_pca$x))
acetylprot_pca_df <- merge(acetylprot_pca_df, stage_df, by = "sample")

#Transcriptomics
rna_df <- py_to_r(luad$get_transcriptomics(source = "broad"))
rna_df <- as.data.frame(apply(rna_df, 2, function(x) { x[is.nan(x)] <- mean(x, na.rm = TRUE); x }))
rna_df_v <- rna_df[, apply(rna_df, 2, var, na.rm = TRUE) > 0]
rna_pca <- prcomp(rna_df_v, scale. = TRUE, center = TRUE)
rna_pca_df <- data.frame(PC1 = rna_pca$x[, 1], PC2 = rna_pca$x[, 2], sample = rownames(rna_pca$x))
rna_pca_df <- merge(rna_pca_df, stage_df, by = "sample")

#Combined histogram grid
dev.off()
par(mfrow = c(3, 3), mar = c(3, 3, 2, 1), cex.main = 0.9)
hist(unlist(cnv_df_v), main = "CNV", xlab = "Value", col = "blue", breaks = 50)
hist(unlist(circrna_df_v), main = "circRNA", xlab = "Value", col = "green", breaks = 50)
hist(unlist(mirna_df_v), main = "miRNA", xlab = "Value", col = "yellow", breaks = 50)
hist(unlist(prot_df_v), main = "Proteomics", xlab = "Value", col = "purple", breaks = 50)
hist(unlist(phospho_df_v), main = "Phosphoprot.", xlab = "Value", col = "brown", breaks = 50)
hist(unlist(acetylprot_df_v), main = "Acetylprot.", xlab = "Value", col = "pink", breaks = 50)
hist(unlist(rna_df_v), main = "Transcriptomics", xlab = "Value", col = "red", breaks = 50)
par(mfrow = c(1, 1))

#PCA plots
pca_theme <- theme(
  plot.title = element_text(size = 11, hjust = 0.5),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8),
  axis.title = element_text(size = 9)
)

ggplot(cnv_pca_df, aes(PC1, PC2, color = stage)) + geom_point(size = 2) + ggtitle("CNV PCA") + pca_theme
ggplot(circrna_pca_df, aes(PC1, PC2, color = stage)) + geom_point(size = 2) + ggtitle("circRNA PCA") + pca_theme
ggplot(mirna_pca_df, aes(PC1, PC2, color = stage)) + geom_point(size = 2) + ggtitle("miRNA PCA") + pca_theme
ggplot(prot_pca_df, aes(PC1, PC2, color = stage)) + geom_point(size = 2) + ggtitle("Proteomics PCA") + pca_theme
ggplot(phospho_pca_df, aes(PC1, PC2, color = stage)) + geom_point(size = 2) + ggtitle("Phosphoproteomics PCA") + pca_theme
ggplot(acetylprot_pca_df, aes(PC1, PC2, color = stage)) + geom_point(size = 2) + ggtitle("Acetylproteomics PCA") + pca_theme
ggplot(rna_pca_df, aes(PC1, PC2, color = stage)) + geom_point(size = 2) + ggtitle("Transcriptomics PCA") + pca_theme

#MOFA
library(MOFA2)

clean_names <- function(x) sapply(strsplit(gsub("[()']", "", x), ", "), `[`, 1)

top_var_features <- function(m, n = 5000) {
  v <- apply(m, 1, var, na.rm = TRUE)
  m[order(v, decreasing = TRUE)[1:min(n, nrow(m))], ]
}

mofa_list <- list(
  CNV = t(as.matrix(cnv_df_v)),
  circRNA = t(as.matrix(circrna_df_v)),
  miRNA = t(as.matrix(mirna_df_v)),
  Proteomics = t(as.matrix(prot_df_v)),
  Phosphoproteomics = t(as.matrix(phospho_df_v)),
  Acetylproteomics = t(as.matrix(acetylprot_df_v)),
  Transcriptomics = t(as.matrix(rna_df_v))
)

sample_lists <- lapply(mofa_list, colnames)
sample_lists$stage <- stage_df$sample
common_samples_all <- Reduce(intersect, sample_lists)
length(common_samples_all)

mofa_list_all <- lapply(mofa_list, function(m) m[, common_samples_all])
mofa_list_all <- lapply(mofa_list_all, top_var_features, n = 5000)
for (v in names(mofa_list_all)) rownames(mofa_list_all[[v]]) <- clean_names(rownames(mofa_list_all[[v]]))
sapply(mofa_list_all, ncol)

MOFAobject_all <- create_mofa(mofa_list_all)
plot_data_overview(MOFAobject_all)

data_opts <- get_default_data_options(MOFAobject_all)
model_opts <- get_default_model_options(MOFAobject_all)
model_opts$num_factors <- 15
train_opts <- get_default_training_options(MOFAobject_all)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42

MOFAobject_all <- prepare_mofa(
  object = MOFAobject_all,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile_all <- file.path(getwd(), "luad_mofa_all_views.hdf5")
MOFAobject_all.trained <- run_mofa(MOFAobject_all, outfile_all, use_basilisk = TRUE)

p1 <- plot_variance_explained(MOFAobject_all.trained, x = "view", y = "factor")
p1 + theme(axis.text.x = element_text(size = 8, angle = 40, hjust = 1), axis.text.y = element_text(size = 8))

p2 <- plot_variance_explained(MOFAobject_all.trained, x = "view", y = "factor", plot_total = TRUE)[[2]]
p2 + theme(axis.text.x = element_text(size = 8, angle = 40, hjust = 1))

#Clinical covariates
colnames(clin_df)

clinical_all <- data.frame(
  sample = rownames(clin_df),
  stage = clin_df$tumor_stage_pathological,
  age = clin_df$age,
  histology = clin_df$histologic_type,
  smoking = clin_df$tobacco_smoking_history,
  tumor_size = clin_df$tumor_size
)

samples_metadata(MOFAobject_all.trained) <- merge(
  samples_metadata(MOFAobject_all.trained),
  clinical_all,
  by = "sample", all.x = TRUE
)

for (v in c("stage", "age", "histology", "smoking", "tumor_size")) {
  print(plot_factor(MOFAobject_all.trained, factor = 1:4, color_by = v) + pca_theme)
}

#Correlation heatmaps
mofa_purple <- colorRampPalette(c("white", "#7B68C4", "#3F1F78"))(100)

correlate_factors_with_covariates(
  MOFAobject_all.trained, covariates = "stage", plot = "log_pval",
  color = mofa_purple, fontsize_row = 9, fontsize_col = 9, angle_col = 45
)

correlate_factors_with_covariates(
  MOFAobject_all.trained, covariates = "age", plot = "log_pval",
  color = mofa_purple, fontsize_row = 9, fontsize_col = 9, angle_col = 45
)

correlate_factors_with_covariates(
  MOFAobject_all.trained, covariates = "histology", plot = "log_pval",
  color = mofa_purple, fontsize_row = 9, fontsize_col = 9, angle_col = 45
)

correlate_factors_with_covariates(
  MOFAobject_all.trained, covariates = "smoking", plot = "log_pval",
  color = mofa_purple, fontsize_row = 9, fontsize_col = 9, angle_col = 45
)

correlate_factors_with_covariates(
  MOFAobject_all.trained, covariates = "tumor_size", plot = "log_pval",
  color = mofa_purple, fontsize_row = 9, fontsize_col = 9, angle_col = 45
)

correlate_factors_with_covariates(
  MOFAobject_all.trained,
  covariates = c("stage", "age", "histology", "smoking", "tumor_size"),
  plot = "log_pval",
  color = mofa_purple, fontsize_row = 9, fontsize_col = 9, angle_col = 45
)

#Weight plots
views_all <- names(mofa_list_all)
for (v in views_all) {
  print(plot_weights(MOFAobject_all.trained, view = v, factor = 1, nfeatures = 10))
  print(plot_top_weights(MOFAobject_all.trained, view = v, factor = 1, nfeatures = 10))
}