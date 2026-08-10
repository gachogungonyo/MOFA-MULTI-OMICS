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



#Adding other covariates 
clinical_all <- data.frame(
  sample = rownames(clin_df),
  age = clin_df$age,
  histology = clin_df$histologic_type,
  smoking = clin_df$tobacco_smoking_history,
  tumor_size = clin_df$tumor_size
)

cnv_pca_df <- merge(cnv_pca_df, clinical_all, by = "sample")
circrna_pca_df <- merge(circrna_pca_df, clinical_all, by = "sample")
mirna_pca_df <- merge(mirna_pca_df, clinical_all, by = "sample")
prot_pca_df <- merge(prot_pca_df, clinical_all, by = "sample")
phospho_pca_df <- merge(phospho_pca_df, clinical_all, by = "sample")
acetylprot_pca_df <- merge(acetylprot_pca_df, clinical_all, by = "sample")
rna_pca_df <- merge(rna_pca_df, clinical_all, by = "sample")

#AGE PCA
ggplot(cnv_pca_df, aes(PC1, PC2, color = age)) + geom_point(size = 2) + ggtitle("CNV PCA") + pca_theme
ggplot(circrna_pca_df, aes(PC1, PC2, color = age)) + geom_point(size = 2) + ggtitle("circRNA PCA") + pca_theme
ggplot(mirna_pca_df, aes(PC1, PC2, color = age)) + geom_point(size = 2) + ggtitle("miRNA PCA") + pca_theme
ggplot(prot_pca_df, aes(PC1, PC2, color = age)) + geom_point(size = 2) + ggtitle("Proteomics PCA") + pca_theme
ggplot(phospho_pca_df, aes(PC1, PC2, color = age)) + geom_point(size = 2) + ggtitle("Phosphoproteomics PCA") + pca_theme
ggplot(acetylprot_pca_df, aes(PC1, PC2, color = age)) + geom_point(size = 2) + ggtitle("Acetylproteomics PCA") + pca_theme
ggplot(rna_pca_df, aes(PC1, PC2, color = age)) + geom_point(size = 2) + ggtitle("Transcriptomics PCA") + pca_theme

#Tumor_size PCA
ggplot(cnv_pca_df, aes(PC1, PC2, color = tumor_size)) + geom_point(size = 2) + ggtitle("CNV PCA") + pca_theme
ggplot(circrna_pca_df, aes(PC1, PC2, color = tumor_size)) + geom_point(size = 2) + ggtitle("circRNA PCA") + pca_theme
ggplot(mirna_pca_df, aes(PC1, PC2, color = tumor_size)) + geom_point(size = 2) + ggtitle("miRNA PCA") + pca_theme
ggplot(prot_pca_df, aes(PC1, PC2, color = tumor_size)) + geom_point(size = 2) + ggtitle("Proteomics PCA") + pca_theme
ggplot(phospho_pca_df, aes(PC1, PC2, color = tumor_size)) + geom_point(size = 2) + ggtitle("Phosphoproteomics PCA") + pca_theme
ggplot(acetylprot_pca_df, aes(PC1, PC2, color = tumor_size)) + geom_point(size = 2) + ggtitle("Acetylproteomics PCA") + pca_theme
ggplot(rna_pca_df, aes(PC1, PC2, color = tumor_size)) + geom_point(size = 2) + ggtitle("Transcriptomics PCA") + pca_theme

#Histology PCA
ggplot(cnv_pca_df, aes(PC1, PC2, color = histology)) + geom_point(size = 2) + ggtitle("CNV PCA") + pca_theme
ggplot(circrna_pca_df, aes(PC1, PC2, color = histology)) + geom_point(size = 2) + ggtitle("circRNA PCA") + pca_theme
ggplot(mirna_pca_df, aes(PC1, PC2, color = histology)) + geom_point(size = 2) + ggtitle("miRNA PCA") + pca_theme
ggplot(prot_pca_df, aes(PC1, PC2, color = histology)) + geom_point(size = 2) + ggtitle("Proteomics PCA") + pca_theme
ggplot(phospho_pca_df, aes(PC1, PC2, color = histology)) + geom_point(size = 2) + ggtitle("Phosphoproteomics PCA") + pca_theme
ggplot(acetylprot_pca_df, aes(PC1, PC2, color = histology)) + geom_point(size = 2) + ggtitle("Acetylproteomics PCA") + pca_theme
ggplot(rna_pca_df, aes(PC1, PC2, color = histology)) + geom_point(size = 2) + ggtitle("Transcriptomics PCA") + pca_theme

#Smoking PCA
ggplot(cnv_pca_df, aes(PC1, PC2, color = smoking)) + geom_point(size = 2) + ggtitle("CNV PCA") + pca_theme
ggplot(circrna_pca_df, aes(PC1, PC2, color = smoking)) + geom_point(size = 2) + ggtitle("circRNA PCA") + pca_theme
ggplot(mirna_pca_df, aes(PC1, PC2, color = smoking)) + geom_point(size = 2) + ggtitle("miRNA PCA") + pca_theme
ggplot(prot_pca_df, aes(PC1, PC2, color = smoking)) + geom_point(size = 2) + ggtitle("Proteomics PCA") + pca_theme
ggplot(phospho_pca_df, aes(PC1, PC2, color = smoking)) + geom_point(size = 2) + ggtitle("Phosphoproteomics PCA") + pca_theme
ggplot(acetylprot_pca_df, aes(PC1, PC2, color = smoking)) + geom_point(size = 2) + ggtitle("Acetylproteomics PCA") + pca_theme
ggplot(rna_pca_df, aes(PC1, PC2, color = smoking)) + geom_point(size = 2) + ggtitle("Transcriptomics PCA") + pca_theme

#HEATMAP of all covariates + omics df
clinical_all <- data.frame(
  sample = rownames(clin_df),
  stage = clin_df$tumor_stage_pathological,
  age = clin_df$age,
  histology = clin_df$histologic_type,
  smoking = clin_df$tobacco_smoking_history,
  tumor_size = clin_df$tumor_size
)

cnv_pca_df <- merge(cnv_pca_df[, c("sample", "PC1", "PC2")], clinical_all, by = "sample")
circrna_pca_df <- merge(circrna_pca_df[, c("sample", "PC1", "PC2")], clinical_all, by = "sample")
mirna_pca_df <- merge(mirna_pca_df[, c("sample", "PC1", "PC2")], clinical_all, by = "sample")
prot_pca_df <- merge(prot_pca_df[, c("sample", "PC1", "PC2")], clinical_all, by = "sample")
phospho_pca_df <- merge(phospho_pca_df[, c("sample", "PC1", "PC2")], clinical_all, by = "sample")
acetylprot_pca_df <- merge(acetylprot_pca_df[, c("sample", "PC1", "PC2")], clinical_all, by = "sample")
rna_pca_df <- merge(rna_pca_df[, c("sample", "PC1", "PC2")], clinical_all, by = "sample")

pca_dfs <- list(
  CNV = cnv_pca_df,
  circRNA = circrna_pca_df,
  miRNA = mirna_pca_df,
  Proteomics = prot_pca_df,
  Phosphoproteomics = phospho_pca_df,
  Acetylproteomics = acetylprot_pca_df,
  Transcriptomics = rna_pca_df
)

#creating covariaties object
covariates <- c("stage", "age", "histology", "smoking", "tumor_size")
#removing Nas
sapply(pca_dfs, function(d) sapply(d[covariates], function(x) sum(!is.na(x))))


get_pval <- function(pc_values, covariate_values) {
  if (is.numeric(covariate_values)) {
    test <- cor.test(pc_values, covariate_values)
    return(test$p.value)
  } else {
    x <- as.factor(covariate_values)
    if (nlevels(x) > 1) {
      fit <- aov(pc_values ~ x)
      return(summary(fit)[[1]][["Pr(>F)"]][1])
    } else {
      return(NA)
    }
  }
}

assoc_matrix <- matrix(NA, nrow = length(pca_dfs), ncol = length(covariates),
                       dimnames = list(names(pca_dfs), covariates))

for (view in names(pca_dfs)) {
  df <- pca_dfs[[view]]
  for (cov in covariates) {
    p1 <- get_pval(df$PC1, df[[cov]])
    p2 <- get_pval(df$PC2, df[[cov]])
    assoc_matrix[view, cov] <- max(-log10(p1), -log10(p2), na.rm = TRUE)
  }
}

assoc_matrix

#plotting

library(pheatmap)

pink_gradient <- colorRampPalette(c("white", "#F5A9C5", "#C2185B"))(100)

pheatmap(
  assoc_matrix,
  color = pink_gradient,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  main = "PC1/PC2 Association with Clinical Covariates (max -log10 p)",
  display_numbers = TRUE,
  number_format = "%.1f"
)




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

#variance plot
p1 <- plot_variance_explained(MOFAobject_all.trained, x = "view", y = "factor")
p1 + theme(axis.text.x = element_text(size = 8, angle = 40, hjust = 1), axis.text.y = element_text(size = 8)) +
  scale_fill_gradient(low = "white", high = "#C2185B")

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

#creating matrix
factors <- get_factors(MOFAobject_all.trained, factors = "all")[[1]]
factors <- as.data.frame(factors)
factors$sample <- rownames(factors)

factor_covariates <- merge(factors, clinical_all, by = "sample")

#cleaning
covariates <- c("stage", "age", "histology", "smoking", "tumor_size")
factor_names <- colnames(factors)[colnames(factors) != "sample"]

get_pval <- function(factor_values, covariate_values) {
  if (is.numeric(covariate_values)) {
    p <- cor.test(factor_values, covariate_values)$p.value
  } else {
    x <- as.factor(covariate_values)
    if (nlevels(x) > 1) {
      p <- summary(aov(factor_values ~ x))[[1]][["Pr(>F)"]][1]
    } else {
      p <- NA
    }
  }
  return(p)
}

pval_matrix <- matrix(NA, nrow = length(factor_names), ncol = length(covariates),
                      dimnames = list(factor_names, covariates))

for (f in factor_names) {
  for (cov in covariates) {
    p <- get_pval(factor_covariates[[f]], factor_covariates[[cov]])
    pval_matrix[f, cov] <- -log10(p)
  }
}

pval_matrix

#visualization
library(pheatmap)

pink_gradient <- colorRampPalette(c("white", "#F5A9C5", "#C2185B"))(100)

pheatmap(
  pval_matrix,
  color = pink_gradient,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  main = "Factor-Covariate Associations (-log10 p)",
  display_numbers = TRUE,
  number_format = "%.1f"
)


#Weight plots
views_all <- names(mofa_list_all)
for (v in views_all) {
  print(plot_weights(MOFAobject_all.trained, view = v, factor = 1, nfeatures = 10))
  print(plot_top_weights(MOFAobject_all.trained, view = v, factor = 1, nfeatures = 10))
}



pink_gradient <- colorRampPalette(c("white", "#F5A9C5", "#C2185B"))(100)
pink_shades <- colorRampPalette(c("#FBE4EC", "#F5A9C5", "#E85D9E", "#C2185B", "#8E0E42", "#4A0E28"))

#Correlation network
install.packages("igraph")
library(igraph)

get_top_features <- function(mofa_obj, view, factor = 1, n = 15) {
  w <- get_weights(mofa_obj, views = view, factors = factor)[[1]]
  names(sort(abs(w[,1]), decreasing = TRUE))[1:n]
}

top_feats <- list()
for (v in names(mofa_list_all)) {
  top_feats[[v]] <- get_top_features(MOFAobject_all.trained, v, factor = 1, n = 15)
}

feature_matrix <- do.call(rbind, lapply(names(top_feats), function(v) {
  mofa_list_all[[v]][top_feats[[v]], , drop = FALSE]
}))

cor_mat <- cor(t(feature_matrix), use = "pairwise.complete.obs")
cor_mat[abs(cor_mat) < 0.6] <- 0
diag(cor_mat) <- 0

g <- graph_from_adjacency_matrix(cor_mat, weighted = TRUE, mode = "undirected")
g <- delete_vertices(g, degree(g) == 0)

view_map <- unlist(lapply(names(top_feats), function(v) setNames(rep(v, length(top_feats[[v]])), top_feats[[v]])))
V(g)$view <- view_map[V(g)$name]

n_views <- length(top_feats)
view_colors <- setNames(pink_shades(n_views), names(top_feats))
V(g)$color <- view_colors[V(g)$view]

plot(g, vertex.size = 8, vertex.label.cex = 0.6, edge.width = abs(E(g)$weight) * 3,
     edge.color = ifelse(E(g)$weight > 0, "#C2185B", "#F5A9C5"))
legend("topleft", legend = names(view_colors), col = view_colors, pch = 19, cex = 0.7, bty = "n")

#Alluvial plot
install.packages("ggalluvial")
library(ggalluvial)

alluvial_df <- as.data.frame(table(clinical_all$stage, clinical_all$histology, clinical_all$smoking))
colnames(alluvial_df) <- c("stage", "histology", "smoking", "Freq")
alluvial_df <- alluvial_df[alluvial_df$Freq > 0, ]

n_stage <- nlevels(as.factor(clinical_all$stage))

ggplot(alluvial_df, aes(axis1 = stage, axis2 = histology, axis3 = smoking, y = Freq)) +
  geom_alluvium(aes(fill = stage)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Stage", "Histology", "Smoking"), expand = c(.1, .1)) +
  scale_fill_manual(values = pink_shades(n_stage)) +
  theme_minimal() +
  ggtitle("Sample Distribution Across Clinical Covariates")

#Circos diagram
install.packages("circlize")
library(circlize)

top_cnv <- get_top_features(MOFAobject_all.trained, "CNV", factor = 1, n = 10)
top_rna <- get_top_features(MOFAobject_all.trained, "Transcriptomics", factor = 1, n = 10)

cnv_mat <- mofa_list_all[["CNV"]][top_cnv, , drop = FALSE]
rna_mat <- mofa_list_all[["Transcriptomics"]][top_rna, , drop = FALSE]

common <- intersect(colnames(cnv_mat), colnames(rna_mat))
cross_cor <- cor(t(cnv_mat[, common]), t(rna_mat[, common]), use = "pairwise.complete.obs")

link_df <- as.data.frame(as.table(cross_cor))
colnames(link_df) <- c("cnv_feature", "rna_feature", "correlation")
link_df <- link_df[abs(link_df$correlation) > 0.5, ]

circos.clear()
chordDiagram(
  link_df[, c("cnv_feature", "rna_feature", "correlation")],
  col = ifelse(link_df$correlation > 0, "#C2185B", "#F5A9C5"),
  transparency = 0.3,
  annotationTrack = "grid",
  preAllocateTracks = 1
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
}, bg.border = NA)
title("CNV \u2013 Transcriptomics Feature Correlations")

