library(reticulate)
library(ggplot2)
library(dplyr)

# load cptac data
py_require("cptac")
cptac <- import("cptac")
luad <- cptac$Luad()

# pull clinical data
clin_df <- py_to_r(luad$get_clinical(source = "mssm"))

# build stage df
stage_df <- data.frame(sample = rownames(clin_df),
                       stage = clin_df$tumor_stage_pathological)

# add other covariates
age_df <- data.frame(sample = rownames(clin_df), age = clin_df$age)
gender_df <- data.frame(sample = rownames(clin_df), gender = clin_df$sex)
histology_df <- data.frame(sample = rownames(clin_df),
                           histology = clin_df$histologic_type)
smoking_df <- data.frame(sample = rownames(clin_df),
                         smoking_history = clin_df$tobacco_smoking_history)
tumor_size_df <- data.frame(sample = rownames(clin_df),
                            tumor_size = clin_df$tumor_size_cm)
TTT_df <- data.frame(sample = rownames(clin_df),
                     TTT = clin_df$`Recurrence-free survival, days`)
treatedAfter_df <- data.frame(sample = rownames(clin_df),
                              treatedAfter = clin_df$`Recurrence status (1, yes; 0, no)` == 1)
TTD_df <- data.frame(sample = rownames(clin_df),
                     TTD = clin_df$`Overall survival, days`)
Died_df <- data.frame(sample = rownames(clin_df),
                      Died = clin_df$`Survival status (1, dead; 0, alive)` == 1)

# merge covariates together
stage_df <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE),
                   list(stage_df, age_df, gender_df, histology_df,
                        smoking_df, tumor_size_df, TTT_df,
                        treatedAfter_df, TTD_df, Died_df))

# clean covariate types
stage_df$age <- as.numeric(stage_df$age)
stage_df$stage <- factor(stage_df$stage,
                         levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                         ordered = TRUE)
stage_df$gender <- as.factor(stage_df$gender)
stage_df$smoking_history <- as.factor(stage_df$smoking_history)

# collapse histology groups
stage_df$histology_grouped <- dplyr::case_when(
  grepl("squamous", stage_df$histology, ignore.case = TRUE) &
    grepl("adeno", stage_df$histology, ignore.case = TRUE) ~ "Adenosquamous",
  grepl("squamous", stage_df$histology, ignore.case = TRUE) ~ "Squamous",
  grepl("adenocarcinoma", stage_df$histology, ignore.case = TRUE) ~ "Adenocarcinoma",
  stage_df$histology == "not specified" ~ NA_character_,
  TRUE ~ "Other"
)
stage_df$histology_grouped <- as.factor(stage_df$histology_grouped)

# helper to drop normal samples
drop_normals <- function(df) df[!grepl("\\.N$", rownames(df)), , drop = FALSE]

# pull cnv omics data
cnv_df <- py_to_r(luad$get_CNV(source = "washu"))
cnv_df <- as.data.frame(apply(cnv_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
cnv_df_v <- cnv_df[, apply(cnv_df, 2, var, na.rm = TRUE) > 0]

# pull circular rna data
circrna_df <- py_to_r(luad$get_circular_RNA(source = "bcm"))
circrna_df <- as.data.frame(apply(circrna_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
circrna_df_v <- circrna_df[, apply(circrna_df, 2, var, na.rm = TRUE) > 0]
circrna_df_v <- drop_normals(circrna_df_v)

# pull mirna data
mirna_df <- py_to_r(luad$get_miRNA(source = "washu"))
mirna_df <- as.data.frame(apply(mirna_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
mirna_df_v <- mirna_df[, apply(mirna_df, 2, var, na.rm = TRUE) > 0]
mirna_df_v <- drop_normals(mirna_df_v)

# pull proteomics data
prot_df <- py_to_r(luad$get_proteomics(source = "umich"))
prot_df <- as.data.frame(apply(prot_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
prot_df_v <- prot_df[, apply(prot_df, 2, var, na.rm = TRUE) > 0]
prot_df_v <- drop_normals(prot_df_v)

# pull phosphoproteomics data
phospho_df <- py_to_r(luad$get_phosphoproteomics(source = "umich"))
phospho_df <- as.data.frame(apply(phospho_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
phospho_df_v <- phospho_df[, apply(phospho_df, 2, var, na.rm = TRUE) > 0]
phospho_df_v <- drop_normals(phospho_df_v)

# pull acetylproteomics data
acetylprot_df <- py_to_r(luad$get_acetylproteomics(source = "umich"))
acetylprot_df <- as.data.frame(apply(acetylprot_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
acetylprot_df_v <- acetylprot_df[, apply(acetylprot_df, 2, var, na.rm = TRUE) > 0]
acetylprot_df_v <- drop_normals(acetylprot_df_v)

# pull transcriptomics data
rna_df <- py_to_r(luad$get_transcriptomics(source = "broad"))
rna_df <- as.data.frame(apply(rna_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
rna_df_v <- rna_df[, apply(rna_df, 2, var, na.rm = TRUE) > 0]
rna_df_v <- drop_normals(rna_df_v)

# omics views list
omics_list <- list(
  CNV               = cnv_df_v,
  circRNA           = circrna_df_v,
  miRNA             = mirna_df_v,
  Proteomics        = prot_df_v,
  Phosphoproteomics = phospho_df_v,
  Acetylproteomics  = acetylprot_df_v,
  Transcriptomics   = rna_df_v
)

covariates <- c("stage", "age", "gender", "smoking_history",
                "tumor_size", "histology_grouped", "Died")

# run pca per omics view
pca_results <- list()

for (omics_name in names(omics_list)) {
  df <- omics_list[[omics_name]]
  pca <- prcomp(df, scale. = TRUE, center = TRUE)
  pca_df <- data.frame(
    sample = rownames(df),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  pca_df <- merge(pca_df, stage_df, by = "sample")
  pca_results[[omics_name]] <- pca_df
}

# pca plot per view covariate
for (omics_name in names(pca_results)) {
  pca_df <- pca_results[[omics_name]]
  for (cov in covariates) {
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[cov]])) +
      geom_point(size = 2, alpha = 0.8) +
      labs(
        title = paste(omics_name, "PCA colored by", cov),
        x = "PC1", y = "PC2"
      ) +
      theme_bw()
    print(p)
  }
}

# histograms for numeric covariates
numeric_covariates <- c("age", "tumor_size", "TTT", "TTD")

for (cov in numeric_covariates) {
  p <- ggplot(stage_df, aes(x = .data[[cov]])) +
    geom_histogram(bins = 20, fill = "#8E24AA", color = "white", na.rm = TRUE) +
    labs(title = paste("distribution of", cov), x = cov, y = "count") +
    theme_bw()
  print(p)
}

# bar plots for categorical covariates
categorical_covariates <- c("stage", "gender", "smoking_history",
                            "histology_grouped", "treatedAfter", "Died")

for (cov in categorical_covariates) {
  p <- ggplot(stage_df, aes(x = .data[[cov]])) +
    geom_bar(fill = "#8E24AA", na.rm = TRUE) +
    labs(title = paste("distribution of", cov), x = cov, y = "count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}