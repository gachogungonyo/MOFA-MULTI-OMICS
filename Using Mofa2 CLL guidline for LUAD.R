#  loading libraries and data
library(reticulate)
library(MOFA2)

py_require("cptac")
cptac <- import("cptac")
luad <- cptac$Luad()
luad$list_data_sources()

# clinical data pulling first for covariates
clin_df <- py_to_r(luad$get_clinical(source = "mssm"))
str(clin_df)
colnames(clin_df)

# building stage df with sample column
stage_df <- data.frame(sample = rownames(clin_df),
                       stage = clin_df$tumor_stage_pathological)

# adding cll style covariates
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

# merging into one metadata table
stage_df <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE),
                   list(stage_df, age_df, gender_df, histology_df,
                        smoking_df, tumor_size_df, TTT_df,
                        treatedAfter_df, TTD_df, Died_df))

# cleaning covariate types 
stage_df$age <- as.numeric(stage_df$age)
stage_df$stage <- factor(stage_df$stage,
                         levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
                         ordered = TRUE)
stage_df$gender <- as.factor(stage_df$gender)
stage_df$smoking_history <- as.factor(stage_df$smoking_history)

#collapsing histology into broader groups
stage_df$histology_grouped <- dplyr::case_when(
  grepl("squamous", stage_df$histology, ignore.case = TRUE) &
    grepl("adeno", stage_df$histology, ignore.case = TRUE) ~ "Adenosquamous",
  grepl("squamous", stage_df$histology, ignore.case = TRUE) ~ "Squamous",
  grepl("adenocarcinoma", stage_df$histology, ignore.case = TRUE) ~ "Adenocarcinoma",
  stage_df$histology == "not specified" ~ NA_character_,
  TRUE ~ "Other"
)
stage_df$histology_grouped <- as.factor(stage_df$histology_grouped)
table(stage_df$histology_grouped, useNA = "ifany")

#creating omics dfs

#cnv 
cnv_df <- py_to_r(luad$get_CNV(source = "washu"))
cnv_df <- as.data.frame(apply(cnv_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
cnv_df_v <- cnv_df[, apply(cnv_df, 2, var, na.rm = TRUE) > 0]

#circular rna 
circrna_df <- py_to_r(luad$get_circular_RNA(source = "bcm"))
circrna_df <- as.data.frame(apply(circrna_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
circrna_df_v <- circrna_df[, apply(circrna_df, 2, var, na.rm = TRUE) > 0]

#mirna 
mirna_df <- py_to_r(luad$get_miRNA(source = "washu"))
mirna_df <- as.data.frame(apply(mirna_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
mirna_df_v <- mirna_df[, apply(mirna_df, 2, var, na.rm = TRUE) > 0]

#proteomics 
prot_df <- py_to_r(luad$get_proteomics(source = "umich"))
prot_df <- as.data.frame(apply(prot_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
prot_df_v <- prot_df[, apply(prot_df, 2, var, na.rm = TRUE) > 0]

# phosphoproteomics 
phospho_df <- py_to_r(luad$get_phosphoproteomics(source = "umich"))
phospho_df <- as.data.frame(apply(phospho_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
phospho_df_v <- phospho_df[, apply(phospho_df, 2, var, na.rm = TRUE) > 0]

# acetylproteomics 
acetylprot_df <- py_to_r(luad$get_acetylproteomics(source = "umich"))
acetylprot_df <- as.data.frame(apply(acetylprot_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
acetylprot_df_v <- acetylprot_df[, apply(acetylprot_df, 2, var, na.rm = TRUE) > 0]

#transcriptomics 
rna_df <- py_to_r(luad$get_transcriptomics(source = "broad"))
rna_df <- as.data.frame(apply(rna_df, 2, function(x) {
  x[is.nan(x)] <- mean(x, na.rm = TRUE); x
}))
rna_df_v <- rna_df[, apply(rna_df, 2, var, na.rm = TRUE) > 0]

#error six views include normal adjacent samples
#solution droping rows ending in dot n suffix
drop_normals <- function(df) df[!grepl("\\.N$", rownames(df)), , drop = FALSE]
circrna_df_v <- drop_normals(circrna_df_v)
mirna_df_v <- drop_normals(mirna_df_v)
prot_df_v <- drop_normals(prot_df_v)
phospho_df_v <- drop_normals(phospho_df_v)
acetylprot_df_v <- drop_normals(acetylprot_df_v)
rna_df_v <- drop_normals(rna_df_v)

# clean feature names #
clean_feature_names <- function(x) {
  x <- as.character(x)
  x <- sub("^\\('([^']+)'.*$", "\\1", x)   # tuple strings -> gene symbol
  x <- make.unique(x)
  x
}

for (nm in c("cnv_df", "circrna_df", "mirna_df", "prot_df",
             "phospho_df", "acetylprot_df", "rna_df")) {
  df <- get(nm)
  colnames(df) <- clean_feature_names(colnames(df))
  assign(nm, df)
}

# variance-filtered _v frames 
zero_var_filter <- function(df) df[, apply(df, 2, var, na.rm = TRUE) > 0]

cnv_df_v        <- zero_var_filter(cnv_df)
circrna_df_v    <- zero_var_filter(circrna_df)
mirna_df_v      <- zero_var_filter(mirna_df)
prot_df_v       <- zero_var_filter(prot_df)
phospho_df_v    <- zero_var_filter(phospho_df)
acetylprot_df_v <- zero_var_filter(acetylprot_df)
rna_df_v        <- zero_var_filter(rna_df)

# variance elbow plots per view
library(ggplot2)

omics_list <- list(
  CNV               = cnv_df_v,
  circRNA           = circrna_df_v,
  miRNA             = mirna_df_v,
  Proteomics        = prot_df_v,
  Phosphoproteomics = phospho_df_v,
  Acetylproteomics  = acetylprot_df_v,
  Transcriptomics   = rna_df_v
)

for (omics_name in names(omics_list)) {
  df <- omics_list[[omics_name]]
  vars_sorted <- sort(apply(df, 2, var, na.rm = TRUE), decreasing = TRUE)
  plot_df <- data.frame(rank = seq_along(vars_sorted), variance = vars_sorted)
  
  p <- ggplot(plot_df, aes(x = rank, y = variance)) +
    geom_line() +
    scale_y_log10() +
    labs(title = paste(omics_name, "variance by feature rank"),
         x = "feature rank", y = "variance") +
    theme_bw()
  print(p)
}

# triming to top variable features per view
top_variable <- function(df, n = 7000) {
  vars <- apply(df, 2, var, na.rm = TRUE)
  keep <- order(vars, decreasing = TRUE)[1:min(n, ncol(df))]
  df[, keep]
}

cnv_df_v        <- top_variable(cnv_df_v, 7000)
rna_df_v        <- top_variable(rna_df_v, 7000)
phospho_df_v    <- top_variable(phospho_df_v, 7000)
prot_df_v       <- top_variable(prot_df_v, 7000)
acetylprot_df_v <- top_variable(acetylprot_df_v, 7000)

# 4.2 transpose each view to features by samples
to_feature_by_sample <- function(df) {
  colnames(df) <- make.unique(colnames(df))
  m <- t(as.matrix(df))
  colnames(m) <- rownames(df)
  m
}

LUAD_data <- list(
  CNV               = to_feature_by_sample(cnv_df_v),
  circRNA           = to_feature_by_sample(circrna_df_v),
  miRNA             = to_feature_by_sample(mirna_df_v),
  Proteomics        = to_feature_by_sample(prot_df_v),
  Phosphoproteomics = to_feature_by_sample(phospho_df_v),
  Acetylproteomics  = to_feature_by_sample(acetylprot_df_v),
  Transcriptomics   = to_feature_by_sample(rna_df_v)
)

# correcting col names and aligning samples
all_samples <- Reduce(intersect, lapply(LUAD_data, colnames))

LUAD_data <- lapply(LUAD_data, function(mat) mat[, all_samples, drop = FALSE])

sapply(LUAD_data, ncol)   
sapply(LUAD_data, nrow)   


# 3 create the mofa object
MOFAobject <- create_mofa(LUAD_data)
MOFAobject

# 3.1 plot data overview
plot_data_overview(MOFAobject)


# 3.2.1 define data options
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE

# 3.2.2 define model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15

# 3.2.3 define training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

# 3.3 prepare the mofa object
MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# 3.3 error mofapy2 not detected in python binary
# 3.3 solution run mofa with use basilisk true
MOFAobject <- run_mofa(MOFAobject, outfile = "MOFA2_LUAD.hdf5",
                       use_basilisk = TRUE)
saveRDS(MOFAobject, "MOFA2_Object_LUAD.rds")

#saving files
dir.create("Mofa_dataframes", showWarnings = FALSE)

df_names <- Filter(function(x) is.data.frame(get(x, envir = .GlobalEnv)), ls(envir = .GlobalEnv))
#as rds
for (name in df_names) {
  saveRDS(get(name), file = file.path("Mofa_dataframes", paste0(name, ".rds")))
}

# 4 overview of the trained mofa model

# 4.1 slots
slotNames(MOFAobject)
names(MOFAobject@data)
dim(MOFAobject@data$Transcriptomics$group1)
names(MOFAobject@expectations)
dim(MOFAobject@expectations$Z$group1)
dim(MOFAobject@expectations$W$Transcriptomics)

# 4.3 correlation between factors
plot_factor_cor(MOFAobject)

# 4.4.1 plot variance explained by factor
plot_variance_explained(MOFAobject, max_r2 = 15)

# 4.4.2 plot total variance explained per view
plot_variance_explained(MOFAobject, plot_total = TRUE)[[2]]

# 4.4 identify the most weighted factors
r2 <- get_variance_explained(MOFAobject)$r2_per_factor$group1
total_r2_per_factor <- rowSums(r2)
sort(total_r2_per_factor, decreasing = TRUE)
n_views_active <- rowSums(r2 > 1)
sort(n_views_active, decreasing = TRUE)






# 5 characterisation of factor 1

# 5.1 association analysis using built in function
correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("stage", "age", "gender", "smoking_history",
                                                 "tumor_size", "histology_grouped", "Died"),
                                  plot = "log_pval"
)

# 5.2 plot factor values
plot_factor(MOFAobject,
            factors = 1,
            color_by = "Factor1")
plot_factor(MOFAobject,
            factors = 2,
            color_by = "Factor1")
plot_factor(MOFAobject,
            factors = 3,
            color_by = "Factor1")
plot_factor(MOFAobject,
            factors = 4,
            color_by = "Factor1")



#Identifying features in
r2["Factor1", ]

r2["Factor2", ]  

r2["Factor3", ]   
   
# 5.3.1 plot weights 
#factor 1
plot_data_scatter(MOFAobject,
                  view = "CNV",
                  factor = 1,
                  features = 4,
                  sign = "positive",
                  color_by = "stage") + labs(y = "CNV")
plot_data_heatmap(MOFAobject,
                  view = "CNV",
                  factor = 1,
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row")

plot_data_scatter(MOFAobject,
                  view = "Proteomics",
                  factor = 1,
                  features = 4,
                  sign = "positive",
                  color_by = "stage") + labs(y = "Protein abundance")
plot_data_heatmap(MOFAobject,
                  view = "Proteomics",
                  factor = 1,
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row")

plot_data_scatter(MOFAobject,
                  view = "Phosphoproteomics",
                  factor = 1,
                  features = 4,
                  sign = "positive",
                  color_by = "stage") + labs(y = "Phosphosite abundance")
plot_data_heatmap(MOFAobject,
                  view = "Phosphoproteomics",
                  factor = 1,
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row")
#factor 2
plot_data_scatter(MOFAobject,
                  view = "CNV",
                  factor = 2,
                  features = 4,
                  sign = "positive",
                  color_by = "stage") + labs(y = "CNV")
plot_data_heatmap(MOFAobject,
                  view = "CNV",
                  factor = 2,
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row")
#factor 3
plot_data_scatter(MOFAobject,
                  view = "Acetylproteomics",
                  factor = 3,
                  features = 4,
                  sign = "positive",
                  color_by = "stage") + labs(y = "Acetylsite abundance")
plot_data_heatmap(MOFAobject,
                  view = "Acetylproteomics",
                  factor = 3,
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row")

