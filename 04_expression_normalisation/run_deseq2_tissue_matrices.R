# ============================================================
# DESeq2 normalisation and VST for tissue-specific
# Swiss-Prot expression matrices across species
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

input_directory <- "tissue_swissprot_count_matrices"

# Primary comparative analysis:
analysis_set <- "all_annotated"

# To run the broader annotated-feature sensitivity analysis later,
# change the line above to:
# analysis_set <- "all_annotated"
#analysis_set <- "shared_all_species"

matrix_suffix <- paste0("_swissprot_counts_", analysis_set, ".csv")
output_directory <- paste0("deseq2_tissue_outputs_", analysis_set)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

tissues <- c(
  "brain",
  "embryo",
  "ovary",
  "testes",
  "trophotaenia"
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

read_count_matrix <- function(tissue) {

  matrix_file <- file.path(
    input_directory,
    paste0(tissue, matrix_suffix)
  )

  if (!file.exists(matrix_file)) {
    stop("Count matrix not found: ", matrix_file)
  }

  count_data <- read.csv(
    matrix_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  count_data <- as.matrix(count_data)
  mode(count_data) <- "numeric"

  if (any(!is.finite(count_data))) {
    stop("Non-finite count values found in: ", matrix_file)
  }

  if (any(count_data < 0)) {
    stop("Negative count values found in: ", matrix_file)
  }

  # Salmon NumReads values are fractional estimates.
  # DESeqDataSetFromMatrix requires integer-valued counts.
  count_data <- round(count_data)
  storage.mode(count_data) <- "integer"

  return(count_data)
}


read_metadata <- function(tissue) {

  metadata_file <- file.path(
    input_directory,
    paste0(tissue, "_sample_metadata.tsv")
  )

  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }

  metadata <- read.delim(
    metadata_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "species",
    "sample",
    "sample_code",
    "tissue",
    "sex",
    "quant_file"
  )

  missing_columns <- setdiff(required_columns, colnames(metadata))

  if (length(missing_columns) > 0) {
    stop(
      "Missing metadata columns for ", tissue, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  rownames(metadata) <- metadata$sample

  metadata$species <- factor(metadata$species)
  metadata$sex <- factor(metadata$sex)
  metadata$tissue <- factor(metadata$tissue)

  return(metadata)
}


choose_design <- function(tissue, metadata) {

  metadata$species <- droplevels(metadata$species)
  metadata$sex <- droplevels(metadata$sex)

  if (tissue == "brain") {
    # Brain includes both female brain (B) and male brain (Br)
    return(~ species + sex)
  } else {
    # Within ovary/testes/embryo/trophotaenia, sex is not
    # an informative comparison variable.
    return(~ species)
  }
}


write_matrix <- function(matrix_object, output_file) {

  write.table(
    matrix_object,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )
}


# ------------------------------------------------------------
# Containers for outputs
# ------------------------------------------------------------

dds_list <- list()
normalized_counts_list <- list()
vst_counts_list <- list()
results_names_list <- list()
summary_rows <- list()

# ------------------------------------------------------------
# Process each tissue
# ------------------------------------------------------------

for (tissue in tissues) {

  cat("\n========================================\n")
  cat("Processing tissue:", tissue, "\n")
  cat("========================================\n")

  count_data <- read_count_matrix(tissue)
  metadata <- read_metadata(tissue)

  # ----------------------------------------------------------
  # Match matrix columns and metadata rows
  # ----------------------------------------------------------

  missing_from_metadata <- setdiff(colnames(count_data), rownames(metadata))
  missing_from_matrix <- setdiff(rownames(metadata), colnames(count_data))

  if (length(missing_from_metadata) > 0) {
    stop(
      "Samples in count matrix but absent from metadata for ", tissue, ": ",
      paste(missing_from_metadata, collapse = ", ")
    )
  }

  if (length(missing_from_matrix) > 0) {
    stop(
      "Samples in metadata but absent from count matrix for ", tissue, ": ",
      paste(missing_from_matrix, collapse = ", ")
    )
  }

  # Reorder metadata rows to exactly match count matrix columns
  metadata <- metadata[colnames(count_data), , drop = FALSE]

  stopifnot(identical(colnames(count_data), rownames(metadata)))

  # ----------------------------------------------------------
  # Minimal filtering
  # ----------------------------------------------------------

  original_features <- nrow(count_data)

  # Remove Swiss-Prot features with zero counts in every sample.
  # DESeq2 also applies independent filtering during results testing.
  keep <- rowSums(count_data) > 0
  count_data <- count_data[keep, , drop = FALSE]

  retained_features <- nrow(count_data)

  cat("Samples:", ncol(count_data), "\n")
  cat("Features before zero-count filtering:", original_features, "\n")
  cat("Features retained:", retained_features, "\n")

  # ----------------------------------------------------------
  # Design formula
  # ----------------------------------------------------------

  design_formula <- choose_design(tissue, metadata)

  cat("Design:", paste(deparse(design_formula), collapse = " "), "\n")

  # Confirm that the model matrix is full rank
  design_matrix <- model.matrix(design_formula, data = metadata)

  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop(
      "The DESeq2 design is not full rank for tissue: ", tissue,
      ". Check species/sex representation in metadata."
    )
  }

  # ----------------------------------------------------------
  # Create and run DESeq2 object
  # ----------------------------------------------------------

  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = metadata,
    design = design_formula
  )

  dds <- DESeq(dds)

  dds_list[[tissue]] <- dds
  results_names_list[[tissue]] <- resultsNames(dds)

  # ----------------------------------------------------------
  # Normalised counts
  # ----------------------------------------------------------

  normalized_counts <- counts(dds, normalized = TRUE)
  normalized_counts_list[[tissue]] <- normalized_counts

  normalized_file <- file.path(
    output_directory,
    paste0(tissue, "_deseq2_normalized_counts.tsv")
  )

  write_matrix(normalized_counts, normalized_file)

  # ----------------------------------------------------------
  # Variance-stabilising transformation
  # ----------------------------------------------------------

  # blind = FALSE preserves the design-based dispersion trend.
  # This is appropriate for PCA/clustering when strong biological
  # species or sex effects are expected.
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  vst_counts <- assay(vsd)

  vst_counts_list[[tissue]] <- vst_counts

  vst_file <- file.path(
    output_directory,
    paste0(tissue, "_vst_counts.tsv")
  )

  write_matrix(vst_counts, vst_file)

  # ----------------------------------------------------------
  # Save sample metadata used for this analysis
  # ----------------------------------------------------------

  metadata_file <- file.path(
    output_directory,
    paste0(tissue, "_metadata_used.tsv")
  )

  write.table(
    as.data.frame(colData(dds)),
    file = metadata_file,
    sep = "\t",
    quote = FALSE,
    row.names = TRUE,
    col.names = NA
  )

  # Save individual DESeq2 object
  saveRDS(
    dds,
    file = file.path(output_directory, paste0(tissue, "_dds.rds"))
  )

  # ----------------------------------------------------------
  # Summary information
  # ----------------------------------------------------------

  summary_rows[[tissue]] <- data.frame(
    tissue = tissue,
    matrix_type = analysis_set,
    number_of_features_before_filtering = original_features,
    number_of_features_after_filtering = retained_features,
    number_of_samples = ncol(count_data),
    number_of_species = length(unique(metadata$species)),
    design = paste(deparse(design_formula), collapse = " "),
    stringsAsFactors = FALSE
  )

  cat("Saved normalised counts:", normalized_file, "\n")
  cat("Saved VST counts:", vst_file, "\n")
}

# ------------------------------------------------------------
# Combine and save summary objects
# ------------------------------------------------------------

summary_table <- do.call(rbind, summary_rows)

write.table(
  summary_table,
  file = file.path(output_directory, "deseq2_tissue_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  dds_list,
  file = file.path(output_directory, "all_tissue_dds_objects.rds")
)

saveRDS(
  normalized_counts_list,
  file = file.path(output_directory, "all_tissue_normalized_counts.rds")
)

saveRDS(
  vst_counts_list,
  file = file.path(output_directory, "all_tissue_vst_counts.rds")
)

capture.output(
  results_names_list,
  file = file.path(output_directory, "deseq2_results_names.txt")
)

# ------------------------------------------------------------
# Display final summary
# ------------------------------------------------------------

cat("\n\n========================================\n")
cat("DESeq2 processing complete\n")
cat("========================================\n\n")

print(summary_table)

cat("\nResults names available for contrasts:\n")
print(results_names_list)

cat("\nOutput directory:\n")
cat(output_directory, "\n")
