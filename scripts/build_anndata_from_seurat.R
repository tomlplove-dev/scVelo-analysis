#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript build_anndata_from_seurat.R <seurat_file> <out_dir> [spliced_dir] [unspliced_dir] [cell_prefix]\n")
  cat("  seurat_file: Seurat object file (.rda/.rdata or .qs)\n")
  quit(status = 1)
}

seurat_file <- args[1]
out_dir <- args[2]
has_matrices <- (length(args) >= 4 && args[3] != "none" && args[4] != "none")
spliced_dir <- if (has_matrices) args[3] else NA
unspliced_dir <- if (has_matrices) args[4] else NA
cell_prefix <- ifelse(length(args) >= 5, args[5], "")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(R.utils))

fail <- function(msg) {
  cat("ERROR:", msg, "\n")
  quit(status = 1)
}

load_mbc_from_qs <- function(path) {
  if (!requireNamespace("qs", quietly = TRUE)) {
    fail("qs package not installed. Run: install.packages('qs')")
  }

  obj <- qs::qread(path)
  if (inherits(obj, "Seurat")) {
    return(list(mBC = obj, source = "qread object"))
  }

  if (is.list(obj) && "mBC" %in% names(obj) && inherits(obj[["mBC"]], "Seurat")) {
    return(list(mBC = obj[["mBC"]], source = "qread list$mBC"))
  }

  fail(".qs file must contain a Seurat object (directly or as list element named mBC)")
}

load_mbc_from_rda <- function(path) {
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = load_env)
  if (length(loaded_names) == 0) {
    fail(".rda file is empty")
  }

  has_named_mbc <- "mBC" %in% loaded_names && inherits(get("mBC", envir = load_env), "Seurat")
  if (has_named_mbc) {
    return(list(mBC = get("mBC", envir = load_env), source = "mBC"))
  }

  seurat_candidates <- loaded_names[vapply(
    loaded_names,
    function(x) inherits(get(x, envir = load_env), "Seurat"),
    logical(1)
  )]

  if (length(seurat_candidates) == 0) {
    fail(".rda file does not contain a Seurat object")
  }

  if (length(seurat_candidates) > 1) {
    fail(paste0(
      ".rda has multiple Seurat objects (",
      paste(seurat_candidates, collapse = ", "),
      "). Rename your target object to mBC to avoid stale annotation."
    ))
  }

  return(list(
    mBC = get(seurat_candidates[1], envir = load_env),
    source = seurat_candidates[1]
  ))
}

file_ext <- tolower(tools::file_ext(seurat_file))
if (file_ext == "qs") {
  loaded <- load_mbc_from_qs(seurat_file)
} else if (file_ext %in% c("rda", "rdata")) {
  loaded <- load_mbc_from_rda(seurat_file)
} else {
  fail("Unsupported file extension. Use .rda/.rdata or .qs")
}

# Always normalize to variable name `mBC`.
mBC <- loaded$mBC
cat("Loaded Seurat object as mBC from:", loaded$source, "\n")
cat("Cells:", ncol(mBC), "Genes:", nrow(mBC), "\n")
cat("Meta columns:", paste(colnames(mBC@meta.data), collapse = ", "), "\n")

if (!("orig.ident" %in% colnames(mBC@meta.data))) {
  fail("mBC@meta.data$orig.ident is missing. Please update Seurat metadata and rerun.")
}
if (!("seurat_clusters" %in% colnames(mBC@meta.data))) {
  cat("WARNING: seurat_clusters missing in mBC@meta.data; using Idents(mBC) as fallback\n")
  mBC@meta.data$seurat_clusters <- as.character(Idents(mBC))
}

if (has_matrices) {
  emat <- Read10X(spliced_dir)
  nmat <- Read10X(unspliced_dir)
  if (nchar(cell_prefix) > 0) {
    colnames(emat) <- paste0(cell_prefix, colnames(emat))
    colnames(nmat) <- paste0(cell_prefix, colnames(nmat))
  }

  cell_IDs <- intersect(colnames(emat), rownames(mBC@meta.data))
  gene_IDs <- intersect(rownames(emat), rownames(mBC))
  if (length(cell_IDs) == 0) {
    fail("No overlapping cells between spliced matrix and mBC metadata")
  }
  if (length(gene_IDs) == 0) {
    fail("No overlapping genes between spliced matrix and mBC RNA assay")
  }
  cat("Intersected cells:", length(cell_IDs), "genes:", length(gene_IDs), "\n")

  emat <- emat[gene_IDs, cell_IDs]
  nmat <- nmat[gene_IDs, cell_IDs]
  for (layer in c("spliced", "unspliced")) {
    mat <- if (layer == "spliced") emat else nmat
    d <- file.path(out_dir, layer)
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    Matrix::writeMM(mat, file = file.path(d, "matrix.mtx"))
    write.table(colnames(mat), file = file.path(d, "barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    feat_df <- data.frame(rownames(mat), rownames(mat), "Gene Expression")
    write.table(feat_df, file = file.path(d, "features.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    R.utils::gzip(file.path(d, "matrix.mtx"), destname = file.path(d, "matrix.mtx.gz"), overwrite = TRUE)
    R.utils::gzip(file.path(d, "barcodes.tsv"), destname = file.path(d, "barcodes.tsv.gz"), overwrite = TRUE)
    R.utils::gzip(file.path(d, "features.tsv"), destname = file.path(d, "features.tsv.gz"), overwrite = TRUE)
  }

  dfvar <- data.frame(gene_IDs = gene_IDs)
  rownames(dfvar) <- gene_IDs
  write.table(dfvar, file = file.path(out_dir, "var.tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
} else {
  cell_IDs <- rownames(mBC@meta.data)
  cat("Metadata-only mode. Exporting", length(cell_IDs), "cells\n")
}

# Explicitly pull annotation each run from mBC@meta.data.
orig_ident <- as.character(mBC@meta.data[cell_IDs, "orig.ident", drop = TRUE])
seurat_clusters <- as.character(mBC@meta.data[cell_IDs, "seurat_clusters", drop = TRUE])

orig_labels <- unique(orig_ident)
cat("orig.ident labels from mBC@meta.data$orig.ident:", paste(orig_labels, collapse = ", "), "\n")
write.table(
  orig_labels,
  file = file.path(out_dir, "orig_ident_levels.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

dfobs <- data.frame(mBC@meta.data[cell_IDs, , drop = FALSE], stringsAsFactors = FALSE)
dfobs$orig.ident <- orig_ident
dfobs$seurat_clusters <- seurat_clusters
write.table(dfobs, file = file.path(out_dir, "obs.tsv"), sep = "\t", quote = FALSE, row.names = TRUE)

for (red in c("FItSNE", "umap", "pca")) {
  if (!is.null(mBC@reductions[[red]])) {
    emb <- mBC@reductions[[red]]@cell.embeddings
    shared_cells <- intersect(cell_IDs, rownames(emb))
    if (length(shared_cells) == 0) {
      cat("WARNING:", red, "embedding has no overlapping cells; skipped\n")
      next
    }
    emb <- emb[shared_cells, , drop = FALSE]
    write.table(emb, file = file.path(out_dir, paste0("emb_", red, ".tsv")), sep = "\t", quote = FALSE, row.names = TRUE)
  }
}

cat("Done. Wrote metadata and embeddings to", out_dir, "\n")
