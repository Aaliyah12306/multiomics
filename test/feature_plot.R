# WholeTissueData_cluster_0810: five-gene projection and expression proportions
# Edit these parameters to run interactively in RStudio.
input_file <- "INPUT_DIRECTORY/WholeTissueData_cluster_0810.rds"
output_dir <- "OUTPUT_DIRECTORY"
genes <- c("UCHL1", "STMN2", "PRPH", "NDRG4", "VIP")
assay_name <- "RNA"
reduction_name <- "harmony.umap"
resolution_column <- NULL # Auto-select only when exactly one res0.1 column exists.

library(Seurat)
library(ggplot2)
library(patchwork)
stopifnot(as.character(getRversion()) == "4.4.1")
stopifnot(file.exists(input_file), dir.exists(output_dir))
output_names <- c("neuron_5genes_harmony_umap.pdf", "expression_proportions.csv",
                  "validation.txt", "SUCCESS.txt")
stopifnot(!any(file.exists(file.path(output_dir, output_names))))
x <- readRDS(input_file)
stopifnot(inherits(x, "Seurat"), assay_name %in% Assays(x),
          reduction_name %in% Reductions(x))
DefaultAssay(x) <- assay_name
candidates <- grep("res[._]?0\\.1$", colnames(x[[]]), value=TRUE)
if (is.null(resolution_column)) {
  if (length(candidates) != 1L) stop("Specify res0.1 column; candidates: ", paste(candidates, collapse=", "))
  resolution_column <- candidates[1]
}
stopifnot(resolution_column %in% candidates)
Idents(x) <- x[[resolution_column, drop=TRUE]]
cells <- colnames(x)
stopifnot(all(genes %in% rownames(x[[assay_name]])))
data_layers <- Layers(x[[assay_name]], search="^data$")
if (length(data_layers) != 1L) stop("Expected one normalized data layer; inspect layers before proceeding")
v <- FetchData(x, vars=genes, cells=cells, layer="data", clean=FALSE)
if (!setequal(rownames(v), cells)) stop("Expression data does not cover all cells")
v <- v[cells, genes, drop=FALSE]
stopifnot(all(is.finite(as.matrix(v))), all(as.matrix(v) >= 0))
embedding <- Embeddings(x, reduction=reduction_name)
stopifnot(setequal(rownames(embedding), cells), all(is.finite(embedding)))
# Denominator is all cells in the object; expression means normalized RNA data > 0.
stats <- data.frame(gene=genes, n_total=length(cells),
                    n_expressing=as.integer(colSums(v > 0)),
                    percent_expressing=100*colMeans(v > 0),
                    threshold=0, assay=assay_name, layer="data")
write.csv(stats, file.path(output_dir, "expression_proportions.csv"), row.names=FALSE)
p <- FeaturePlot(x, features=genes, cells=cells, reduction=reduction_name,
                 order=TRUE, raster=TRUE, raster.dpi=c(1000,1000),
                 pt.size=0.1, ncol=3) +
  plot_annotation(title="Neuronal markers | harmony.umap | resolution 0.1")
ggsave(file.path(output_dir, "neuron_5genes_harmony_umap.pdf"),
       p, width=18, height=11, units="in", limitsize=TRUE)
writeLines(c(paste("Input:", input_file), paste("Total cells:", length(cells)),
             paste("Resolution column:", resolution_column),
             "No reclustering or cell-type filtering performed.",
             "PDF includes rasterized points; labels are retained in PDF.",
             capture.output(sessionInfo())), file.path(output_dir,"validation.txt"))
writeLines("Analysis and output writes completed", file.path(output_dir,"SUCCESS.txt"))
