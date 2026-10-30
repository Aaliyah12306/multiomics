suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
cat("Actual R:", R.version.string, "\n")
args <- commandArgs(TRUE)
out <- args[1]
input <- Sys.getenv('INPUT_RDS')
stopifnot(nzchar(input), file.exists(input))
cat('Loading input:', input, '\n')
x <- readRDS(input)
stopifnot(inherits(x, 'Seurat'), 'harmony.umap' %in% Reductions(x))
rescols <- grep('res[._]0\\.1$', colnames(x[[]]), value=TRUE)
if (!length(rescols)) stop('No res0.1 metadata found; do not guess clustering')
cat('res0.1 metadata:', paste(rescols, collapse=', '), '\n')
if (!'RNA' %in% Assays(x)) stop('RNA assay missing; requires review')
DefaultAssay(x) <- 'RNA'
genes <- c('UCHL1','STMN2','PRPH','NDRG4','VIP')
if (!all(genes %in% rownames(x[['RNA']]))) stop('Missing requested gene(s): ',paste(setdiff(genes,rownames(x[['RNA']])),collapse=', '))
# Validate normalized expression availability without recomputing normalization.
v <- FetchData(x, vars=genes, slot='data', clean=FALSE)
if (!all(genes %in% colnames(v)) || any(!is.finite(as.matrix(v)))) stop('Invalid expression values')
write.csv(data.frame(gene=genes,min=sapply(v,min),max=sapply(v,max)),file.path(out,'expression_range.csv'),row.names=FALSE)
rm(v)
# Existing harmony.umap is reused; resolution does not change expression values.
for (g in genes) {
  p <- FeaturePlot(x, features=g, reduction='harmony.umap', slot='data', order=TRUE, raster=TRUE) + ggtitle(paste(g, '| harmony.umap | res0.1 object'))
  ggsave(file.path(out,paste0('FeaturePlot_',g,'.png')),p,width=7,height=6,dpi=300)
  ggsave(file.path(out,paste0('FeaturePlot_',g,'.pdf')),p,width=7,height=6)
}
p <- FeaturePlot(x, features=genes, reduction='harmony.umap',slot='data',order=TRUE,raster=TRUE,ncol=3)
ggsave(file.path(out,'FeaturePlot_neuron_genes.png'),p,width=18,height=10,dpi=250)
ggsave(file.path(out,'FeaturePlot_neuron_genes.pdf'),p,width=18,height=10)
writeLines(c(paste('Cells:',ncol(x)),paste('Assay: RNA'),paste('Resolution columns:',paste(rescols,collapse=', ')),capture.output(sessionInfo())), file.path(out,'validation.txt'))
writeLines('All five genes plotted successfully',file.path(out,'SUCCESS.txt'))
