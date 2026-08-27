#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(edgeR)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork=TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork=FALSE)
results_root <- Sys.getenv("M28_RESULTS_DIR", file.path(project_root, "results"))
out <- file.path(results_root, "01_DONOR_ROBUSTNESS")
counts_path <- file.path(out, "RAWCOUNT_PSEUDOBULK_COUNTS.tsv.gz")
meta_path <- file.path(out, "RAWCOUNT_PSEUDOBULK_SAMPLE_METADATA.tsv")

tab <- read.delim(counts_path, check.names=FALSE, stringsAsFactors=FALSE)
meta <- read.delim(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(identical(colnames(tab)[-(1:2)], meta$sample))

cts <- as.matrix(tab[, -(1:2), drop=FALSE])
storage.mode(cts) <- "integer"
rownames(cts) <- tab$gene_id
meta$donor_id <- factor(meta$donor_id)
meta$Region <- factor(meta$Region, levels=c("Peripheral", "Macular"))
design <- model.matrix(~ donor_id + Region, data=meta)
if (!"RegionMacular" %in% colnames(design)) stop("RegionMacular coefficient absent")

y <- DGEList(counts=cts, genes=data.frame(gene_id=tab$gene_id, gene=tab$gene, stringsAsFactors=FALSE))
keep <- filterByExpr(y, design=design)
yf <- y[keep, , keep.lib.sizes=FALSE]
yf <- calcNormFactors(yf, method="TMM")
yf <- estimateDisp(yf, design=design, robust=TRUE)
fit <- glmQLFit(yf, design=design, robust=TRUE)
qlf <- glmQLFTest(fit, coef="RegionMacular")
res <- topTags(qlf, n=Inf, sort.by="none")$table
res$gene_id <- rownames(res)
res <- res[, c("gene_id", "gene", "logFC", "logCPM", "F", "PValue", "FDR")]
res <- res[order(res$FDR, -res$logFC), ]
write.table(res, file.path(out, "RAWCOUNT_PSEUDOBULK_ALL_GENES.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

logcpm <- cpm(yf, log=TRUE, prior.count=0.25)
logcpm_tab <- data.frame(gene_id=rownames(logcpm), gene=yf$genes$gene, logcpm, check.names=FALSE)
logcpm_connection <- gzfile(file.path(out, "RAWCOUNT_TMM_LOGCPM.tsv.gz"), open="wt")
write.table(logcpm_tab, logcpm_connection, sep="\t", quote=FALSE, row.names=FALSE)
close(logcpm_connection)

norm <- data.frame(
  sample=colnames(yf),
  library_size=yf$samples$lib.size,
  norm_factor=yf$samples$norm.factors,
  effective_library_size=yf$samples$lib.size * yf$samples$norm.factors
)
write.table(norm, file.path(out, "RAWCOUNT_TMM_NORMALIZATION_FACTORS.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

capture.output(sessionInfo(), file=file.path(out, "EDGER_SESSION_INFO.txt"))
