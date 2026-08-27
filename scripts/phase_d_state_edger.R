#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(edgeR))

script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork=TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork=FALSE)
results_root <- Sys.getenv("M28_RESULTS_DIR", file.path(project_root, "results"))
out <- file.path(results_root, "02_RPE_STATE_DECOMPOSITION")
tab <- read.delim(file.path(out, "STATE_RAWCOUNT_PSEUDOBULK_COUNTS.tsv.gz"), check.names=FALSE, stringsAsFactors=FALSE)
meta <- read.delim(file.path(out, "STATE_PSEUDOBULK_SAMPLE_METADATA.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(identical(colnames(tab)[-(1:2)], meta$sample))
all_counts <- as.matrix(tab[, -(1:2), drop=FALSE])
storage.mode(all_counts) <- "integer"
rownames(all_counts) <- tab$gene_id

all_results <- list()
for (threshold in c(20, 10)) {
  for (state_name in unique(meta$state)) {
    state_meta <- meta[meta$state == state_name, ]
    wide <- reshape(state_meta[, c("donor_id", "Region", "n_cells")], idvar="donor_id", timevar="Region", direction="wide")
    eligible <- wide$donor_id[wide$n_cells.Macular >= threshold & wide$n_cells.Peripheral >= threshold]
    if (length(eligible) < 8) next
    use <- meta$state == state_name & meta$donor_id %in% eligible
    m <- meta[use, ]
    ord <- order(m$donor_id, factor(m$Region, levels=c("Peripheral", "Macular")))
    m <- m[ord, ]
    cts <- all_counts[, match(m$sample, meta$sample), drop=FALSE]
    m$donor_id <- factor(m$donor_id)
    m$Region <- factor(m$Region, levels=c("Peripheral", "Macular"))
    design <- model.matrix(~ donor_id + Region, data=m)
    y <- DGEList(counts=cts, genes=data.frame(gene_id=tab$gene_id, gene=tab$gene, stringsAsFactors=FALSE))
    keep <- filterByExpr(y, design=design)
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y, method="TMM")
    y <- estimateDisp(y, design=design, robust=TRUE)
    fit <- glmQLFit(y, design=design, robust=TRUE)
    qlf <- glmQLFTest(fit, coef="RegionMacular")
    res <- topTags(qlf, n=Inf, sort.by="none")$table
    res$gene_id <- rownames(res)
    res$state <- state_name
    res$threshold_cells <- threshold
    res$paired_donor_n <- length(eligible)
    res <- res[, c("state", "threshold_cells", "paired_donor_n", "gene_id", "gene", "logFC", "logCPM", "F", "PValue", "FDR")]
    all_results[[paste(state_name, threshold, sep="_")]] <- res
  }
}
combined <- do.call(rbind, all_results)
combined_connection <- gzfile(file.path(out, "STATE_RAWCOUNT_EDGER_ALL_GENES.tsv.gz"), open="wt")
write.table(combined, combined_connection, sep="\t", quote=FALSE, row.names=FALSE)
close(combined_connection)
capture.output(sessionInfo(), file=file.path(out, "STATE_EDGER_SESSION_INFO.txt"))
