suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
})
register(MulticoreParam(16))

base <- "/home/seyoung/academic_symposium/fibrosis"
out_dir <- file.path(base, "results/deseq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(file.path(base, "geo_data/GSE135251/GSE135251_samples.tsv"), stringsAsFactors = FALSE)
samples$exclude <- samples$exclude == "True"
stopifnot(nrow(samples) == 216, sum(samples$exclude) == 2)

read_counts <- function(f) {
  read.delim(gzfile(file.path(base, "geo_data/GSE135251/RAW", f)), header = FALSE,
             col.names = c("gene", "count"), colClasses = c("character", "integer"))
}
tabs <- lapply(samples$file, read_counts)
genes <- tabs[[1]]$gene
stopifnot(all(vapply(tabs, function(t) identical(t$gene, genes), logical(1))),
          !any(grepl("^__", genes)))
counts <- do.call(cbind, lapply(tabs, `[[`, "count"))
dimnames(counts) <- list(genes, samples$gsm)

run_contrast <- function(name, case_idx, ref_idx, n_case, n_ref) {
  stopifnot(sum(case_idx) == n_case, sum(ref_idx) == n_ref,
            !any(case_idx & ref_idx), !any(samples$exclude[case_idx | ref_idx]))
  use <- case_idx | ref_idx
  cd <- data.frame(group = factor(ifelse(case_idx[use], "case", "ref"), levels = c("ref", "case")),
                   row.names = samples$gsm[use])
  dds <- DESeqDataSetFromMatrix(counts[, use], cd, ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= min(n_case, n_ref), ]
  dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)

  res <- results(dds, name = "group_case_vs_ref", alpha = 0.05)
  shr <- lfcShrink(dds, coef = "group_case_vs_ref", type = "apeglm", parallel = TRUE, quiet = TRUE)
  tab <- data.frame(ensembl_id = rownames(res), baseMean = res$baseMean,
                    log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                    log2FC_shrunk = shr$log2FoldChange, stat = res$stat,
                    pvalue = res$pvalue, padj = res$padj)
  tab$deg <- !is.na(tab$padj) & tab$padj < 0.05 & abs(tab$log2FC_shrunk) > 0.5
  tab <- tab[order(tab$pvalue), ]
  write.table(tab, file.path(out_dir, paste0("GSE135251_", name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(dds, file.path(out_dir, paste0("GSE135251_", name, "_dds.rds")))

  cat(sprintf("%-28s case=%3d ref=%3d | genes tested=%5d | padj<0.05=%5d | DEG up=%4d down=%4d\n",
              name, n_case, n_ref, nrow(tab), sum(tab$padj < 0.05, na.rm = TRUE),
              sum(tab$deg & tab$log2FC_shrunk > 0), sum(tab$deg & tab$log2FC_shrunk < 0)))
}

st <- samples$fibrosis_stage
nafld <- samples$group_in_paper != "control"
f0_nafld <- st == 0 & nafld
f01_nafld <- st %in% c(0, 1) & nafld

# Main = F4 vs NAFLD F0 only: pooling healthy controls into the reference makes its purity
# cohort-dependent (8/46 here vs 31/66 in GSE162694). Old main kept as sens3 (plan.md §11, 2026-09-17 3rd).
run_contrast("main_F4_vs_F0",       st == 4,                       f0_nafld,                            14, 38)
run_contrast("sens1_F4_vs_F0F1",    st == 4,                       f01_nafld,                           14, 85)
run_contrast("sens2_F3F4_vs_F0F1",  st %in% c(3, 4),               f01_nafld,                           68, 85)
run_contrast("sens3_F4_vs_F0ctrl",  samples$main_contrast == "F4", samples$main_contrast == "reference", 14, 46)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
