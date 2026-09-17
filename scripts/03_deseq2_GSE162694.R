suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
})
register(MulticoreParam(16))

base <- "/home/seyoung/academic_symposium/fibrosis"
out_dir <- file.path(base, "results/deseq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(file.path(base, "geo_data/GSE162694/GSE162694_samples.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(samples) == 143)

counts_raw <- read.csv(file.path(base, "geo_data/GSE162694/GSE162694_raw_counts.csv.gz"),
                        row.names = 1, check.names = FALSE)
stopifnot(!any(grepl("^__", rownames(counts_raw))))
stopifnot(setequal(colnames(counts_raw), samples$col_id))
counts_raw <- as.matrix(counts_raw)
mode(counts_raw) <- "integer"
counts <- counts_raw[, samples$col_id]
colnames(counts) <- samples$gsm

run_contrast <- function(name, case_idx, ref_idx, n_case, n_ref) {
  stopifnot(sum(case_idx) == n_case, sum(ref_idx) == n_ref, !any(case_idx & ref_idx))
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
  write.table(tab, file.path(out_dir, paste0("GSE162694_", name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(dds, file.path(out_dir, paste0("GSE162694_", name, "_dds.rds")))

  cat(sprintf("%-28s case=%3d ref=%3d | genes tested=%5d | padj<0.05=%5d | DEG up=%4d down=%4d\n",
              name, n_case, n_ref, nrow(tab), sum(tab$padj < 0.05, na.rm = TRUE),
              sum(tab$deg & tab$log2FC_shrunk > 0), sum(tab$deg & tab$log2FC_shrunk < 0)))
}

st <- samples$fibrosis_stage
f4 <- st == "4"
ctrl <- st == "normal"
f0 <- st == ("0")
f01 <- st %in% c("0", "1")
f3f4 <- st %in% c("3", "4")
f0ctrl <- st %in% c("0", "normal")

run_contrast("main_F4_vs_ctrl",    f4,   ctrl,   12, 31)
run_contrast("main_F4_vs_F0",      f4,   f0,     12, 35)
run_contrast("sens1_F4_vs_F0F1",   f4,   f01,    12, 65)
run_contrast("sens2_F3F4_vs_F0F1", f3f4, f01,    20, 65)
run_contrast("sens3_F4_vs_F0ctrl", f4,   f0ctrl, 12, 66)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "GSE162694_sessionInfo.txt"))
