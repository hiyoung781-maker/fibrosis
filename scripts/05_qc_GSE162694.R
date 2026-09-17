suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(ggplot2)
})
register(MulticoreParam(16))
set.seed(162694)

# Diagnose why GSE162694 yields ~5x more DEG than GSE135251 for the same contrast.
# Hypothesis: sample block 548nash1119+ (8/12 F4, 0 F0, 0 normal) is a hidden batch confounded with F4.

base <- "/home/seyoung/academic_symposium/fibrosis"
out_dir <- file.path(base, "results/qc/GSE162694")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim(file.path(base, "geo_data/GSE162694/GSE162694_samples.tsv"), stringsAsFactors = FALSE)
mat_lines <- readLines(gzfile(file.path(base, "geo_data/GSE162694/GSE162694_series_matrix.txt.gz")))
nas_line <- grep("^!Sample_characteristics_ch1.*nas score", mat_lines, value = TRUE)
samples$nas <- sub("^nas score: ", "", gsub('"', "", strsplit(nas_line, "\t")[[1]][-1]))
samples$nas_NA <- samples$nas == "NA"
samples$num <- as.integer(sub("548nash", "", samples$col_id))
samples$block <- ifelse(samples$num >= 1119, "1119plus", "other")
samples$stage <- factor(samples$fibrosis_stage, levels = c("normal", "0", "1", "2", "3", "4"))
stopifnot(nrow(samples) == 143, sum(samples$block == "1119plus") == 21)

counts <- as.matrix(read.csv(file.path(base, "geo_data/GSE162694/GSE162694_raw_counts.csv.gz"),
                             row.names = 1, check.names = FALSE))[, samples$col_id]
mode(counts) <- "integer"
colnames(counts) <- samples$gsm

# ---- 1. per-sample library QC ----
samples$lib_size <- colSums(counts)
samples$detected_ge10 <- colSums(counts >= 10)
samples$top100_frac <- apply(counts, 2, function(x) sum(sort(x, decreasing = TRUE)[1:100]) / sum(x))

dds_all <- DESeqDataSetFromMatrix(counts, samples, ~ 1)
dds_all <- estimateSizeFactors(dds_all)
samples$size_factor <- sizeFactors(dds_all)

# ---- 2. PCA on all 143 samples (blind VST, top 500 variable genes) ----
dds_f <- dds_all[rowSums(counts >= 10) >= 12, ]
vsd <- vst(dds_f, blind = TRUE)
rv <- matrixStats::rowVars(assay(vsd))
top <- order(rv, decreasing = TRUE)[1:500]
pca <- prcomp(t(assay(vsd)[top, ]))
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
for (k in 1:5) samples[[paste0("PC", k)]] <- pca$x[, k]

write.table(samples, file.path(out_dir, "sample_qc_pca.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# variance of each PC explained by each covariate (univariate R^2)
covs <- c("stage", "block", "nas_NA", "age", "sex", "lib_size", "detected_ge10", "top100_frac")
r2 <- do.call(rbind, lapply(1:5, function(k) {
  data.frame(PC = paste0("PC", k), pct_var = round(pve[k], 1),
             t(sapply(covs, function(cv) round(summary(lm(samples[[paste0("PC", k)]] ~ samples[[cv]]))$r.squared, 3))))
}))
write.table(r2, file.path(out_dir, "pc_covariate_R2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

theme_set(theme_bw(base_size = 11))
p1 <- ggplot(samples, aes(PC1, PC2, colour = stage, shape = block)) + geom_point(size = 2.4, alpha = 0.85) +
  scale_shape_manual(values = c("1119plus" = 17, other = 16)) +
  labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2]),
       title = "GSE162694 PCA (VST, top 500 variable genes)", subtitle = "triangle = sample block 548nash1119+")
ggsave(file.path(out_dir, "pca_stage_block.png"), p1, width = 7, height = 5.5, dpi = 150)

qc_long <- rbind(
  data.frame(samples[, c("stage", "block")], metric = "library size (M reads)", value = samples$lib_size / 1e6),
  data.frame(samples[, c("stage", "block")], metric = "genes with count >= 10", value = samples$detected_ge10),
  data.frame(samples[, c("stage", "block")], metric = "top-100 gene fraction", value = samples$top100_frac),
  data.frame(samples[, c("stage", "block")], metric = "DESeq2 size factor", value = samples$size_factor))
p2 <- ggplot(qc_long, aes(stage, value, colour = block)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") + labs(title = "GSE162694 library QC by stage and sample block")
ggsave(file.path(out_dir, "library_qc_by_stage_block.png"), p2, width = 10, height = 7, dpi = 150)

# ---- 3. DE tests with the same pipeline as 03_deseq2_GSE162694.R ----
run_de <- function(case_idx, ref_idx) {
  use <- case_idx | ref_idx
  cd <- data.frame(group = factor(ifelse(case_idx[use], "case", "ref"), levels = c("ref", "case")),
                   row.names = samples$gsm[use])
  dds <- DESeqDataSetFromMatrix(counts[, use], cd, ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= min(sum(case_idx), sum(ref_idx)), ]
  dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)
  res <- results(dds, name = "group_case_vs_ref", alpha = 0.05)
  shr <- lfcShrink(dds, coef = "group_case_vs_ref", type = "apeglm", parallel = TRUE, quiet = TRUE)
  deg <- !is.na(res$padj) & res$padj < 0.05 & abs(shr$log2FoldChange) > 0.5
  list(dds = dds, tab = data.frame(ensembl_id = rownames(res), log2FC_shrunk = shr$log2FoldChange,
                                   pvalue = res$pvalue, padj = res$padj, deg = deg),
       summary = c(case = sum(case_idx), ref = sum(ref_idx), tested = nrow(res),
                   padj05 = sum(res$padj < 0.05, na.rm = TRUE),
                   up = sum(deg & shr$log2FoldChange > 0), down = sum(deg & shr$log2FoldChange < 0),
                   median_disp = median(dispersions(dds)[mcols(dds)$baseMean > 10], na.rm = TRUE)))
}

st <- samples$fibrosis_stage
blk <- samples$block == "1119plus"
de_rows <- list()

# (a) same stage, different block: pure block effect
f2_blk <- run_de(st == "2" & blk, st == "2" & !blk)
de_rows[["F2 block1119+ vs F2 other"]] <- f2_blk$summary
write.table(f2_blk$tab, file.path(out_dir, "DE_F2_block1119plus_vs_F2_other.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# (b) null reference: random 10 vs 17 split of non-block F0/F1 samples (same n as test a)
pool <- which(st %in% c("0", "1") & !blk)
for (i in 1:3) {
  pick <- sample(pool, 27)
  case_idx <- seq_along(st) %in% pick[1:10]
  ref_idx <- seq_along(st) %in% pick[11:27]
  de_rows[[sprintf("random split F0/F1 10 vs 17 (seed %d)", i)]] <- run_de(case_idx, ref_idx)$summary
}

# (c) main contrast with and without the block
de_rows[["main F4 vs F0 (all)"]] <- run_de(st == "4", st == "0")$summary
de_rows[["F4 vs F0, F4 outside block only"]] <- run_de(st == "4" & !blk, st == "0")$summary

de_summary <- data.frame(test = names(de_rows), do.call(rbind, de_rows), row.names = NULL, check.names = FALSE)

# dispersion reference from the discovery cohort
d135 <- readRDS(file.path(base, "results/deseq2/GSE135251_main_F4_vs_F0_dds.rds"))
disp135 <- median(dispersions(d135)[mcols(d135)$baseMean > 10], na.rm = TRUE)
write.table(de_summary, file.path(out_dir, "DE_diagnostic_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. does the block effect explain the main-contrast DEG? ----
main <- read.delim(file.path(base, "results/deseq2/GSE162694_main_F4_vs_F0.tsv"))
m <- merge(main[, c("ensembl_id", "log2FC_shrunk", "deg")], f2_blk$tab[, c("ensembl_id", "log2FC_shrunk", "deg")],
           by = "ensembl_id", suffixes = c("_main", "_block"))
overlap <- c(
  main_deg = sum(m$deg_main), block_deg = sum(m$deg_block), both = sum(m$deg_main & m$deg_block),
  both_same_direction = sum(m$deg_main & m$deg_block & sign(m$log2FC_shrunk_main) == sign(m$log2FC_shrunk_block)),
  spearman_lfc_all = round(cor(m$log2FC_shrunk_main, m$log2FC_shrunk_block, method = "spearman"), 3),
  spearman_lfc_main_deg = round(cor(m$log2FC_shrunk_main[m$deg_main], m$log2FC_shrunk_block[m$deg_main], method = "spearman"), 3))
write.table(data.frame(metric = names(overlap), value = overlap), file.path(out_dir, "main_vs_block_effect_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p3 <- ggplot(m, aes(log2FC_shrunk_block, log2FC_shrunk_main, colour = deg_main)) +
  geom_point(size = 0.6, alpha = 0.4) + geom_abline(linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#e34948"), name = "main DEG") +
  labs(x = "log2FC  F2 block1119+ vs F2 other  (block effect, same stage)",
       y = "log2FC  main F4 vs F0",
       title = sprintf("Block effect vs main contrast (Spearman rho all genes = %.2f)", overlap["spearman_lfc_all"]))
ggsave(file.path(out_dir, "block_effect_vs_main_lfc.png"), p3, width = 7, height = 6, dpi = 150)

# ---- 5. NAS-missing group (PC1 R^2 = 0.64) as the technical axis ----
na <- samples$nas_NA
nas_rows <- list()
f0_na <- run_de(st == "0" & na, st == "0" & !na)
nas_rows[["F0 NAS-missing vs F0 NAS-present"]] <- f0_na$summary
write.table(f0_na$tab, file.path(out_dir, "DE_F0_NASmissing_vs_F0_NASpresent.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
f1_na <- run_de(st == "1" & na, st == "1" & !na)
nas_rows[["F1 NAS-missing vs F1 NAS-present"]] <- f1_na$summary
nas_rows[["F4 vs F0, NAS-present only"]] <- run_de(st == "4" & !na, st == "0" & !na)$summary
nas_rows[["F4 vs F0, NAS-missing only"]] <- run_de(st == "4" & na, st == "0" & na)$summary

# diagnostic only (not adopted): same main contrast with NAS-missing as covariate
use <- st %in% c("0", "4")
cd <- data.frame(nas_missing = factor(na[use]), group = factor(ifelse(st[use] == "4", "case", "ref"), levels = c("ref", "case")),
                 row.names = samples$gsm[use])
dds_adj <- DESeqDataSetFromMatrix(counts[, use], cd, ~ nas_missing + group)
dds_adj <- dds_adj[rowSums(counts(dds_adj) >= 10) >= 12, ]
dds_adj <- DESeq(dds_adj, parallel = TRUE, quiet = TRUE)
res_adj <- results(dds_adj, name = "group_case_vs_ref", alpha = 0.05)
shr_adj <- lfcShrink(dds_adj, coef = "group_case_vs_ref", type = "apeglm", parallel = TRUE, quiet = TRUE)
deg_adj <- !is.na(res_adj$padj) & res_adj$padj < 0.05 & abs(shr_adj$log2FoldChange) > 0.5
nas_rows[["F4 vs F0, ~ nas_missing + group (diagnostic)"]] <- c(
  case = sum(st == "4"), ref = sum(st == "0"), tested = nrow(res_adj), padj05 = sum(res_adj$padj < 0.05, na.rm = TRUE),
  up = sum(deg_adj & shr_adj$log2FoldChange > 0), down = sum(deg_adj & shr_adj$log2FoldChange < 0),
  median_disp = median(dispersions(dds_adj)[mcols(dds_adj)$baseMean > 10], na.rm = TRUE))
write.table(data.frame(ensembl_id = rownames(res_adj), log2FC_shrunk = shr_adj$log2FoldChange, pvalue = res_adj$pvalue,
                       padj = res_adj$padj, deg = deg_adj),
            file.path(out_dir, "DE_F4_vs_F0_adjusted_nas_missing_DIAGNOSTIC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

nas_summary <- data.frame(test = names(nas_rows), do.call(rbind, nas_rows), row.names = NULL, check.names = FALSE)
write.table(nas_summary, file.path(out_dir, "DE_nas_missing_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

m2 <- merge(main[, c("ensembl_id", "log2FC_shrunk", "deg")], f0_na$tab[, c("ensembl_id", "log2FC_shrunk", "deg")],
            by = "ensembl_id", suffixes = c("_main", "_nas"))
nas_overlap <- c(
  main_deg = sum(m2$deg_main), nasF0_deg = sum(m2$deg_nas), both = sum(m2$deg_main & m2$deg_nas),
  both_same_direction = sum(m2$deg_main & m2$deg_nas & sign(m2$log2FC_shrunk_main) == sign(m2$log2FC_shrunk_nas)),
  spearman_lfc_all = round(cor(m2$log2FC_shrunk_main, m2$log2FC_shrunk_nas, method = "spearman"), 3),
  spearman_lfc_main_deg = round(cor(m2$log2FC_shrunk_main[m2$deg_main], m2$log2FC_shrunk_nas[m2$deg_main], method = "spearman"), 3))
write.table(data.frame(metric = names(nas_overlap), value = nas_overlap), file.path(out_dir, "main_vs_nas_missing_effect_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p4 <- ggplot(m2, aes(log2FC_shrunk_nas, log2FC_shrunk_main, colour = deg_main)) +
  geom_point(size = 0.6, alpha = 0.4) + geom_abline(linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#e34948"), name = "main DEG") +
  labs(x = "log2FC  F0 NAS-missing vs F0 NAS-present  (technical axis, same stage)",
       y = "log2FC  main F4 vs F0",
       title = sprintf("NAS-missing effect vs main contrast (Spearman rho all genes = %.2f)", nas_overlap["spearman_lfc_all"]))
ggsave(file.path(out_dir, "nas_missing_effect_vs_main_lfc.png"), p4, width = 7, height = 6, dpi = 150)

p5 <- ggplot(samples, aes(PC1, PC2, colour = stage, shape = nas_NA)) + geom_point(size = 2.4, alpha = 0.85) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), name = "NAS missing") +
  labs(x = sprintf("PC1 (%.1f%%)", pve[1]), y = sprintf("PC2 (%.1f%%)", pve[2]),
       title = "GSE162694 PCA coloured by stage, shape = NAS score missing")
ggsave(file.path(out_dir, "pca_stage_nas_missing.png"), p5, width = 7, height = 5.5, dpi = 150)

# ---- console report ----
options(width = 200)
cat("\n=== sample blocks x stage ===\n"); print(table(samples$block, samples$stage))
cat("\n=== PC variance explained by covariates (R^2) ===\n"); print(r2, row.names = FALSE)
cat("\n=== library QC medians by block ===\n")
print(aggregate(cbind(lib_M = lib_size / 1e6, detected_ge10, top100_frac, size_factor) ~ block, samples, median))
cat("\n=== DE diagnostics (same DESeq2 pipeline, DEG = padj<0.05 & |shrunk LFC|>0.5) ===\n")
print(de_summary, row.names = FALSE)
cat(sprintf("\nmedian dispersion (baseMean>10): GSE135251 main F4 vs F0 = %.3f\n", disp135))
cat("\n=== main F4 vs F0 DEG vs block effect ===\n"); print(overlap)
cat("\n=== NAS-missing diagnostics ===\n"); print(nas_summary, row.names = FALSE)
cat("\n=== main F4 vs F0 DEG vs NAS-missing effect (within F0) ===\n"); print(nas_overlap)
cat("\noutputs:", out_dir, "\n")
