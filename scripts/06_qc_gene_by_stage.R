suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# Per-sample expression-by-stage QC for genes whose DESeq2 direction needs checking against the
# plan.md §2 assumption ("F4에서 강하게 증가"). Checks raw/CPM/normalized expression by stage,
# zero inflation, co-expression with COL1A1, and replication in GSE162694.
#
# Usage: Rscript 06_qc_gene_by_stage.R <panel> SYMBOL=ENSG [SYMBOL=ENSG ...]
#   Rscript 06_qc_gene_by_stage.R POSTN POSTN=ENSG00000133110
#   Rscript 06_qc_gene_by_stage.R LOX LOX=ENSG00000113083 LOXL1=ENSG00000129038 LOXL2=ENSG00000134013 \
#     LOXL3=ENSG00000115318 LOXL4=ENSG00000138131

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2, all(grepl("^[A-Za-z0-9-]+=ENSG[0-9]{11}$", args[-1])))
panel <- args[1]
kv <- strsplit(args[-1], "=")
targets <- setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))

base <- "/home/seyoung/academic_symposium/fibrosis"
out_dir <- file.path(base, "results/qc", paste0(panel, "_check"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref_markers <- c(COL1A1 = "ENSG00000108821", COL3A1 = "ENSG00000168542",
                 LUM = "ENSG00000139329", THBS2 = "ENSG00000186340", ALB = "ENSG00000163631")
genes <- c(targets, ref_markers[!names(ref_markers) %in% names(targets)])
stage_levels <- c("control", "0", "1", "2", "3", "4")

# ---- GSE135251 ----
s1 <- read.delim(file.path(base, "geo_data/GSE135251/GSE135251_samples.tsv"), stringsAsFactors = FALSE)
s1 <- s1[s1$exclude != "True", ]
read_counts <- function(f) {
  read.delim(gzfile(file.path(base, "geo_data/GSE135251/RAW", f)), header = FALSE,
             col.names = c("gene", "count"), colClasses = c("character", "integer"))
}
tabs <- lapply(s1$file, read_counts)
c1 <- do.call(cbind, lapply(tabs, `[[`, "count"))
dimnames(c1) <- list(tabs[[1]]$gene, s1$gsm)
s1$stage <- factor(ifelse(s1$group_in_paper == "control", "control", s1$fibrosis_stage), stage_levels)
stopifnot(ncol(c1) == 214, all(genes %in% rownames(c1)))

# ---- GSE162694 ----
s2 <- read.delim(file.path(base, "geo_data/GSE162694/GSE162694_samples.tsv"), stringsAsFactors = FALSE)
c2 <- as.matrix(read.csv(file.path(base, "geo_data/GSE162694/GSE162694_raw_counts.csv.gz"),
                         row.names = 1, check.names = FALSE))[, s2$col_id]
mode(c2) <- "integer"
colnames(c2) <- s2$gsm
s2$stage <- factor(ifelse(s2$fibrosis_stage == "normal", "control", s2$fibrosis_stage), stage_levels)
stopifnot(ncol(c2) == 143, all(genes %in% rownames(c2)))

long_table <- function(counts, samples, cohort) {
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(counts, samples, ~ 1))
  norm <- counts(dds, normalized = TRUE)
  lib <- colSums(counts)
  do.call(rbind, lapply(names(genes), function(g) {
    id <- genes[[g]]
    data.frame(cohort = cohort, gsm = samples$gsm, stage = samples$stage, gene = g, ensembl_id = id,
               raw = counts[id, ], cpm = counts[id, ] / lib * 1e6, norm = norm[id, ],
               lib_size = lib, size_factor = sizeFactors(dds))
  }))
}
d <- rbind(long_table(c1, s1, "GSE135251"), long_table(c2, s2, "GSE162694"))
d$gene <- factor(d$gene, names(genes))
write.table(d, file.path(out_dir, paste0(panel, "_per_sample.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- summaries ----
summ <- aggregate(cbind(raw, norm) ~ cohort + gene + stage, d, function(x) round(median(x), 1))
names(summ)[4:5] <- c("median_raw", "median_norm")
summ$n <- aggregate(raw ~ cohort + gene + stage, d, length)$raw
summ$frac_zero <- round(aggregate(raw ~ cohort + gene + stage, d, function(x) mean(x == 0))$raw, 3)
write.table(summ, file.path(out_dir, paste0(panel, "_median_by_stage.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

col1a1 <- d[d$gene == "COL1A1", c("cohort", "gsm", "norm")]
names(col1a1)[3] <- "norm_COL1A1"
stats <- do.call(rbind, lapply(split(d, list(d$cohort, d$gene), drop = TRUE), function(x) {
  f4 <- x[x$stage == "4", ]; f0 <- x[x$stage == "0", ]
  nafld <- x[x$stage != "control", ]
  xc <- merge(x, col1a1)
  data.frame(cohort = x$cohort[1], gene = x$gene[1], ensembl_id = x$ensembl_id[1],
             n_F4 = nrow(f4), n_F0 = nrow(f0),
             median_raw_F4 = median(f4$raw), median_raw_F0 = median(f0$raw),
             median_norm_F4 = median(f4$norm), median_norm_F0 = median(f0$norm),
             log2_median_ratio = log2(median(f4$norm) / median(f0$norm)),
             wilcox_p = wilcox.test(f4$norm, f0$norm)$p.value,
             spearman_rho_stage = cor(nafld$norm, as.integer(as.character(nafld$stage)), method = "spearman"),
             spearman_rho_COL1A1 = cor(xc$norm, xc$norm_COL1A1, method = "spearman"),
             frac_zero_F0 = mean(f0$raw == 0), frac_zero_F4 = mean(f4$raw == 0))
}))
write.table(stats, file.path(out_dir, paste0(panel, "_F4_vs_F0_stats.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
print(stats[, -3], digits = 3, row.names = FALSE)

# ---- plots ----
stage_cols <- c(control = "#9aa0a6", "0" = "#4e79a7", "1" = "#76b7b2", "2" = "#edc948", "3" = "#f28e2b", "4" = "#e15759")
theme_set(theme_bw(base_size = 11) + theme(legend.position = "none", strip.background = element_rect(fill = "grey95")))
box_layers <- list(geom_boxplot(outlier.shape = NA, alpha = 0.6),
                   geom_jitter(width = 0.2, size = 0.7, alpha = 0.6),
                   scale_fill_manual(values = stage_cols))

td <- d[d$gene %in% names(targets), ]
td_long <- rbind(transform(td, scale = "raw count", value = log2(raw + 1)),
                 transform(td, scale = "CPM", value = log2(cpm + 1)),
                 transform(td, scale = "DESeq2 normalized", value = log2(norm + 1)))
td_long$scale <- factor(td_long$scale, c("raw count", "CPM", "DESeq2 normalized"))
p1 <- ggplot(td_long, aes(stage, value, fill = stage)) + box_layers +
  facet_grid(gene ~ cohort + scale, scales = "free_y") +
  labs(title = sprintf("%s per-sample expression by fibrosis stage (all log2 +1)", panel),
       x = "Fibrosis stage (control = histologically normal liver)", y = NULL)
ggsave(file.path(out_dir, paste0(panel, "_boxplot_by_stage.png")), p1,
       width = 16, height = 1.2 + 2.6 * length(targets), dpi = 150, limitsize = FALSE)

md <- d[d$gene != "ALB", ]
p2 <- ggplot(md, aes(stage, log2(norm + 1), fill = stage)) + box_layers +
  facet_grid(cohort ~ gene, scales = "free_y") +
  labs(title = sprintf("%s vs canonical fibrosis markers (DESeq2 normalized, log2 +1)", panel),
       x = "Fibrosis stage", y = NULL)
ggsave(file.path(out_dir, paste0(panel, "_vs_fibrosis_markers.png")), p2,
       width = 2.6 * nlevels(droplevels(md$gene)), height = 6, dpi = 150, limitsize = FALSE)

sd <- merge(td, col1a1)
p3 <- ggplot(sd, aes(log2(norm_COL1A1 + 1), log2(norm + 1), colour = stage)) +
  geom_point(size = 1, alpha = 0.8) +
  facet_grid(cohort ~ gene) +
  scale_colour_manual(values = stage_cols) +
  theme(legend.position = "right") +
  labs(title = sprintf("%s vs COL1A1 co-expression", panel), x = "COL1A1 log2(norm +1)", y = "log2(norm +1)")
ggsave(file.path(out_dir, paste0(panel, "_vs_COL1A1_scatter.png")), p3,
       width = 2 + 3 * length(targets), height = 6, dpi = 150)
