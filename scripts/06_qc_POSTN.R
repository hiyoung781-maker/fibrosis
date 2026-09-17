suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# POSTN (ENSG00000133110) comes out DOWN in GSE135251 F4 vs F0 (shrunk log2FC -1.29), opposite to
# plan.md §2 negative-control family 2 ("F4에서 강하게 증가"). Check raw/normalized per-sample
# expression by stage, count-level artefacts, and whether the direction replicates in GSE162694.

base <- "/home/seyoung/academic_symposium/fibrosis"
out_dir <- file.path(base, "results/qc/POSTN_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

POSTN <- "ENSG00000133110"
markers <- c(POSTN = POSTN, COL1A1 = "ENSG00000108821", COL3A1 = "ENSG00000168542",
             LUM = "ENSG00000139329", THBS2 = "ENSG00000186340", ALB = "ENSG00000163631")
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
stopifnot(ncol(c1) == 214, all(markers %in% rownames(c1)))

# ---- GSE162694 ----
s2 <- read.delim(file.path(base, "geo_data/GSE162694/GSE162694_samples.tsv"), stringsAsFactors = FALSE)
c2 <- as.matrix(read.csv(file.path(base, "geo_data/GSE162694/GSE162694_raw_counts.csv.gz"),
                         row.names = 1, check.names = FALSE))[, s2$col_id]
colnames(c2) <- s2$gsm
s2$stage <- factor(ifelse(s2$fibrosis_stage == "normal", "control", s2$fibrosis_stage), stage_levels)
stopifnot(ncol(c2) == 143, all(markers %in% rownames(c2)))

long_table <- function(counts, samples, cohort) {
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(round(counts), samples, ~ 1))
  norm <- counts(dds, normalized = TRUE)
  lib <- colSums(counts)
  do.call(rbind, lapply(names(markers), function(g) {
    id <- markers[[g]]
    data.frame(cohort = cohort, gsm = samples$gsm, stage = samples$stage, gene = g,
               raw = counts[id, ], cpm = counts[id, ] / lib * 1e6, norm = norm[id, ],
               lib_size = lib, size_factor = sizeFactors(dds))
  }))
}
d <- rbind(long_table(c1, s1, "GSE135251"), long_table(c2, s2, "GSE162694"))
write.table(d, file.path(out_dir, "POSTN_markers_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- summaries ----
summ <- aggregate(cbind(raw, norm) ~ cohort + gene + stage, d, function(x) round(median(x), 1))
names(summ)[4:5] <- c("median_raw", "median_norm")
summ$n <- aggregate(raw ~ cohort + gene + stage, d, length)$raw
write.table(summ, file.path(out_dir, "POSTN_markers_median_by_stage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

stats <- do.call(rbind, lapply(split(d, list(d$cohort, d$gene), drop = TRUE), function(x) {
  f4 <- x$norm[x$stage == "4"]; f0 <- x$norm[x$stage == "0"]
  nafld <- x[x$stage != "control", ]
  data.frame(cohort = x$cohort[1], gene = x$gene[1], n_F4 = length(f4), n_F0 = length(f0),
             median_F4 = median(f4), median_F0 = median(f0),
             log2_median_ratio = log2(median(f4) / median(f0)),
             wilcox_p = wilcox.test(f4, f0)$p.value,
             spearman_rho_stage = cor(nafld$norm, as.integer(as.character(nafld$stage)), method = "spearman"),
             frac_zero = mean(x$raw == 0))
}))
write.table(stats, file.path(out_dir, "POSTN_markers_F4_vs_F0_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(stats, digits = 3, row.names = FALSE)

# POSTN vs COL1A1 co-expression within each cohort (does POSTN track the fibrogenic program at all?)
wide <- reshape(d[d$gene %in% c("POSTN", "COL1A1"), c("cohort", "gsm", "stage", "gene", "norm")],
                idvar = c("cohort", "gsm", "stage"), timevar = "gene", direction = "wide")
for (co in unique(wide$cohort)) {
  w <- wide[wide$cohort == co, ]
  cat(sprintf("%s  Spearman POSTN~COL1A1 = %.3f\n", co, cor(w$norm.POSTN, w$norm.COL1A1, method = "spearman")))
}

# ---- plots ----
stage_cols <- c(control = "#9aa0a6", "0" = "#4e79a7", "1" = "#76b7b2", "2" = "#edc948", "3" = "#f28e2b", "4" = "#e15759")
theme_set(theme_bw(base_size = 11) + theme(legend.position = "none", strip.background = element_rect(fill = "grey95")))

pd <- d[d$gene == "POSTN", ]
pd_long <- rbind(
  transform(pd, scale = "raw count (log10 +1)", value = log10(raw + 1)),
  transform(pd, scale = "CPM (log2 +1)", value = log2(cpm + 1)),
  transform(pd, scale = "DESeq2 normalized (log2 +1)", value = log2(norm + 1)))
pd_long$scale <- factor(pd_long$scale, unique(pd_long$scale))
p1 <- ggplot(pd_long, aes(stage, value, fill = stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  facet_grid(scale ~ cohort, scales = "free_y") +
  scale_fill_manual(values = stage_cols) +
  labs(title = "POSTN (ENSG00000133110) per-sample expression by fibrosis stage",
       x = "Fibrosis stage (control = histologically normal liver)", y = NULL)
ggsave(file.path(out_dir, "POSTN_boxplot_by_stage.png"), p1, width = 10, height = 9, dpi = 150)

md <- d[d$gene != "ALB", ]
md$gene <- factor(md$gene, names(markers)[names(markers) != "ALB"])
p2 <- ggplot(md, aes(stage, log2(norm + 1), fill = stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_grid(cohort ~ gene, scales = "free_y") +
  scale_fill_manual(values = stage_cols) +
  labs(title = "POSTN vs canonical fibrosis markers (DESeq2 normalized, log2 +1)", x = "Fibrosis stage", y = NULL)
ggsave(file.path(out_dir, "POSTN_vs_fibrosis_markers.png"), p2, width = 13, height = 6, dpi = 150)

p3 <- ggplot(wide, aes(log2(norm.COL1A1 + 1), log2(norm.POSTN + 1), colour = stage)) +
  geom_point(size = 1.2, alpha = 0.8) +
  facet_wrap(~ cohort) +
  scale_colour_manual(values = stage_cols) +
  theme(legend.position = "right") +
  labs(title = "POSTN vs COL1A1 co-expression", x = "COL1A1 log2(norm +1)", y = "POSTN log2(norm +1)")
ggsave(file.path(out_dir, "POSTN_vs_COL1A1_scatter.png"), p3, width = 10, height = 4.5, dpi = 150)
