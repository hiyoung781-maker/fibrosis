#!/usr/bin/env python3
"""GATE 2 -- independent replication of GSE135251 DEG in GSE162694.

Run with: ~/miniconda3/envs/geneformer/bin/python scripts/04_gate2_replication.py
(the `geo` env lacks scipy; `geneformer` env has pandas/numpy/scipy/matplotlib)

Inputs (already computed independently per-cohort by 01_deseq2_GSE135251.R
and 03_deseq2_GSE162694.R, no batch merging / no ComBat):
  results/deseq2/GSE135251_main_F4_vs_F0.tsv   (discovery cohort)
  results/deseq2/GSE162694_main_F4_vs_F0.tsv   (replication cohort)

plan.md GATE 2 criteria (pre-registered, do not adjust):
  - metric scope: top 500 DEG by discovery (GSE135251) p-value
  - PASS if direction concordance >= 70% (binomial test vs 0.5, p < 1e-10)
    AND Spearman rho >= 0.3
  - RRHO reported as a supplementary/visual check, not a numeric gate criterion
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest, hypergeom
import matplotlib.pyplot as plt

BASE = "/home/seyoung/academic_symposium/fibrosis"
DESEQ_DIR = os.path.join(BASE, "results/deseq2")
OUT_DIR = os.path.join(BASE, "results/gate2")
os.makedirs(OUT_DIR, exist_ok=True)

TOP_N = 500
CONCORDANCE_THRESH = 0.70
BINOM_P_THRESH = 1e-10
SPEARMAN_THRESH = 0.30

COL_UP = "#e34948"
COL_DOWN = "#2a78d6"
COL_CONC = "#2a9d5c"
COL_DISC = "#c3c2b7"
COL_TEXT = "#0b0b0b"
COL_SEC_TEXT = "#52514e"
COL_GRID = "#e1e0d9"
COL_SURFACE = "#fcfcfb"

# ---------------------------------------------------------------------------
# 1. load & merge
# ---------------------------------------------------------------------------
a = pd.read_csv(os.path.join(DESEQ_DIR, "GSE135251_main_F4_vs_F0.tsv"), sep="\t").set_index("ensembl_id")
b = pd.read_csv(os.path.join(DESEQ_DIR, "GSE162694_main_F4_vs_F0.tsv"), sep="\t").set_index("ensembl_id")

common = a.index.intersection(b.index)
merged = pd.DataFrame({
    "lfc_discovery": a.loc[common, "log2FC_shrunk"],
    "p_discovery": a.loc[common, "pvalue"],
    "padj_discovery": a.loc[common, "padj"],
    "lfc_replication": b.loc[common, "log2FC_shrunk"],
    "p_replication": b.loc[common, "pvalue"],
    "padj_replication": b.loc[common, "padj"],
})
merged.index.name = "ensembl_id"

n_a, n_b, n_common = len(a), len(b), len(common)

# ---------------------------------------------------------------------------
# 2. top-500 direction concordance + binomial test
# ---------------------------------------------------------------------------
top500 = merged.sort_values("p_discovery").head(TOP_N).copy()
top500["concordant"] = np.sign(top500["lfc_discovery"]) == np.sign(top500["lfc_replication"])
n_conc = int(top500["concordant"].sum())
n_top = len(top500)
concordance = n_conc / n_top
binom = binomtest(n_conc, n_top, 0.5, alternative="greater")

rho_top500, rho_top500_p = spearmanr(top500["lfc_discovery"], top500["lfc_replication"])
rho_all, rho_all_p = spearmanr(merged["lfc_discovery"], merged["lfc_replication"])

top500.sort_values("p_discovery").to_csv(os.path.join(OUT_DIR, "GATE2_top500_direction_concordance.csv"))
merged.to_csv(os.path.join(OUT_DIR, "GATE2_merged_common_genes.csv"))

# ---------------------------------------------------------------------------
# 3. RRHO (rank-rank hypergeometric overlap), signed score = sign(lfc)*-log10(p)
#    full genome-wide two-sided map: rank 1 = most up-regulated ... rank N = most down-regulated
# ---------------------------------------------------------------------------
def signed_score(lfc, p):
    return -np.log10(p.clip(lower=1e-300)) * np.sign(lfc)

merged["score_discovery"] = signed_score(merged["lfc_discovery"], merged["p_discovery"])
merged["score_replication"] = signed_score(merged["lfc_replication"], merged["p_replication"])

# rank 1 = highest (most up-regulated) score
rank_a = merged["score_discovery"].rank(ascending=False, method="first").astype(int)
rank_b = merged["score_replication"].rank(ascending=False, method="first").astype(int)

N = n_common
STEP = max(1, N // 150)  # ~150x150 grid
steps = np.arange(STEP, N + 1, STEP)
if steps[-1] != N:
    steps = np.append(steps, N)

# arr[k] = rank_b of the gene whose rank_a == k+1  (index by discovery rank)
order_by_a = rank_a.sort_values().index
arr = rank_b.loc[order_by_a].to_numpy()

rrho_logp = np.zeros((len(steps), len(steps)))
for ii, i in enumerate(steps):
    sub = np.sort(arr[:i])
    counts_le_j = np.searchsorted(sub, steps, side="right")
    for jj, j in enumerate(steps):
        overlap = int(counts_le_j[jj])
        pval = hypergeom.sf(overlap - 1, N, i, j)
        rrho_logp[ii, jj] = -np.log10(max(pval, 1e-300))

rrho_df = pd.DataFrame(rrho_logp, index=steps, columns=steps)
rrho_df.to_csv(os.path.join(OUT_DIR, "GATE2_RRHO_neglog10p_grid.csv"))

# corner-summary table (classic reporting quadrants) for quick numeric read
def quadrant_overlap(rank_x, rank_y, k, ascending_x=True, ascending_y=True):
    rx = rank_x if ascending_x else (N + 1 - rank_x)
    ry = rank_y if ascending_y else (N + 1 - rank_y)
    genes = (rx <= k) & (ry <= k)
    overlap = int(genes.sum())
    pval = hypergeom.sf(overlap - 1, N, k, k)
    return overlap, pval

quad_rows = []
for k in (100, 300, 500, 1000):
    ov, pv = quadrant_overlap(rank_a, rank_b, k, True, True)
    quad_rows.append({"quadrant": "up_vs_up", "k": k, "overlap": ov, "expected": k * k / N, "hypergeom_p": pv})
    ov, pv = quadrant_overlap(rank_a, rank_b, k, False, False)
    quad_rows.append({"quadrant": "down_vs_down", "k": k, "overlap": ov, "expected": k * k / N, "hypergeom_p": pv})
    ov, pv = quadrant_overlap(rank_a, rank_b, k, True, False)
    quad_rows.append({"quadrant": "up_vs_down (discordant)", "k": k, "overlap": ov, "expected": k * k / N, "hypergeom_p": pv})
    ov, pv = quadrant_overlap(rank_a, rank_b, k, False, True)
    quad_rows.append({"quadrant": "down_vs_up (discordant)", "k": k, "overlap": ov, "expected": k * k / N, "hypergeom_p": pv})
quad_df = pd.DataFrame(quad_rows)
quad_df.to_csv(os.path.join(OUT_DIR, "GATE2_RRHO_quadrant_summary.csv"), index=False)

# ---------------------------------------------------------------------------
# 4. GATE 2 verdict
# ---------------------------------------------------------------------------
pass_concordance = (concordance >= CONCORDANCE_THRESH) and (binom.pvalue < BINOM_P_THRESH)
pass_spearman = rho_top500 >= SPEARMAN_THRESH
gate2_pass = pass_concordance and pass_spearman

summary = {
    "n_tested_GSE135251": n_a,
    "n_tested_GSE162694": n_b,
    "n_common_genes": n_common,
    "top_n": n_top,
    "top500_n_concordant": n_conc,
    "top500_direction_concordance": concordance,
    "top500_binomial_p_vs_0.5": binom.pvalue,
    "criterion_concordance_ge_0.70_and_p_lt_1e-10": pass_concordance,
    "top500_spearman_rho": rho_top500,
    "top500_spearman_p": rho_top500_p,
    "criterion_spearman_rho_ge_0.30": pass_spearman,
    "allgenes_spearman_rho": rho_all,
    "allgenes_spearman_p": rho_all_p,
    "GATE2_VERDICT": "PASS" if gate2_pass else "FAIL",
}
pd.Series(summary).to_csv(os.path.join(OUT_DIR, "GATE2_summary_metrics.csv"), header=False)

report_lines = [
    "GATE 2 -- Independent Replication (GSE135251 discovery vs GSE162694 replication)",
    "=" * 78,
    f"Genes tested in GSE135251 (discovery): {n_a:,}",
    f"Genes tested in GSE162694 (replication): {n_b:,}",
    f"Common genes (both tested): {n_common:,}",
    "",
    f"[Top {n_top} DEG by GSE135251 discovery p-value]",
    f"  Direction concordance: {n_conc}/{n_top} = {concordance:.4f}",
    f"  Binomial test vs 0.5 (one-sided, greater): p = {binom.pvalue:.3e}",
    f"  -> concordance criterion (>=0.70, p<1e-10): {'PASS' if pass_concordance else 'FAIL'}",
    "",
    f"  Spearman rho (log2FC_shrunk, top {n_top}): {rho_top500:.4f} (p={rho_top500_p:.3e})",
    f"  -> spearman criterion (rho>=0.30): {'PASS' if pass_spearman else 'FAIL'}",
    "",
    f"  [supplementary] Spearman rho over all {n_common:,} common genes: {rho_all:.4f} (p={rho_all_p:.3e})",
    "",
    "[RRHO quadrant summary] (see GATE2_RRHO_quadrant_summary.csv, GATE2_RRHO_heatmap.png for full map)",
]
for _, row in quad_df.iterrows():
    report_lines.append(
        f"  {row['quadrant']:<26} top{int(row['k']):<5} overlap={int(row['overlap']):4d} "
        f"(expected {row['expected']:.1f})  p={row['hypergeom_p']:.3e}"
    )
report_lines += [
    "",
    "=" * 78,
    f"GATE 2 VERDICT: {'PASS' if gate2_pass else 'FAIL'}",
    "=" * 78,
]
report_txt = "\n".join(report_lines)
with open(os.path.join(OUT_DIR, "GATE2_report.txt"), "w") as f:
    f.write(report_txt + "\n")
print(report_txt)

# ---------------------------------------------------------------------------
# 5. plots
# ---------------------------------------------------------------------------
# 5a. log2FC scatter (all common genes, top500 highlighted, concordant/discordant colored)
fig, ax = plt.subplots(figsize=(7, 7), dpi=150)
fig.patch.set_facecolor(COL_SURFACE)
ax.set_facecolor(COL_SURFACE)

rest = merged.drop(top500.index)
ax.scatter(rest["lfc_discovery"], rest["lfc_replication"], s=4, c=COL_DISC, alpha=0.25,
           linewidths=0, label=f"other tested genes (n={len(rest):,})")
conc_pts = top500[top500["concordant"]]
disc_pts = top500[~top500["concordant"]]
ax.scatter(disc_pts["lfc_discovery"], disc_pts["lfc_replication"], s=14, c=COL_DOWN, alpha=0.8,
           linewidths=0, label=f"top{TOP_N} discordant (n={len(disc_pts)})")
ax.scatter(conc_pts["lfc_discovery"], conc_pts["lfc_replication"], s=14, c=COL_CONC, alpha=0.8,
           linewidths=0, label=f"top{TOP_N} concordant (n={len(conc_pts)})")

lim = max(merged["lfc_discovery"].abs().max(), merged["lfc_replication"].abs().max()) * 1.05
ax.axhline(0, color=COL_GRID, linewidth=0.8)
ax.axvline(0, color=COL_GRID, linewidth=0.8)
ax.plot([-lim, lim], [-lim, lim], color=COL_SEC_TEXT, linestyle="--", linewidth=0.8, alpha=0.6)
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_xlabel("log2FC (shrunk) -- GSE135251 (discovery)", color=COL_TEXT)
ax.set_ylabel("log2FC (shrunk) -- GSE162694 (replication)", color=COL_TEXT)
ax.set_title(
    f"GATE 2 replication: top{TOP_N} direction concordance = {concordance:.1%}, "
    f"Spearman rho = {rho_top500:.3f}\nVerdict: {'PASS' if gate2_pass else 'FAIL'}",
    color=COL_TEXT, fontsize=11,
)
ax.grid(True, color=COL_GRID, linewidth=0.5)
ax.set_axisbelow(True)
for spine in ax.spines.values():
    spine.set_color(COL_GRID)
ax.tick_params(colors=COL_SEC_TEXT)
legend = ax.legend(loc="upper left", frameon=False, fontsize=8)
for text in legend.get_texts():
    text.set_color(COL_SEC_TEXT)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "GATE2_logFC_scatter.png"), facecolor=fig.get_facecolor())
plt.close(fig)

# 5b. RRHO heatmap
fig, ax = plt.subplots(figsize=(7.5, 6.5), dpi=150)
fig.patch.set_facecolor(COL_SURFACE)
im = ax.pcolormesh(steps, steps, rrho_logp.T, shading="auto", cmap="inferno")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("-log10(hypergeometric p)", color=COL_TEXT)
cbar.ax.yaxis.set_tick_params(color=COL_SEC_TEXT)
plt.setp(cbar.ax.get_yticklabels(), color=COL_SEC_TEXT)
ax.set_xlabel("rank in GSE135251 (discovery)\n1=most up-regulated -> N=most down-regulated", color=COL_TEXT)
ax.set_ylabel("rank in GSE162694 (replication)\n1=most up-regulated -> N=most down-regulated", color=COL_TEXT)
ax.set_title(f"GATE 2 RRHO map (N={N:,} common genes, step={STEP})", color=COL_TEXT, fontsize=11)
ax.tick_params(colors=COL_SEC_TEXT)
for spine in ax.spines.values():
    spine.set_color(COL_GRID)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "GATE2_RRHO_heatmap.png"), facecolor=fig.get_facecolor())
plt.close(fig)

print(f"\nSaved to {OUT_DIR}/:")
for fn in [
    "GATE2_top500_direction_concordance.csv",
    "GATE2_merged_common_genes.csv",
    "GATE2_RRHO_neglog10p_grid.csv",
    "GATE2_RRHO_quadrant_summary.csv",
    "GATE2_summary_metrics.csv",
    "GATE2_report.txt",
    "GATE2_logFC_scatter.png",
    "GATE2_RRHO_heatmap.png",
]:
    print(f"  {fn}")
