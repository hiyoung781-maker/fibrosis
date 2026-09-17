import glob
import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

BASE = '/home/seyoung/academic_symposium/fibrosis'
DEG_FILE = f'{BASE}/results/deseq2/GSE135251_main_F4_vs_F0.tsv'
OUT_DIR = f'{BASE}/results/gate1'
# Highest-evidence tier wins, so tier B means ligand/pocket evidence without any clinical precedent.
TIERS = {
    'A': ['Approved Drug', 'Advanced Clinical', 'Phase 1 Clinical'],
    'B': ['Structure with Ligand', 'High-Quality Ligand', 'High-Quality Pocket'],
    'C': ['Med-Quality Pocket', 'Druggable Family'],
}
PASS_B, PASS_ABC = 10, 20
CONTROLS = {
    'positive': ['THRB', 'PNPLA3', 'HSD17B13', 'PPARA', 'PPARD', 'PPARG'],
    'negative': ['LOXL2', 'MAP3K5', 'CCR2', 'CCR5', 'CASP1', 'CASP3', 'CASP8', 'SERPINH1'],
}


def sm_tier(tract):
    sm = [x['id'] for x in tract if x['modality'] == 'SM' and x['value']] if tract is not None else []
    tier = next((t for t, cats in TIERS.items() if set(cats) & set(sm)), 'none')
    return tier, ';'.join(sm)


files = sorted(glob.glob(f'{BASE}/opentarget/target/*.parquet'))
ot = pa.concat_tables(
    [pq.read_table(f, columns=['id', 'approvedSymbol', 'biotype', 'tractability']) for f in files]
).to_pandas()
ot['sm_tier'], ot['sm_categories'] = zip(*(sm_tier(t) for t in ot.tractability))
# Annotation only, never used for ranking: strong = HQ Ligand or HQ Pocket, weak = any other SM evidence.
ot['sm_evidence'] = [
    '' if not c else 'strong' if ('High-Quality Ligand' in c or 'High-Quality Pocket' in c) else 'weak'
    for c in ot.sm_categories
]

res = pd.read_csv(DEG_FILE, sep='\t')
tab = res.merge(ot.drop(columns='tractability'), left_on='ensembl_id', right_on='id', how='left').drop(columns='id')
tab['sm_tier'] = tab.sm_tier.fillna('none')
assert len(tab) == len(res)

os.makedirs(OUT_DIR, exist_ok=True)
tab.to_csv(f'{OUT_DIR}/GSE135251_main_tractability.tsv', sep='\t', index=False)

deg = tab[tab.deg]
print(f'tested genes: {len(tab)} | not in Open Targets: {tab.approvedSymbol.isna().sum()}')
print(f'DEG: {len(deg)} (protein_coding {(deg.biotype == "protein_coding").sum()})')
print('\nDEG by SM tier:')
print(pd.crosstab(deg.sm_tier, deg.log2FC_shrunk.gt(0).map({True: 'up', False: 'down'}), margins=True).to_string())

n_b = (deg.sm_tier == 'B').sum()
n_abc = deg.sm_tier.isin(['A', 'B', 'C']).sum()
verdict = 'PASS' if n_b >= PASS_B and n_abc >= PASS_ABC else 'FAIL'
print(f'\nGATE 1: tier B = {n_b} (need >= {PASS_B}), tier A-C = {n_abc} (need >= {PASS_ABC}) -> {verdict}')

print('\ntier B DEG by sm_evidence:', deg[deg.sm_tier == 'B'].sm_evidence.value_counts().to_dict())

cols = ['approvedSymbol', 'log2FC_shrunk', 'padj', 'baseMean', 'sm_evidence', 'sm_categories']
for tier in ['B', 'A']:
    sub = deg[deg.sm_tier == tier].sort_values('padj')
    sub[cols].to_csv(f'{OUT_DIR}/GSE135251_main_DEG_tier{tier}.tsv', sep='\t', index=False)
    print(f'\n-- tier {tier} DEG (top 25 of {len(sub)} by padj) --')
    print(sub[cols].head(25).to_string(index=False, float_format=lambda x: f'{x:.3g}'))

print('\n-- control genes --')
by_sym = tab.set_index('approvedSymbol')
for kind, genes in CONTROLS.items():
    for g in genes:
        if g in by_sym.index:
            r = by_sym.loc[g]
            print(f'{kind:<8} {g:<9} tier={r.sm_tier:<4} DEG={r.deg!s:<5} log2FC_shrunk={r.log2FC_shrunk:+.2f}')
        else:
            print(f'{kind:<8} {g:<9} not tested (low expression)')
