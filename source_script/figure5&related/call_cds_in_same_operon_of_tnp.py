import pandas as pd
import numpy as np
import os

BASE_DIR = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats"
ANNOT_FILE = os.path.join(BASE_DIR, "allpp_annot/allpharokka_annot.tsv")
META_FILE = os.path.join(BASE_DIR, "allpp_checkv_quality/hqphage_metadata.tsv")
OUTPUT_FILE = "Transposase_Adjacent_Operon_Genes.tsv"
df_meta = pd.read_csv(META_FILE, sep='\t', low_memory=False)

df_meta['ictv_Class'] = df_meta['ictv_Class'].fillna('').astype(str)
caudo_phages = df_meta[df_meta['ictv_Class'].str.contains('Caudo', case=False)]['Phage_Name'].unique()
df_cds = pd.read_csv(ANNOT_FILE, sep='\t', low_memory=False)

df_cds['contig_clean'] = df_cds['contig'].str.replace(r'_(\d+)-(\d+)$', r'_\1_\2', regex=True)
caudo_phages_clean = [str(x).replace('-', '_') for x in caudo_phages] 
df_cds = df_cds[df_cds['contig'].isin(caudo_phages) | df_cds['contig_clean'].isin(caudo_phages_clean)].copy()

df_cds['phys_start'] = df_cds[['start', 'stop']].min(axis=1)
df_cds = df_cds.sort_values(by=['contig', 'phys_start']).reset_index(drop=True)

df_cds['annot_lower'] = df_cds['annot'].fillna('').astype(str).str.lower()
df_cds['is_transposase'] = df_cds['annot_lower'].str.contains('transposase|transposon')
adjacent_genes = []

for contig, group in df_cds.groupby('contig'):
    if not group['is_transposase'].any():
        continue
        
    group = group.reset_index(drop=True)
    n_genes = len(group)
    tnp_indices = group.index[group['is_transposase']].tolist()
    
    for idx in tnp_indices:
        strand = group.loc[idx, 'frame']
        
        if idx > 0:
            prev_idx = idx - 1
            if group.loc[prev_idx, 'frame'] == strand and not group.loc[prev_idx, 'is_transposase']:
                adjacent_genes.append(group.loc[prev_idx])
                
        if idx < n_genes - 1:
            next_idx = idx + 1
            if group.loc[next_idx, 'frame'] == strand and not group.loc[next_idx, 'is_transposase']:
                adjacent_genes.append(group.loc[next_idx])
df_adj = pd.DataFrame(adjacent_genes).drop_duplicates(subset=['gene'])
df_unique_per_contig = df_adj.drop_duplicates(subset=['contig', 'phrog', 'category', 'annot'])
df_stats = df_unique_per_contig.groupby(['phrog', 'category', 'annot']).size().reset_index(name='Contig_Count')
total_valid_contigs = df_adj['contig'].nunique() 
df_stats['Genome_Frequency(%)'] = (df_stats['Contig_Count'] / total_valid_contigs) * 100
df_stats = df_stats.sort_values(by='Contig_Count', ascending=False).reset_index(drop=True)
df_stats.to_csv(OUTPUT_FILE, sep='\t', index=False)

top10_to_print = df_stats[['phrog', 'category', 'annot', 'Count', 'Frequency(%)']].head(10)
top10_to_print['Frequency(%)'] = top10_to_print['Frequency(%)'].round(2).astype(str) + '%'

print(top10_to_print.to_string(index=False))

