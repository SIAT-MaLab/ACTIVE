import pandas as pd
import numpy as np
import os

BASE_DIR = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats"
ANNOT_DIR = os.path.join(BASE_DIR, "allpp_annot")
PHAROKKA_FILE = os.path.join(ANNOT_DIR, "allpharokka_annot.tsv")
PFAM_FILE = os.path.join(ANNOT_DIR, "pfamscan/all_merged_best_hits.tsv")
ISESCAN_RAW_FILE = os.path.join(ANNOT_DIR, "isescan_results/result_best_hits_cleaned.tsv")
COG_ANNOT_FILE = os.path.join(ANNOT_DIR, "COG2842_result.final.tsv")
PC_MAP_FILE = os.path.join(ANNOT_DIR, "allpp_PC_map.tsv")
METADATA_FILE = os.path.join(BASE_DIR, "merged_phage_stats_taxonomy.tsv")


DISCREPANCY_FILE = "annotation_method_discrepancies.tsv"
OUTPUT_FILE = "Phage_Lifestyle.tsv"

df_cds = pd.read_csv(PHAROKKA_FILE, sep='\t', low_memory=False)

df_cds['phys_start'] = df_cds[['start', 'stop']].min(axis=1)
df_cds = df_cds.sort_values(by=['contig', 'phys_start']).reset_index(drop=True)
df_cds['annot_lower'] = df_cds['annot'].fillna('').astype(str).str.lower()

integrase_pfams = {'PF07508', 'PF00239', 'PF00589'}
transposase_pfams = {'PF00665', 'PF13333', 'PF13683'}

pfam_gene_annot = {}
try:
    with open(PFAM_FILE, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split()
            if len(parts) > 5:
                gene_id = parts[0]
                pfam_id = parts[4].split('.')[0]
                if pfam_id in integrase_pfams:
                    pfam_gene_annot[gene_id] = 'integrase'
                elif pfam_id in transposase_pfams:
                    pfam_gene_annot[gene_id] = 'transposase'
except Exception as e:
    print(f"Warning: Error reading Pfam file: {e}")

pharokka_gene_annot = {}
for _, row in df_cds.iterrows():
    gene_id = row['gene']
    annot = row['annot_lower']
    is_int = 'integrase' in annot
    is_trans = 'transposase' in annot or 'transposon' in annot

    if is_int and not is_trans:
        pharokka_gene_annot[gene_id] = 'integrase'
    elif is_trans and not is_int:
        pharokka_gene_annot[gene_id] = 'transposase'
    elif is_int and is_trans:
        pharokka_gene_annot[gene_id] = 'both'

all_annotated_genes = set(pfam_gene_annot.keys()).union(set(pharokka_gene_annot.keys()))
comparison_records = []
reconciled_gene_annot = {}
gene_source_map = {}

for gene in all_annotated_genes:
    p_annot = pfam_gene_annot.get(gene, 'None')
    ph_annot = pharokka_gene_annot.get(gene, 'None')
    contig = gene.rsplit('_CDS_', 1)[0]

    status = 'Consistent'
    final_annot = 'None'
    source = 'None'

    if p_annot != 'None' and ph_annot == 'None':
        status = 'Pfam_Only'
        final_annot = p_annot
        source = 'Pfam'
    elif p_annot == 'None' and ph_annot != 'None':
        status = 'Pharokka_Only'
        final_annot = ph_annot
        source = 'Pharokka'
    elif p_annot != ph_annot:
        status = 'Conflict (Same CDS, Diff Annot)'
        final_annot = ph_annot
        source = 'Pharokka'
    else: 
        final_annot = p_annot
        source = 'Both'

    reconciled_gene_annot[gene] = final_annot
    gene_source_map[gene] = source

    if status != 'Consistent':
        comparison_records.append({
            'Contig': contig,
            'Gene': gene,
            'Pfam_Annot': p_annot,
            'Pharokka_Annot': ph_annot,
            'Discrepancy_Type': status,
            'Adopted_Annot': final_annot
        })

df_comparison = pd.DataFrame(comparison_records)

isescan_valid_hits = set()
try:
    if os.path.exists(ISESCAN_RAW_FILE):
        df_ise = pd.read_csv(ISESCAN_RAW_FILE, sep='\t', header=None, comment='#', usecols=range(10))
        subset_ise = df_ise[ (df_ise[6] <= 1e-9) & (df_ise[7] >= 60) ]
        isescan_valid_hits = set(subset_ise[0].values)
except Exception as e:
    print(f"Warning: Error processing ISEScan file: {e}")

def map_final_roles(row):
    gene = row['gene']
    annot = row['annot_lower']
    final_role = reconciled_gene_annot.get(gene, 'None')
    source = gene_source_map.get(gene, 'None')
    if final_role == 'None' and 'hypothetical' in annot and gene in isescan_valid_hits:
        final_role = 'transposase'
        source = 'ISEScan'

    is_tnp = (final_role in ['transposase', 'both'])
    is_int = (final_role in ['integrase', 'both'])
    tnp_source = source if is_tnp else 'None'

    return pd.Series([is_tnp, is_int, tnp_source])

df_cds[['Is_Transposase', 'Is_Integrase', 'Tnp_Source']] = df_cds.apply(map_final_roles, axis=1)

cog_target_reps = set() 
try:
    if os.path.exists(COG_ANNOT_FILE):
        with open(COG_ANNOT_FILE, 'r') as f:
            for line in f:
                if not line.strip(): continue
                parts = line.split('\t')
                if len(parts) >= 1:
                    cog_target_reps.add(parts[0].strip())
          else:
        print(f" no file {COG_ANNOT_FILE}")
except Exception as e:
    print(f"Error reading COG file: {e}")


mu_tnpb_genes = set() try:
    if os.path.exists(PC_MAP_FILE) and len(cog_target_reps) > 0:
       df_map = pd.read_csv(PC_MAP_FILE, sep='\t', low_memory=False)
        matched_rows = df_map[df_map['rep_seq_id'].isin(cog_target_reps)]
        mu_tnpb_genes = set(matched_rows['protein_id'].values)
    elif not os.path.exists(PC_MAP_FILE):
        print(f"no  {PC_MAP_FILE}")  
except Exception as e:
    print(f"Error reading Map file: {e}")
df_cds['Is_Mu_TnpB'] = df_cds['gene'].isin(mu_tnpb_genes)

genome_results = {}
for contig, group in df_cds.groupby('contig'):
    group = group.sort_values('phys_start').reset_index(drop=True)
    n_genes = len(group)
    has_tnp = group['Is_Transposase'].any()
    has_int = group['Is_Integrase'].any()
    if has_int and not has_tnp:
        basic_lifestyle = 'Integrase_lifestyle'
    elif has_tnp and not has_int:
        basic_lifestyle = 'Transposase_lifestyle'
    elif has_int and has_tnp:
        basic_lifestyle = 'Mixed_lifestyle'
    else:
        basic_lifestyle = 'Unannotated'
    encodes_transposase = 'Yes' if has_tnp else 'No'
    tnp_indices = group.index[group['Is_Transposase']].tolist()
    mu_indices = group.index[group['Is_Mu_TnpB']].tolist()
    mu_lifestyle = "Unknown"
    tnp_gene_id = "NA"
    tnp_source = "NA"
    mu_gene_id = "NA"
    found_pair = False
    for idx in tnp_indices:
        current_strand = group.loc[idx, 'frame']
        if idx > 0:
            prev_idx = idx - 1
            if group.loc[prev_idx, 'Is_Mu_TnpB']:
                if group.loc[prev_idx, 'frame'] == current_strand:
                    found_pair = True
                    tnp_gene_id = group.loc[idx, 'gene']
                    tnp_source = group.loc[idx, 'Tnp_Source']
                    mu_gene_id = group.loc[prev_idx, 'gene']
                    break
        if idx < n_genes - 1:
            next_idx = idx + 1
            if group.loc[next_idx, 'Is_Mu_TnpB']:
                if group.loc[next_idx, 'frame'] == current_strand:
                    found_pair = True
                    tnp_gene_id = group.loc[idx, 'gene']
                    tnp_source = group.loc[idx, 'Tnp_Source']
                    mu_gene_id = group.loc[next_idx, 'gene']
                    break

    if found_pair:
        mu_lifestyle = "Mu-like Lifestyle"
    elif len(tnp_indices) > 0 and len(mu_indices) == 0:
        mu_lifestyle = "Transposon Insertion / Other Tnp"
        tnp_gene_id = group.loc[tnp_indices[0], 'gene']
        tnp_source = group.loc[tnp_indices[0], 'Tnp_Source']
    elif len(tnp_indices) == 0 and len(mu_indices) > 0:
        mu_lifestyle = "Orphan Mu_TnpB"
        mu_gene_id = group.loc[mu_indices[0], 'gene']
    elif len(tnp_indices) > 0 and len(mu_indices) > 0:
        mu_lifestyle = "Transposase + Mu_TnpB (Separated)"
        tnp_gene_id = group.loc[tnp_indices[0], 'gene']
        mu_gene_id = group.loc[mu_indices[0], 'gene']

    genome_results[contig] = {
        'Contig': contig,
        'Basic_Lifestyle': basic_lifestyle,
        'Encodes_Transposase': encodes_transposase,
        'Analyzed_Mu_Lifestyle': mu_lifestyle,
        'Tnp_Gene': tnp_gene_id,
        'Tnp_Source': tnp_source,
        'Mu_TnpB_Gene': mu_gene_id
    }

df_res = pd.DataFrame.from_dict(genome_results, orient='index')

df_meta = pd.read_csv(METADATA_FILE, sep='\t', low_memory=False)
df_res['Merge_Key'] = df_res['Contig'].str.replace(r'_(\d+)-(\d+)$', r'_\1_\2', regex=True)
if 'Phage_Name' in df_meta.columns:
    df_meta['Phage_Name'] = df_meta['Phage_Name'].astype(str).str.strip()
    df_final = df_res.merge(df_meta, left_on='Merge_Key', right_on='Phage_Name', how='left')
    unmatched = df_final['Phage_Name'].isna().sum()
    if unmatched > 0:
        print(f"find {unmatched}  Contig not in Metadata ")
else:
    df_final = df_res

cols_to_export = [
    'Contig', 'Basic_Lifestyle', 'Encodes_Transposase', 'Analyzed_Mu_Lifestyle',
    'Tnp_Gene', 'Tnp_Source', 'Mu_TnpB_Gene',
    'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species',
    'vFAM', 'vSUBFAM', 'vGENUS', 'vOTU', 'Activity_Score'
]
cols_final = [c for c in cols_to_export if c in df_final.columns]

df_final[cols_final].to_csv(OUTPUT_FILE, sep='\t', index=False)

# ==============================================================================
print("\n=======================================================")
print("Summary of Basic Lifestyles")
print("=======================================================")
print(df_final['Basic_Lifestyle'].value_counts().to_string())

print("\n=======================================================")
print("Summary of Mu-like Lifestyles")
print("=======================================================")
print(df_final['Analyzed_Mu_Lifestyle'].value_counts().to_string())

