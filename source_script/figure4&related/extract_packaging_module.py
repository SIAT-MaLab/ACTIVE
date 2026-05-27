import pandas as pd
import numpy as np
import re
ANNOT_FILE = "vFAM7_extracted_annotations.tsv"
META_FILE = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/merged_phage_stats_taxonomy.tsv"
OUTPUT_FILE = "packaging_module_structure_vFAM7.csv"

KEY_A = ['HNH endonuclease', 'HNH']
KEY_B = ['terminase small subunit']
KEY_C = ['terminase large subunit']
KEY_PORTAL = ['portal protein']

MAX_AB_DIST = 17  
PORTAL_WINDOW = 5 # ===========================================

def get_gene_info(df, start_idx, end_idx):
    if start_idx + 1 >= end_idx:
        return 0, "None", "None"
    sub_df = df.loc[start_idx + 1 : end_idx - 1]
    return len(sub_df), ";".join(sub_df['annot'].fillna("Unknown")), ";".join(sub_df['category'].fillna("Unknown"))

def analyze_genome(group):
    group = group.reset_index(drop=True)
    potential_modules = []
    a_idxs = group.index[group['annot'].str.contains('|'.join(KEY_A), case=False, regex=True, na=False)].tolist()
    b_idxs = group.index[group['annot'].str.contains('|'.join(KEY_B), case=False, regex=True, na=False)].tolist()
    c_idxs = group.index[group['annot'].str.contains('|'.join(KEY_C), case=False, regex=True, na=False)].tolist()
    portal_idxs = group.index[group['annot'].str.contains('|'.join(KEY_PORTAL), case=False, regex=True, na=False)].tolist()

    for b_idx in b_idxs:
        strand = group.loc[b_idx, 'frame']
        
        # --- 1. A (HNH) ---
        found_a = False
        final_a_idx = None
        
        if strand == '+':
            candidate_as = [a for a in a_idxs if a < b_idx and group.loc[a, 'frame'] == '+']
            a_idx = max(candidate_as) if candidate_as else None
        else:
            candidate_as = [a for a in a_idxs if a > b_idx and group.loc[a, 'frame'] == '-']
            a_idx = min(candidate_as) if candidate_as else None
  
        if a_idx is not None and abs(a_idx - b_idx) <= MAX_AB_DIST:
            found_a = True
            final_a_idx = a_idx

        # --- 2.  C (TerL) ---
        found_c = False
        final_c_idx = None
        
        if strand == '+':
            candidate_cs = [c for c in c_idxs if c > b_idx and group.loc[c, 'frame'] == '+']
            c_idx = min(candidate_cs) if candidate_cs else None
        else:
            candidate_cs = [c for c in c_idxs if c < b_idx and group.loc[c, 'frame'] == '-']
            c_idx = max(candidate_cs) if candidate_cs else None


        if c_idx is not None:
            has_portal = False
            if strand == '+':
                search_range = range(c_idx + 1, min(c_idx + 1 + PORTAL_WINDOW, len(group)))
            else:
                search_range = range(max(0, c_idx - PORTAL_WINDOW), c_idx)

            for p_idx in search_range:
                if p_idx in portal_idxs and group.loc[p_idx, 'frame'] == strand:
                    has_portal = True
                    break
            
            if has_portal:
                found_c = True
                final_c_idx = c_idx


        if not found_a and not found_c:
            continue
        dist_ab = abs(final_a_idx - b_idx) - 1 if found_a else 9999

        potential_modules.append({
            'a_idx': final_a_idx, 
            'b_idx': b_idx, 
            'c_idx': final_c_idx,
            'dist_ab': dist_ab,
            'strand': strand,
            'has_a': found_a,
            'has_c': found_c
        })

    if not potential_modules:
        return None

    def sort_key(mod):
        completeness = 2 if (mod['has_a'] and mod['has_c']) else 1
        return (-completeness, mod['dist_ab'])

    best_mod = sorted(potential_modules, key=sort_key)[0]
    a, b, c = best_mod['a_idx'], best_mod['b_idx'], best_mod['c_idx']
    gap_ab_count, gap_ab_annot, gap_ab_cat = 0, "NA", "NA"
    gap_bc_count, gap_bc_annot, gap_bc_cat = 0, "NA", "NA"
    gene_a_id = "NA"
    gene_c_id = "NA"
    if best_mod['has_a']:
        gene_a_id = group.loc[a, 'gene']
        idx_min_ab, idx_max_ab = (a, b) if a < b else (b, a)
        gap_ab_count, gap_ab_annot, gap_ab_cat = get_gene_info(group, idx_min_ab, idx_max_ab)
    if best_mod['has_c']:
        gene_c_id = group.loc[c, 'gene']
        idx_min_bc, idx_max_bc = (b, c) if b < c else (c, b)
        gap_bc_count, gap_bc_annot, gap_bc_cat = get_gene_info(group, idx_min_bc, idx_max_bc)

    return {
        'Genome': group.loc[0, 'contig'],
        'Module_Strand': best_mod['strand'],
        'Gene_A_ID': gene_a_id,
        'Gene_B_ID': group.loc[b, 'gene'],         'Gene_C_ID': gene_c_id,
        'Gap_HNH_TerS_Count': gap_ab_count if best_mod['has_a'] else np.nan,
        'Gap_TerS_TerL_Count': gap_bc_count if best_mod['has_c'] else np.nan,
        'Gap_HNH_TerS_Annot': gap_ab_annot,
        'Gap_HNH_TerS_Category': gap_ab_cat,
        'Gap_TerS_TerL_Annot': gap_bc_annot,
        'Gap_TerS_TerL_Category': gap_bc_cat
    }

def main():
    print("1. Loading and Analyzing...")
    df_annot = pd.read_csv(ANNOT_FILE, sep='\t').sort_values(['contig', 'start'])
    results = []
    for contig, group in df_annot.groupby('contig'):
        res = analyze_genome(group)
        if res: results.append(res)

    df_res = pd.DataFrame(results)

    print(f"2. Merging Metadata... (Found {len(df_res)} valid modules)")
    df_meta = pd.read_csv(META_FILE, sep='\t')

    def clean_id(f):
        base = str(f).replace('phage_depth_', '').replace('.txt', '')
        return re.sub(r'_([^_]+)$', r'-\1', base)

    df_meta['Join_Key'] = df_meta['File'].apply(clean_id)
    df_final = pd.merge(df_res, df_meta[['Join_Key', 'vFAM', 'vSUBFAM', 'vGENUS', 'vOTU']],
                        left_on='Genome', right_on='Join_Key', how='left').drop(columns=['Join_Key'])

    df_final.to_csv(OUTPUT_FILE, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
