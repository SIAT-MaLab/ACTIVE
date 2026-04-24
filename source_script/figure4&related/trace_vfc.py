import pandas as pd
import numpy as np
from collections import Counter
JOINT_TAX_FILE = "viral_with_vfc_taxonomy.csv"
VFC_META_FILE = "vfc_votu.tsv"
MAPPING_FILE = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/votu_blast/vOTU-Representative.tsv"
FISHER_STATS_FILE = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/active_tax/stats_fisher_vOTU.tsv"
OLD_TAXONOMY_FILE = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/higher_tax/clustering_results/viral_taxonomy.csv"
OUTPUT_FILE = "vfc_full_flow_breakdown.tsv"

def main():
    
    # --- A. VFC Data ---
    df_vfc_meta = pd.read_csv(VFC_META_FILE, sep='\t')
    vfc_id_map = pd.Series(df_vfc_meta.VFC.values, index=df_vfc_meta.votu_id.astype(str)).to_dict()
    vfc_ids_set = set(df_vfc_meta.votu_id.astype(str))

    # --- B. User Mapping ---
    df_map = pd.read_csv(MAPPING_FILE, sep='\t')
    rep_to_votu = pd.Series(df_map.vOTU_ID.values, index=df_map.Representative.astype(str).str.strip()).to_dict()
    
    df_stats = pd.read_csv(FISHER_STATS_FILE, sep='\t')
    votu_sig_map = pd.Series(df_stats.Significance.values, index=df_stats.vOTU).to_dict()

    # --- C. User Old Taxonomy ---
    df_old_tax = pd.read_csv(OLD_TAXONOMY_FILE)
    genome_to_old_vfam = pd.Series(
        df_old_tax.vFAM.values, 
        index=df_old_tax.genome_id.astype(str).str.strip()
    ).to_dict()

    votu_old_fam_map = {}
    for rep_id, votu_id in rep_to_votu.items():
        if rep_id in genome_to_old_vfam:
            votu_old_fam_map[votu_id] = genome_to_old_vfam[rep_id]

    # --- D. Process Joint Clustering ---
    df_joint = pd.read_csv(JOINT_TAX_FILE)
    
    processed_rows = []
    
    for _, row in df_joint.iterrows():
        genome_id = str(row['genome_id']).strip()
        new_vfam = row['vFAM']
        
        source_type = "Unknown"
        vfc_origin = None
        user_votu = None
        user_sig = None
        user_old_fam = None
        
        # Check VFC
        base_id = genome_id.split('_provirus')[0]
        if genome_id in vfc_ids_set:
            source_type = "VFC"
            vfc_origin = vfc_id_map.get(genome_id)
        elif base_id in vfc_ids_set:
            source_type = "VFC"
            vfc_origin = vfc_id_map.get(base_id)
        
        # Check User
        if source_type == "Unknown":
            user_votu = rep_to_votu.get(genome_id)
            if user_votu:
                source_type = "User"
                user_sig = votu_sig_map.get(user_votu, "Unknown")
                user_old_fam = votu_old_fam_map.get(user_votu, "Unknown")
        
        processed_rows.append({
            'new_vfam': new_vfam,
            'type': source_type,
            'vfc_origin': vfc_origin,
            'user_votu': user_votu,
            'user_sig': user_sig,
            'user_old_fam': user_old_fam
        })
        
    df_proc = pd.DataFrame(processed_rows)

    df_only_vfc = df_proc[df_proc['type'] == 'VFC']
    # Group by [VFC_Source, Target_vFAM]
    vfc_flow = df_only_vfc.groupby(['vfc_origin', 'new_vfam']).size().reset_index(name='count')
    vfc_totals = df_only_vfc.groupby('vfc_origin').size().to_dict()
    vfc_flow['total_in_vfc'] = vfc_flow['vfc_origin'].map(vfc_totals)
    vfc_flow['percentage'] = (vfc_flow['count'] / vfc_flow['total_in_vfc']) * 100
    vfc_flow = vfc_flow.sort_values(by=['vfc_origin', 'percentage'], ascending=[True, False])
    
    final_results = []
    
    for _, row in vfc_flow.iterrows():
        vfc_source = row['vfc_origin']
        target_vfam = row['new_vfam']
        count_vfc = row['count']
        total_vfc = row['total_in_vfc']
        pct = row['percentage']

        user_rows = df_proc[(df_proc['new_vfam'] == target_vfam) & (df_proc['type'] == 'User')]
        
        user_total = len(user_rows)
        sig_rows = user_rows[user_rows['user_sig'] == 'Significant']
        sig_count = len(sig_rows)
        sig_list = list(sig_rows['user_votu'].unique())
        old_fams = user_rows['user_old_fam'].dropna().tolist()
        old_fams = [f for f in old_fams if f != 'Unknown']
        if old_fams:
            old_fam_counts = Counter(old_fams).most_common(3)
            old_fam_str = "; ".join([f"{k}({v})" for k, v in old_fam_counts])
        else:
            old_fam_str = "None"
            
        final_results.append({
            'VFC_Source': vfc_source,         
            'Target_vFAM': target_vfam,       
            'Flow_Count': count_vfc,          
            'Total_In_VFC': total_vfc,       
            'Flow_Percentage': round(pct, 2),
            'User_vOTU_Total': user_total,    
            'User_Significant_Count': sig_count,
            'User_Significant_List': ",".join(sig_list) if sig_list else "None",
            'User_Original_Families': old_fam_str
        })
        
    df_final = pd.DataFrame(final_results)

    df_final.to_csv(OUTPUT_FILE, sep='\t', index=False)

    vfc2_rows = df_final[df_final['VFC_Source'] == 2.0]
    if not vfc2_rows.empty:
        print(vfc2_rows[['VFC_Source', 'Target_vFAM', 'Flow_Percentage', 'User_vOTU_Total']].to_string())
    else:
        print(df_final[['VFC_Source', 'Target_vFAM', 'Flow_Percentage']].head(10).to_string())

if __name__ == "__main__":
    main()
