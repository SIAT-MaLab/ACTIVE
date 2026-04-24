import pandas as pd
import os

metadata_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/active_tax/hankyphage/specific_species/2_taxonomy/vFAM1/vFAM1_metadata.tsv"
summary_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/votu_marker_tree/vOTU_selection_summary.tsv"
faa_file = "../vFAM1_prophages_proteins.faa"
output_faa = "vFAM1_All_vOTU_TerL_sequences.faa" 
df_meta = pd.read_csv(metadata_file, sep='\t')
all_votus = df_meta['vOTU'].unique()  
total_votu_count = len(all_votus)

df_summary = pd.read_csv(summary_file, sep='\t')
df_votu_status = df_summary[df_summary['vOTU'].isin(all_votus)].copy()
def normalize_gene_id(gene_id):
    if pd.isna(gene_id):
        return None
    if "_CDS_" in gene_id:
        prefix, suffix = gene_id.split("_CDS_")
        return f"{prefix}_{int(suffix)}"
    return gene_id

df_votu_status['Normalized_ID'] = df_votu_status['Selected_Gene'].apply(normalize_gene_id)

terl_positive = df_votu_status[df_votu_status['Status'] == 'Selected']
pos_count = len(terl_positive)
neg_count = total_votu_count - pos_count
target_ids = set(terl_positive['Normalized_ID'].dropna().unique())
extracted_count = 0
with open(faa_file, 'r') as f_in, open(output_faa, 'w') as f_out:
    can_write = False
    for line in f_in:
        if line.startswith('>'):
            current_id = line.split()[0].lstrip('>')
            if current_id in target_ids:
                can_write = True
                f_out.write(line)
                extracted_count += 1
            else:
                can_write = False
        else:
            if can_write:
                f_out.write(line)


# 6. 保存详细表格
terl_positive.drop(columns=['Normalized_ID']).to_csv("TerL_All_vOTU_Table.tsv", sep='\t', index=False)
print("详细清单已保存至: TerL_All_vOTU_Table.tsv")
