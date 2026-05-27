import pandas as pd
import numpy as np

metadata_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/active_tax/hankyphage/specific_species/2_taxonomy/vFAM1/vFAM1_metadata.tsv"
cross_host_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/active_tax_new_20260511/Detailed_Cross_Host_vOTUs_GlobalFDR.csv"
hanky_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/active_tax/hankyphage/specific_species/1_identify/hanky_votu_summary.tsv"
output_file = "vOTU_tree_annotation_updated.tsv"

df = pd.read_csv(metadata_file, sep='\t')
df['Activity_Score'] = pd.to_numeric(df['Activity_Score'], errors='coerce').fillna(0)

df_hanky_raw = pd.read_csv(hanky_file, sep='\t')
hanky_counts = df_hanky_raw.groupby('vOTU').size().reset_index(name='Hanky_Genomes')
def get_mode(x):
    m = x.mode()
    return m.iloc[0] if not m.empty else "Unknown"

agg_funcs = {
    'Genome': 'count',
    'Activity_Score': lambda x: (x >= 0.7).sum(),
    'Family': get_mode,
    'Genus': get_mode,
    'vSUBFAM': 'first',  
    'vGENUS': 'first'    
}

votu_stats = df.groupby('vOTU').agg(agg_funcs).reset_index()
votu_stats.columns = ['vOTU', 'Total_Genomes', 'Active_Genomes', 'Dominant_Host_Family', 'Dominant_Host_Genus', 'vSUBFAM', 'vGENUS']
votu_stats = pd.merge(votu_stats, hanky_counts, on='vOTU', how='left')
votu_stats['Hanky_Genomes'] = votu_stats['Hanky_Genomes'].fillna(0)
votu_stats['Log10_Genome_Count'] = np.log10(votu_stats['Total_Genomes'] + 1)
votu_stats['Activity_Rate'] = votu_stats['Active_Genomes'] / votu_stats['Total_Genomes']
votu_stats['Hanky_Rate'] = votu_stats['Hanky_Genomes'] / votu_stats['Total_Genomes']
votu_stats['Is_Hankyphage'] = votu_stats['Hanky_Genomes'].apply(lambda x: 'Yes' if x > 0 else 'No')

df_cross = pd.read_csv(cross_host_file)
df_cross_sub = df_cross[['vOTU', 'Phylum', 'Cross_Category', 'Significance_Type', 'Analysis_Group']].copy()
df_cross_sub['Is_SAvOTU'] = df_cross_sub['Analysis_Group'].apply(lambda x: 'Yes' if x == 'SAvOTU' else 'No')

final_table = pd.merge(votu_stats, df_cross_sub[['vOTU', 'Phylum', 'Cross_Category', 'Significance_Type', 'Is_SAvOTU']], on='vOTU', how='left')

final_table['Cross_Category'] = final_table['Cross_Category'].fillna('Single-Host')
final_table['Is_SAvOTU'] = final_table['Is_SAvOTU'].fillna('No')
final_table['Phylum'] = final_table['Phylum'].fillna('Unknown') 
final_table['Significance_Type'] = final_table['Significance_Type'].fillna('N/A') 

cols = ['vOTU', 'vSUBFAM', 'vGENUS', 'Total_Genomes', 'Log10_Genome_Count', 'Active_Genomes', 'Activity_Rate',
        'Hanky_Genomes', 'Hanky_Rate', 'Is_Hankyphage',
        'Phylum', 'Dominant_Host_Family', 'Dominant_Host_Genus',
        'Cross_Category', 'Significance_Type', 'Is_SAvOTU']
final_table = final_table[cols]

final_table.to_csv(output_file, sep='\t', index=False)

print(f"Success! Updated annotation file saved to: {output_file}")
print(final_table.head())
