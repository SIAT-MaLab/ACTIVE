import pandas as pd
import os
import numpy as np

species_dir = "./species"
stats_table_path = "prophage_stats_final.tsv"
output_table = "tree_annotation_log1p.tsv"

print("Reading the total statistics table...")
df = pd.read_csv(stats_table_path, sep='\t')

df['Species'] = df['Species'].str.strip()
df['Clean_Species'] = df['Species'].str.replace(" ", "_")

df['Is_Lysogen'] = df['Prophage_Count'] > 0
df['Has_Active_Prophage'] = (df['Is_Lysogen']) & (df['Active_Prophage_Count'] > 0)

print("Calculating statistical metrics for each species...")

agg_rules = {
    'Genome_ID': 'count',
    'Is_Lysogen': 'sum',
    'Has_Active_Prophage': 'sum',
    'Prophage_Count': 'mean',
    'Phylum': 'first', 'Class': 'first', 'Order': 'first', 'Family': 'first', 'Genus': 'first'
}

df_stats = df.groupby('Clean_Species').agg(agg_rules).reset_index()

df_stats = df_stats.rename(columns={
    'Genome_ID': 'Total_Genomes_N',
    'Is_Lysogen': 'Total_Lysogens_N',
    'Has_Active_Prophage': 'Active_Lysogens_N',
    'Prophage_Count': 'Avg_Prophage_Count'
})

df_stats['Active_Rate_Pct'] = (df_stats['Active_Lysogens_N'] / df_stats['Total_Lysogens_N']) * 100
df_stats['Active_Rate_Pct'] = df_stats['Active_Rate_Pct'].fillna(0) 

df_stats['Log10_Total_Genomes'] = np.log10(df_stats['Total_Genomes_N'] + 1)
df_stats['Log10_Total_Lysogens'] = np.log10(df_stats['Total_Lysogens_N'] + 1)

print("Generating the final table...")

target_species = []
for f in os.listdir(species_dir):
    if f.endswith(".fasta") or f.endswith(".fa"):
        target_species.append(os.path.splitext(f)[0])

df_targets = pd.DataFrame({'Clean_Species': target_species})

final_df = pd.merge(df_targets, df_stats, on='Clean_Species', how='left')

final_df.fillna({
    'Total_Genomes_N': 0, 
    'Total_Lysogens_N': 0, 
    'Log10_Total_Genomes': 0, 
    'Log10_Total_Lysogens': 0,
    'Active_Rate_Pct': 0,
    'Avg_Prophage_Count': 0
}, inplace=True)

output_cols = [
    'Clean_Species',
    'Log10_Total_Genomes',   
    'Log10_Total_Lysogens',  
    'Active_Rate_Pct',       
    'Avg_Prophage_Count',    
    'Phylum', 'Class', 'Order', 'Family', 'Genus'
]

final_output = final_df[output_cols].sort_values(by='Log10_Total_Genomes', ascending=False)

final_output.to_csv(output_table, sep='\t', index=False, float_format='%.2f')

print(f"Processing complete! Results saved to: {output_table}")
