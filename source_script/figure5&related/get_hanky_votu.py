import pandas as pd
import sys

blast_file = "hankyphage_like_hits.tsv"
metadata_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/merged_phage_stats_taxonomy.tsv" # 根据你的相对路径调整

print(f"Loading BLAST hits from {blast_file}...")
try:
    blast_hits = pd.read_csv(blast_file, sep='\t', header=None, usecols=[1], names=['hit_id'])
    target_ids = set(blast_hits['hit_id'].unique())
    print(f"Found {len(target_ids)} unique Hanky-like hits.")
except Exception as e:
    print(f"Error reading BLAST file: {e}")
    sys.exit(1)

print(f"Loading Metadata from {metadata_file}...")

try:
    meta = pd.read_csv(metadata_file, sep='\t')
except Exception as e:
    print(f"Error reading Metadata file: {e}")
    sys.exit(1)

print("Constructing IDs for linking...")
meta['constructed_id'] = (
    meta['Contig_ID'].astype(str) + '_' + 
    meta['Start'].astype(str) + '-' + 
    meta['Stop'].astype(str)
)

matched_df = meta[meta['constructed_id'].isin(target_ids)].copy()

if matched_df.empty:
    print("\n[Warning] No matches found! Please check if the ID formats in BLAST and Metadata are consistent.")
    sys.exit()

votu_counts = matched_df['vOTU'].value_counts()

print("-" * 40)
print("Analysis Result: Hanky-like Phage Taxonomy")
print("-" * 40)
print(f"Total matched hits: {len(matched_df)}")
print(f"Unique vOTUs found: {len(votu_counts)}")
print("\nvOTU Counts (Top 10):")
print(votu_counts.head(10))

output_file = "hanky_votu_summary.tsv"
cols_to_save = ['constructed_id', 'vOTU', 'vFAM', 'Taxon_Lineage', 'Host_Species']
available_cols = [c for c in cols_to_save if c in matched_df.columns]
if not available_cols:
    available_cols = matched_df.columns

matched_df.to_csv(output_file, sep='\t', index=False, columns=available_cols)
print(f"\nDetailed summary saved to: {output_file}")
