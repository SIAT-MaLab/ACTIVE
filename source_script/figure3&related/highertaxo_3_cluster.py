import pandas as pd
import subprocess
import os
import sys

# ===========================================
INPUT_ABC = "genome_similarity.abc"
OUTPUT_DIR = "clustering_results"

RANKS = [
    ('vFAM',    0.055),
    ('vSUBFAM', 0.32),
    ('vGENUS',  0.65),
    ('vSUBGEN', 0.80)
]

MCL_BIN = "mcl"
# ===========================================

def check_mcl_installed():
    try:
        subprocess.run([MCL_BIN, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("Error: 'mcl' command not found. Please install MCL first (sudo apt install mcl or conda install mcl)")
        sys.exit(1)

def main():
    check_mcl_installed()
    
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"1. Reading edge file: {INPUT_ABC} ...")
    
    edges_df = pd.read_csv(INPUT_ABC, sep='\t', names=['source', 'target', 'weight'])
    
    all_genomes = set(edges_df['source']).union(set(edges_df['target']))
    print(f"   - Number of genomes included: {len(all_genomes)}")
    print(f"   - Number of edges included: {len(edges_df)}")

    taxonomy_df = pd.DataFrame(index=list(all_genomes))
    taxonomy_df.index.name = 'genome_id'

    prev_rank_clusters = None 

    for rank_name, threshold in RANKS:
        print(f"\n>>> Processing rank: {rank_name} (Threshold >= {threshold})")
        
        current_edges = edges_df[edges_df['weight'] >= threshold].copy()
        
        if prev_rank_clusters is not None:
            current_edges['c1'] = current_edges['source'].map(prev_rank_clusters)
            current_edges['c2'] = current_edges['target'].map(prev_rank_clusters)
            
            before_prune = len(current_edges)
            current_edges = current_edges[
                (current_edges['c1'].notna()) & 
                (current_edges['c2'].notna()) & 
                (current_edges['c1'] == current_edges['c2'])
            ]
            after_prune = len(current_edges)
            print(f"   - Taxonomic Pruning: Removed {before_prune - after_prune} edges inconsistent with the previous rank.")
        
        print(f"   - Effective edges in current rank: {len(current_edges)}")

        if current_edges.empty:
            print(f"   Warning: No edges meet the criteria for rank {rank_name}. All genomes will be processed as Singletons.")
            taxonomy_df[rank_name] = [f"{rank_name}_S_{i}" for i in range(len(taxonomy_df))]
            prev_rank_clusters = taxonomy_df[rank_name].to_dict()
            continue

        temp_abc = os.path.join(OUTPUT_DIR, f"temp_{rank_name}.abc")
        current_edges[['source', 'target', 'weight']].to_csv(temp_abc, sep='\t', header=False, index=False)
        
        mcl_out = os.path.join(OUTPUT_DIR, f"{rank_name}.mcl_out")
        
        cmd = [MCL_BIN, temp_abc, "--abc", "-I", "2.0", "-o", mcl_out]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL) 
        
        current_rank_map = {} 
        
        with open(mcl_out, 'r') as f:
            for cluster_idx, line in enumerate(f):
                members = line.strip().split('\t')
                cluster_id = f"{rank_name}_{cluster_idx + 1}"
                
                for genome in members:
                    current_rank_map[genome] = cluster_id
        
        taxonomy_df[rank_name] = taxonomy_df.index.map(current_rank_map)
        
        null_mask = taxonomy_df[rank_name].isna()
        singleton_count = null_mask.sum()
        if singleton_count > 0:
            print(f"   - Found {singleton_count} Singletons (nodes with no connections).")
            singleton_ids = [f"{rank_name}_S_{i+1}" for i in range(singleton_count)]
            taxonomy_df.loc[null_mask, rank_name] = singleton_ids
            
            for genome, s_id in zip(taxonomy_df.index[null_mask], singleton_ids):
                current_rank_map[genome] = s_id

        prev_rank_clusters = current_rank_map
        
        if os.path.exists(temp_abc):
            os.remove(temp_abc)


    final_output = os.path.join(OUTPUT_DIR, "viral_taxonomy.csv")
    taxonomy_df.to_csv(final_output)
    print(f"\nAll completed! Taxonomy results saved to: {final_output}")
    print("Preview of first 5 rows:")
    print(taxonomy_df.head())

if __name__ == "__main__":
    main()
