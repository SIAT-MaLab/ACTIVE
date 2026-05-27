import pandas as pd
import networkx as nx
import subprocess
import tempfile
import os
import sys
import numpy as np

ANI_FILE = "result_ani.tsv"
QUALITY_FILE = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_checkv_quality/quality_summary.tsv"
OUTPUT_FILE = "vOTU_clusters_uhgv_final.tsv"

THREADS = 20    
MCL_BINARY = "mcl"
MCL_INFLATION = 2.0
MIN_ANI = 95.0
MIN_AF = 85.0

print(f"--- Step 1: Loading Quality Data & Defining Tiers ---")
try:
    use_cols = ['contig_id', 'contig_length', 'completeness', 'checkv_quality', 'viral_genes']
    qdf = pd.read_csv(QUALITY_FILE, sep='\t')
    qdf.set_index('contig_id', inplace=True)
except Exception as e:
    print(f"Error reading quality file: {e}")
    sys.exit(1)

qdf['checkv_quality'] = qdf['checkv_quality'].fillna('')
is_complete = qdf['checkv_quality'] == 'Complete'
is_high_quality = (~is_complete) & (qdf['completeness'] > 90.0)

tiers = {
    1: set(qdf[is_complete].index),
    2: set(qdf[is_high_quality].index),
    3: set()
}
tiers[3] = set(qdf.index) - tiers[1] - tiers[2]

print(f"  Tier 1 (Complete):       {len(tiers[1])}")
print(f"  Tier 2 (High Quality):   {len(tiers[2])}")
print(f"  Tier 3 (Fragments):      {len(tiers[3])}")

print(f"\n--- Step 2: Processing Edges (Vectorized) ---")

try:
    ani_df = pd.read_csv(ANI_FILE, sep='\t')
except Exception as e:
    print(f"Error reading ANI file: {e}")
    sys.exit(1)

ani_df = ani_df[ani_df['qname'] != ani_df['tname']]
len_map = qdf['contig_length'].to_dict()
valid_ids = set(qdf.index)
ani_df = ani_df[ani_df['qname'].isin(valid_ids) & ani_df['tname'].isin(valid_ids)].copy()

ani_df['q_len'] = ani_df['qname'].map(len_map)
ani_df['t_len'] = ani_df['tname'].map(len_map)

ani_df['af'] = np.where(ani_df['q_len'] < ani_df['t_len'], ani_df['qcov'], ani_df['tcov'])

mask = (ani_df['pid'] >= MIN_ANI) & (ani_df['af'] >= MIN_AF)
edges_df = ani_df[mask].copy()

edges_df['weight'] = (edges_df['pid'] * edges_df['af']) / 100.0
del ani_df

print("  Building NetworkX Graph for clustering...")
G = nx.Graph()
G.add_nodes_from(qdf.index)
G.add_weighted_edges_from(edges_df[['qname', 'tname', 'weight']].to_records(index=False))

print(f"  Graph ready: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

def run_mcl_parallel(node_list):
    if not node_list: return []

    subgraph = G.subgraph(node_list)
    if subgraph.number_of_edges() == 0:
        return [[n] for n in node_list]

    fd, tmp_abc = tempfile.mkstemp(suffix='.abc')
    os.close(fd)

    with open(tmp_abc, 'w') as f:
        for u, v, d in subgraph.edges(data=True):
            f.write(f"{u}\t{v}\t{d['weight']}\n")

    tmp_out = tmp_abc + ".out"

    cmd = [MCL_BINARY, tmp_abc, "--abc", "-I", str(MCL_INFLATION), "-te", str(THREADS), "-o", tmp_out]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    clusters = []
    clustered_nodes = set()
    if os.path.exists(tmp_out):
        with open(tmp_out, 'r') as f:
            for line in f:
                mems = line.strip().split('\t')
                clusters.append(mems)
                clustered_nodes.update(mems)
        os.remove(tmp_out)
    os.remove(tmp_abc)

    for n in node_list:
        if n not in clustered_nodes:
            clusters.append([n])

    return clusters

print(f"\n--- Step 3: Iterative Clustering (Matrix Accelerated) ---")

final_votu_map = {} # node -> vOTU_ID
votu_idx = 1

# --- Round 1: Tier 1 (De Novo) ---
print(f"  [Round 1] Clustering Tier 1 ({len(tiers[1])} nodes)...")
if tiers[1]:
    c_t1 = run_mcl_parallel(list(tiers[1]))
    for mems in c_t1:
        vid = f"vOTU_{votu_idx}"
        votu_idx += 1
        for m in mems: final_votu_map[m] = vid

def process_tier_vectorized(target_nodes_set, tier_name):
    global votu_idx
    print(f"  [Round {tier_name}] Processing {len(target_nodes_set)} nodes...")

    if not target_nodes_set: return

    existing_pool = set(final_votu_map.keys())

    assignable_edges = edges_df[edges_df['qname'].isin(target_nodes_set) & edges_df['tname'].isin(existing_pool)]
    assignable_clusters = run_mcl_parallel(list(assignable_edges['qname'].unique()))

    for cluster in assignable_clusters:
        vid = f"vOTU_{votu_idx}"
        votu_idx += 1
        for m in cluster:
            final_votu_map[m] = vid

    unassigned_nodes = list(target_nodes_set - set(final_votu_map.keys()))
    if unassigned_nodes:
        print(f"    -> Clustering {len(unassigned_nodes)} unassigned nodes de novo...")
        new_clusters = run_mcl_parallel(unassigned_nodes)
        for mems in new_clusters:
            vid = f"vOTU_{votu_idx}"
            votu_idx += 1
            for m in mems:
                final_votu_map[m] = vid

process_tier_vectorized(tiers[2], "2 - High Qual")
process_tier_vectorized(tiers[3], "3 - Fragments")

print(f"\n--- Step 4: Selecting Representatives ---")
votu_dict = {}
for n, v in final_votu_map.items():
    if v not in votu_dict: votu_dict[v] = []
    votu_dict[v].append(n)

qdf['expected_len'] = np.where(qdf['completeness'] > 0, qdf['contig_length'] / (qdf['completeness']/100), qdf['contig_length'])
qdf['is_complete'] = qdf['checkv_quality'] == 'Complete'
qdf_fast = qdf[['contig_length', 'is_complete', 'expected_len', 'viral_genes']].to_dict('index')

output_records = []

for vid, members in votu_dict.items():
    m_info = []
    lens = []
    for m in members:
        d = qdf_fast[m]
        m_info.append((m, d['contig_length'], d['is_complete'], d['expected_len'], d.get('viral_genes', 0)))
        lens.append(d['contig_length'])

    median_len = np.median(lens)
    rep = None

    cands_a = [x for x in m_info if x[2] and x[1] >= median_len]
    if cands_a:
        cands_a.sort(key=lambda x: x[1], reverse=True)
        rep = cands_a[0][0]
    else:
        m_info.sort(key=lambda x: (abs(x[1] - x[3]), -x[4]))
        rep = m_info[0][0]

    for m in members:
        tier = "3"
        if m in tiers[1]: tier = "1"
        elif m in tiers[2]: tier = "2"

        output_records.append(f"{vid}\t{rep}\t{m}\t{qdf_fast[m]['contig_length']}\t{qdf.at[m, 'checkv_quality']}\tTier{tier}\n")

print(f"--- Step 5: Writing to {OUTPUT_FILE} ---")
with open(OUTPUT_FILE, 'w') as f:
    f.write("vOTU_ID\tRepresentative\tMember_ID\tLength\tCheckV_Quality\tTier\n")
    f.writelines(output_records)

print(f"Done. Total vOTUs: {votu_idx-1}")

