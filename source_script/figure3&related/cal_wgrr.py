import pandas as pd
import argparse
import os
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description="Calculate wGRR matrix from MMseqs2 output.")
    parser.add_argument("-i", "--input_tsv", required=True, help="Input MMseqs2 all-vs-all TSV file")
    parser.add_argument("-f", "--input_faa", required=True, help="Input raw FAA file to count total proteins")
    parser.add_argument("-o", "--output_csv", required=True, help="Output wGRR edge list CSV")
    parser.add_argument("--threshold", type=float, default=0.05, help="Minimum wGRR to keep (default: 0.05)")
    args = parser.parse_args()

    protein_counts = defaultdict(int)
    if not os.path.exists(args.input_faa):
        print(f"missing {args.input_faa}")
        return

    with open(args.input_faa, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()[1:]
                genome_id = header.split("_CDS_")[0]
                protein_counts[genome_id] += 1

    df = pd.read_csv(args.input_tsv, sep='\t', usecols=[0, 1, 2], names=["query", "target", "pident"])
    df = df[df['query'] != df['target']].copy()
    df['query_genome'] = df['query'].apply(lambda x: x.split("_CDS_")[0])
    df['target_genome'] = df['target'].apply(lambda x: x.split("_CDS_")[0])
    df = df[df['query_genome'] != df['target_genome']]
    df['pident_frac'] = df['pident'] / 100.0
    df.drop(columns=['query', 'target', 'pident'], inplace=True)
    wgrr_sum = df.groupby(['query_genome', 'target_genome'])['pident_frac'].max().reset_index()
    wgrr_data = wgrr_sum.groupby(['query_genome', 'target_genome'])['pident_frac'].sum().reset_index()
    def calc_wgrr(row):
        gen_A = row['query_genome']
        gen_B = row['target_genome']
        min_genes = min(protein_counts.get(gen_A, 0), protein_counts.get(gen_B, 0))
        
        if min_genes == 0:
            return 0.0
        return row['pident_frac'] / min_genes

    wgrr_data['wGRR'] = wgrr_data.apply(calc_wgrr, axis=1)

    wgrr_data['wGRR'] = wgrr_data['wGRR'].clip(upper=1.0)
    wgrr_filtered = wgrr_data[wgrr_data['wGRR'] >= args.threshold]
    wgrr_filtered.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
