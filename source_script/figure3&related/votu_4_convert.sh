awk -F'\t' 'NR>1 {print $3 "\t" $1 "\t" $2}' vOTU_clusters_uhgv_final.tsv > genome_to_vOTU.tsv

