
PROT_FILE="/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_annot/allpp_proteins.faa"

for hmm_file in *.hmm; do
    base_name=$(basename "$hmm_file" .hmm)
    hmmsearch --cpu 24 --noali -E 1e-5 --domtblout "${base_name}.domtblout" "$hmm_file" "$PROT_FILE" > /dev/null
    grep -v "^#" "${base_name}.domtblout" | \
    sort -k1,1 -k7,7g | \
    awk '!seen[$1]++' > "${base_name}.best_hit.tsv"
done

awk '{
    if ($0 == "" || /^==>/) next; 
    
    query = $1;
    evalue = $7 + 0; 
    
      if (!(query in min_e) || evalue < min_e[query]) {
        min_e[query] = evalue;
        best_line[query] = $0;
    }
} 
END {

    for (q in best_line) {
        print best_line[q];
    }
}' PF*.best_hit.tsv > all_merged_best_hits.tsv
