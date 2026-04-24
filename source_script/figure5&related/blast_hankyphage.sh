blastn \
  -query BK010646.fna \
  -db ../votu_blast/virus_db \
  -evalue 1e-5 \
  -max_hsps 1 \
  -max_target_seqs 50000 \
  -outfmt 6 \
  > hankyphage_vs_allphage.blastn.tsv
awk '$4 >= 4283' hankyphage_vs_allphage.blastn.tsv \
  > hankyphage_like_hits.tsv

