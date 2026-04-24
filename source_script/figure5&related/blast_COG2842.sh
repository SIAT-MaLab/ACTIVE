rpsblast -query mmseqs_out_rep_seq.faa -db /home/huangyan/databaseHY/CDD/raw/COG2842_db -evalue 1e-5 -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" -out COG2842_result.tsv -num_threads 24
awk '$9 <= 1e-9 && $10 > 50' COG2842_result.tsv > COG2842_result.final.tsv
