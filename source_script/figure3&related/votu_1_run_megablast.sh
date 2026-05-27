makeblastdb -in viruse.fasta -dbtype nucl -out virus_db

blastn \
    -query viruse.fasta \
    -db virus_db \
    -task megablast \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 2000 \
    -num_threads 24 \
    -out blast_output.tsv
