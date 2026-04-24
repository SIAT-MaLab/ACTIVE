diamond makedb --in vOUT_genomad_annot/vOTU_annotate/vOTU_proteins.faa -d viral_db
diamond blastp \
    --db viral_db \
    --query vOUT_genomad_annot/vOTU_annotate/vOTU_proteins.faa \
    --out "all_vs_all.tsv" \
    --outfmt 6 qseqid sseqid bitscore \
    --very-sensitive \
    --max-target-seqs 100000 \
    --evalue 1e-3 \
    --threads 26
