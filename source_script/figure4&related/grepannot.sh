
LIST="vFAM7.list"
ANNOT="vFAM7_Heliusvirales_1032_annot/pharokka_cds_final_merged_output.tsv"
OUTPUT="vFAM7_extracted_annotations.tsv"

awk -F'\t' 'NR==FNR{ids[$1]; next} {gsub(/\r/,"",$5)} (FNR==1 || $5 in ids)' "$LIST" "$ANNOT" > "$OUTPUT"
cut -f5 vFAM7_Heliusvirales_1032_annot/pharokka_cds_final_merged_output.tsv | grep -v "contig" | head -n 3

ID_SAMPLE=$(head -n 1 vFAM7.list)

grep -w "$ID_SAMPLE" vFAM7_Heliusvirales_1032_annot/pharokka_cds_final_merged_output.tsv | head -n 1

cut -f5 vFAM7_extracted_annotations.tsv | grep -v "contig" | sort | uniq > found_ids.tmp

comm -23 <(sort vFAM7.list) found_ids.tmp > missing_ids.txt

