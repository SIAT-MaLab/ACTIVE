#!/bin/bash


COORD_FILE="hmp2_healthy_b1_coordinate.tsv"  # prophage coordinate
SEQ_FILE="all_hqpp_extracted.fna"            # prophage seq
OUTPUT_FILE="virome_pp.fna"                  
TEMP_LIST="target_ids.list"                  


RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'
ls *_virome_*.fastq.gz | sed 's/_virome.*//' | sort | uniq > current_samples.txt
SAMPLE_COUNT=$(wc -l < current_samples.txt)
awk 'NR==FNR {valid_samples[$1]=1; next} 
{
    split($1, parts, "-");
    sample_name = parts[1];
    if (sample_name in valid_samples) {
        printf "%s_%s_%s\n", $1, $2, $3
    }
}' current_samples.txt "$COORD_FILE" > "$TEMP_LIST"

ID_COUNT=$(wc -l < "$TEMP_LIST")

seqkit grep -f "$TEMP_LIST" "$SEQ_FILE" -o "$OUTPUT_FILE"


rm current_samples.txt "$TEMP_LIST"
