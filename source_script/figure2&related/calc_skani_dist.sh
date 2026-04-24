
GENOME_DIR="/home/huangyan/experiment/SPA/3/gut_isolate/allgenome"
INPUT_FILE="prophage_stats_final.tsv"
OUTPUT_FILE="species_skani_distances.tsv"
TMP_DIR="skani_work_final"
THREADS=16
mkdir -p "$TMP_DIR"
rm -f "$TMP_DIR"/*.txt "$TMP_DIR"/*.map
echo -e "Genus\tSpecies_A\tSpecies_B\tGenomic_Dist\tANI" > "$OUTPUT_FILE"

echo "Step 1: Reading metadata and organizing genome lists..."
tail -n +2 "$INPUT_FILE" | cut -f1,10,11 | while IFS=$'\t' read -r gid genus species; do
    if [[ -z "$gid" || -z "$genus" || "$species" == "s__" || -z "$species" ]]; then
        continue
    fi
    fpath="${GENOME_DIR}/${gid}.fasta"
    if [[ -f "$fpath" ]]; then
        echo "$fpath" >> "${TMP_DIR}/${genus}.txt"
        echo -e "${fpath}\t${species}" >> "${TMP_DIR}/${genus}.map"
    fi

done

echo "Step 2: Calculating Genomic Distances with skani..."
total_genera=$(ls "$TMP_DIR"/*.txt 2>/dev/null | wc -l)
current_idx=0

for list_file in "$TMP_DIR"/*.txt; do
    [[ -e "$list_file" ]] || continue
    ((current_idx++))

    genus=$(basename "$list_file" .txt)
    count=$(wc -l < "$list_file")
    if [[ "$count" -lt 2 ]]; then
        continue
    fi

    echo "[$current_idx/$total_genera] Processing Genus: $genus ($count genomes)..."
    raw_out="${TMP_DIR}/${genus}.skani"

    skani dist -t "$THREADS" --ql "$list_file" --rl "$list_file" -o "$raw_out" --medium >/dev/null 2>&1

    if [[ ! -f "$raw_out" ]]; then
        echo "    [Warning] skani failed to generate output for $genus"
        continue
    fi


    map_file="${TMP_DIR}/${genus}.map"
    awk -v genus="$genus" '
    BEGIN { OFS="\t"; }

    NR==FNR {
        path = $1;
        species_name = substr($0, index($0, "\t") + 1);
        path_to_sp[path] = species_name;
        next;
    }

    {
        if ($1 == "Ref_file" || $1 == "Ref") next;

        ref_p = $1; qry_p = $2; ani = $3;
        
        sp_a = path_to_sp[ref_p];
        sp_b = path_to_sp[qry_p];

        if (sp_a != "" && sp_b != "" && sp_a < sp_b) {
            dist = 1 - (ani / 100);
            if (dist < 0) dist = 0;
            print genus, sp_a, sp_b, dist, ani;
        }
    }
    ' "$map_file" "$raw_out" >> "$OUTPUT_FILE"
done

echo "------------------------------------------------"
echo "Mission Accomplished!"
echo "Final results saved to: $OUTPUT_FILE"
echo "Total species pairs calculated: $(tail -n +2 "$OUTPUT_FILE" | wc -l)"
