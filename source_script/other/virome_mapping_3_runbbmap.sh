STATS_SUFFIX="_coverage.txt"
for dir in */; do
    sample="${dir%/}"
    ref_file="${sample}/${sample}_virome.fna"
    r1_file="${sample}/${sample}_virome_1.fastq.gz"
    r2_file="${sample}/${sample}_virome_2.fastq.gz"
    output_stats="${sample}/${sample}${STATS_SUFFIX}"

        bbmap.sh $MEM_FLAG \
            ref="$ref_file" \
            in="$r1_file" \
            in2="$r2_file" \
            out=/dev/null \
            covstats="$output_stats" \
            minid=0.90 \
            threads=28 \
        
        if [ $? -eq 0 ]; then
            echo -e "success $output_stats"
        else
            echo -e "error $sample"
        fi
done
