#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
echo ">>> Script directory: $SCRIPT_DIR"

show_help() {
cat <<'EOF'
ACTIVE - Prophage Induction Activity Testing Pipeline (Standalone Edition)

This pipeline tests whether sequencing depth in prophage regions is significantly higher
than the bacterial host genome background. Supports Illumina short reads (BBMap) and
ONT / PacBio CLR / PacBio HiFi long reads (Minimap2), including hybrid datasets.
All Python, R statistical, and visualization routines are embedded with no external scripts required.

Genome Naming Convention (strictly required for competitive mode):
  fasta  : sample01_bin.12.fasta
  contig : sample01_bin.12_1
  Column 1 of the prophage coordinate file must exactly match the sequence ID in the FASTA
  header (token before the first whitespace). Contig names across bins of the same sample
  must be globally unique in competitive mode.

Usage:
  bash active.sh -g <genome_dir> -s <seq_dir> -c <phage_coords> -o <out_dir>
             [-t <threads>] [-r <r_threads>] [-i <index_tsv>] [-m] [--score <float>] [-q]
             [--mask <bp>] [--prodigal-meta] [--map-mode per-bin|competitive]
             [--cov <float>] [--cutoff <float>] [--conservative-off] [--iter <int>]
             [--plot] [--plot-all] [--plot-overview] [--buffer <bp>] [--min-contig <bp>]
             [--rmbam] [--rmdepth] [--rmprodigal] [--rmmgs] [--clean] [--rmall]

Required Arguments:
  -g <dir>    Reference host genome directory; files must end with .fasta
  -s <dir>    Sequencing reads directory:
                Illumina PE : <sample>_1_clean.fastq[.gz] + <sample>_2_clean.fastq[.gz] (or .fq[.gz])
                Long reads  : <sample>_ont.fastq[.gz]        (ONT,         minimap2 -ax map-ont)
                              <sample>_pacbio.fastq[.gz]     (PacBio CLR,  minimap2 -ax map-pb)
                              <sample>_pbhifi.fastq[.gz]     (PacBio HiFi, minimap2 -ax map-hifi)
                Hybrid      : When <sample> has both Illumina PE and exactly one long-read type,
                              both are aligned independently and merged into a single BAM
  -c <file>   Prophage coordinate file (headerless tab-delimited TSV: contig  start  end, 1-based closed)
  -o <dir>    Output root directory

Core Alignment & Preprocessing Options:
  -t <int>    Alignment and preprocessing threads [default: 20]
  -r <int>    R statistical parallel threads [default: 20]
  -i <file>   Sample-to-reference index table (TSV with header: sample<TAB>ref)
              sample = read prefix without suffix; ref = genome filename without .fasta
              One sample may map to multiple rows (multiple bins)
  -m          Switch to precise mode [default: fast mode]
                fast    : Runs Prodigal + FetchMGs; uses single-copy marker genes (MGS) as host background
                precise : Skips gene calling; uses whole-genome non-prophage regions as host background
  --score <float>
              FetchMGs bit-score threshold (effective in fast mode only) [default: 300]
  -q          Enable samtools depth -Q 20 to filter low-MAPQ alignments [default: off]
  --mask <int>
              Mask n bp at both ends of every contig from depth extraction and statistics
              [default: 0, recommended: 100]. Contigs with length <= 2*n are excluded.
  --prodigal-meta
              Run Prodigal in -p meta mode (effective in fast mode only) [default: single]
  --map-mode <per-bin|competitive>
              Alignment mode [default: per-bin]
                per-bin     : Align reads to each reference genome independently
                competitive : Concatenate all bins linked to the same sample so reads compete once;
                              downstream depth extraction and statistics remain per-genome

Statistical & Activity Prediction Options:
  --cov <float>
              Minimum effective prophage breadth of coverage required for hypothesis testing
              [default: 0.75]. Regions below this threshold are assigned prediction = low_coverage
  --cutoff <float>
              Cliff's Delta effect-size threshold for calling a prophage active [default: 0.8]
              Recommended range: 0.8 - 0.95. Lower (more permissive) values detect weaker induction
              signals; higher (stricter) values increase false negatives but isolate stronger signals
  --conservative-off
              Disable conservative non-parametric-only mode [default: off].
              By default, Jarque-Bera normality testing is skipped and Brunner-Munzel + Cliff's Delta
              are applied to all testable regions. When --conservative-off is specified, Jarque-Bera
              normality testing is performed first; normally distributed pairs undergo Student's/Welch's
              t-test + Cohen's D (assigned active if p < 0.01 and Cohen's D >= 0.7, else dormant)
  --iter <int>
              Subsampling iterations for Cliff's Delta confidence intervals in precise mode [default: 500]

Visualization Options:
  --plot      Generate PDF plots for prophages predicted as active
                fast mode    : Outputs MGS vs Phage depth density plots only
                precise mode : Outputs local flanking window + density combined plots by default
  --plot-all  Generate PDF plots for all prophages (active, dormant, and low_coverage)
  --plot-overview
              (precise mode only) Additionally generate full-genome and host-contig overview depth
              profiles (100 bp binned). Automatically enables --plot if --plot-all is not set
  --buffer <int>
              (precise mode only) Host flanking length (bp) shown on each side of a prophage in
              local window plots [default: 30000]
  --min-contig <int>
              (--plot-overview only) Minimum contig length (bp) included in full-genome overview
              plots [default: 8000]

Intermediate File Cleanup Options:
  --rmbam      Remove BAM, BAI, and competitive concatenated FASTA files right after depth extraction
  --rmdepth    Remove per-base depth (.txt) files after statistics and plotting finish
  --rmprodigal Remove Prodigal output directory upon completion (fast mode only)
  --rmmgs      Remove FetchMGs and MGS coordinate intermediate files upon completion (fast mode only)
  --clean      Shorthand for --rmbam --rmdepth
  --rmall      Remove all intermediate files and directories, keeping only activity_results/
  -h, --help   Show this help message and exit
EOF
}

threads=20
r_threads=20
index_table=""
use_precise=false
mgs_score=300
use_mapq=false
prodigal_mode="single"
map_mode="per-bin"
mask_bp=0
coverage_threshold=0.75
cutoff=0.8
conservative_off=false
n_iter=500

do_plot=false
plot_all=false
plot_overview=false
plot_buffer=30000
plot_min_contig=8000

rm_bam=false
rm_depth=false
rm_prodigal=false
rm_mgs=false
rm_all=false

filtered_args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prodigal-meta)
            prodigal_mode="meta"; shift ;;
        --map-mode)
            map_mode="${2:-}"; shift 2 ;;
        --map-mode=*)
            map_mode="${1#*=}"; shift ;;
        --score)
            mgs_score="${2:-}"; shift 2 ;;
        --score=*)
            mgs_score="${1#*=}"; shift ;;
        --mask)
            mask_bp="${2:-}"; shift 2 ;;
        --mask=*)
            mask_bp="${1#*=}"; shift ;;
        --cov)
            coverage_threshold="${2:-}"; shift 2 ;;
        --cov=*)
            coverage_threshold="${1#*=}"; shift ;;
        --cutoff)
            cutoff="${2:-}"; shift 2 ;;
        --cutoff=*)
            cutoff="${1#*=}"; shift ;;
        --conservative-off)
            conservative_off=true; shift ;;
        --iter)
            n_iter="${2:-}"; shift 2 ;;
        --iter=*)
            n_iter="${1#*=}"; shift ;;
        --plot)
            do_plot=true; shift ;;
        --plot-all)
            do_plot=true; plot_all=true; shift ;;
        --plot-overview)
            plot_overview=true; shift ;;
        --buffer)
            plot_buffer="${2:-}"; shift 2 ;;
        --buffer=*)
            plot_buffer="${1#*=}"; shift ;;
        --min-contig)
            plot_min_contig="${2:-}"; shift 2 ;;
        --min-contig=*)
            plot_min_contig="${1#*=}"; shift ;;
        --rmbam)
            rm_bam=true; shift ;;
        --rmdepth)
            rm_depth=true; shift ;;
        --rmprodigal)
            rm_prodigal=true; shift ;;
        --rmmgs)
            rm_mgs=true; shift ;;
        --clean)
            rm_bam=true; rm_depth=true; shift ;;
        --rmall)
            rm_all=true; shift ;;
        --help)
            show_help; exit 0 ;;
        *)
            filtered_args+=("$1"); shift ;;
    esac
done

if [ ${#filtered_args[@]} -gt 0 ]; then
    set -- "${filtered_args[@]}"
else
    set --
fi

while getopts "hg:s:c:o:t:r:i:mq" opt; do
    case $opt in
        h) show_help; exit 0 ;;
        g) in_hostgenome="$OPTARG" ;;
        s) in_seq="$OPTARG" ;;
        c) phage_region="$OPTARG" ;;
        o) output="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        r) r_threads="$OPTARG" ;;
        i) index_table="$OPTARG" ;;
        m) use_precise=true ;;
        q) use_mapq=true ;;
        *) echo "Error: Invalid option. Use -h or --help for usage."; exit 1 ;;
    esac
done

if [ -z "${in_hostgenome:-}" ] || [ -z "${in_seq:-}" ] || [ -z "${phage_region:-}" ] || [ -z "${output:-}" ]; then
    echo "Error: Missing required arguments (-g, -s, -c, -o). Use -h for help."
    exit 1
fi

case "$map_mode" in
    per-bin|competitive) ;;
    *)
        echo "Error: --map-mode must be 'per-bin' or 'competitive' (got: $map_mode)"
        exit 1
        ;;
esac

if ! [[ "$mask_bp" =~ ^[0-9]+$ ]]; then
    echo "Error: --mask must be a non-negative integer (got: $mask_bp)"
    exit 1
fi

if ! [[ "$plot_buffer" =~ ^[0-9]+$ ]] || ! [[ "$plot_min_contig" =~ ^[0-9]+$ ]]; then
    echo "Error: --buffer and --min-contig must be non-negative integers"
    exit 1
fi

if ! awk -v x="$cutoff" -v c="$coverage_threshold" -v s="$mgs_score" 'BEGIN {
    if (x == "" || c == "" || s == "") exit 1
    if (x+0 != x || c+0 != c || s+0 != s) exit 1
    if (c < 0 || c > 1) exit 1
}'; then
    echo "Error: Invalid numeric parameter in --cutoff, --cov (must be 0..1), or --score"
    exit 1
fi

if [ ! -d "$in_hostgenome" ]; then
    echo "Error: Reference genome directory does not exist: $in_hostgenome"
    exit 1
fi
if [ ! -d "$in_seq" ]; then
    echo "Error: Sequencing reads directory does not exist: $in_seq"
    exit 1
fi
if [ ! -f "$phage_region" ]; then
    echo "Error: Prophage coordinate file does not exist: $phage_region"
    exit 1
fi

if [ "$use_precise" = true ]; then
    mode="precise"
else
    mode="fast"
fi

if [ "$plot_overview" = true ]; then
    if [ "$mode" = "precise" ]; then
        do_plot=true
    else
        echo ">>> [Notice] --plot-overview is only applicable to precise (-m) mode and will be ignored in fast mode."
        plot_overview=false
    fi
fi

if [ "$rm_all" = true ]; then
    rm_bam=true
    rm_depth=true
    rm_prodigal=true
    rm_mgs=true
fi

echo "=========================================="
echo "ACTIVE Configuration:"
echo "  Threads=$threads | R_Threads=$r_threads | Index=${index_table:-None}"
echo "  Mode=$mode | Map_Mode=$map_mode"
echo "  End_Mask=${mask_bp}bp | Min_Coverage=${coverage_threshold} | Cliff_Delta_Cutoff=${cutoff}"
echo "  Conservative_Off=${conservative_off} | Permutations=${n_iter}"
echo "  Prodigal=-p $prodigal_mode | FetchMGs_Score=${mgs_score} | MAPQ20_Filter=$use_mapq"
echo "  Plot=${do_plot} (Plot_All=${plot_all}, Overview=${plot_overview}, Buffer=${plot_buffer}bp, Min_Contig=${plot_min_contig}bp)"
echo "=========================================="

if [ "$mode" = "fast" ]; then
    echo ">>> Fast mode active: FetchMGs bit-score > $mgs_score, Prodigal = -p $prodigal_mode"
elif [ "$prodigal_mode" = "meta" ]; then
    echo ">>> [Notice] precise mode skips Prodigal; --prodigal-meta is ignored"
fi

if [ "$map_mode" = "competitive" ]; then
    echo ">>> Competitive mapping enabled: all bins per sample will be concatenated for alignment"
fi

if [ "$mask_bp" -gt 0 ]; then
    echo ">>> Contig end masking enabled: ${mask_bp} bp masked at each end (contigs <= $((mask_bp * 2)) bp excluded)"
fi

[ "$use_mapq" = true ] && echo ">>> Enabled samtools depth -Q 20"
[ "$rm_all" = true ] && echo ">>> Cleanup policy: --rmall (retaining only activity_results/)"
[ "$rm_bam" = true ] && [ "$rm_all" = false ] && echo ">>> Cleanup policy: --rmbam"
[ "$rm_depth" = true ] && [ "$rm_all" = false ] && echo ">>> Cleanup policy: --rmdepth"
[ "$rm_prodigal" = true ] && [ "$rm_all" = false ] && echo ">>> Cleanup policy: --rmprodigal"
[ "$rm_mgs" = true ] && [ "$rm_all" = false ] && echo ">>> Cleanup policy: --rmmgs"

mkdir -p "$output"

declare -A sample_to_ref
if [ -n "$index_table" ] && [ -f "$index_table" ]; then
    echo ">>> Loading sample-to-reference index table: $index_table"
    while IFS=$'\t' read -r sample ref _; do
        sample="${sample//$'\r'/}"
        ref="${ref//$'\r'/}"
        sample="${sample%"${sample##*[![:space:]]}"}"
        ref="${ref%"${ref##*[![:space:]]}"}"
        [ -z "$sample" ] && continue
        [ "$sample" = "sample" ] && continue
        if [ -z "${sample_to_ref["$sample"]:-}" ]; then
            sample_to_ref["$sample"]="$ref"
        else
            sample_to_ref["$sample"]="${sample_to_ref["$sample"]},$ref"
        fi
    done < <(tail -n +2 "$index_table")
else
    echo ">>> No index table provided; matching sample names directly to reference genome names"
fi

make_bed_from_tsv() {
    awk -v OFS='\t' 'NF>=3{ s=$2; e=$3; if(s<1)s=1; if(e>=s) print $1, s-1, e }' "$1" \
    | sort -k1,1 -k2,2n | bedtools merge -i - > "$2"
}

make_nonphage_bed() {
    bedtools subtract -a "$1" -b "$2" > "$3"
}

make_valid_contig_bed() {
    local fai_file="$1"
    local m="$2"
    local out_bed="$3"
    awk -v m="$m" -v OFS='\t' '{
        s = m
        e = $2 - m
        if (e > s) print $1, s, e
    }' "$fai_file" | sort -k1,1 -k2,2n > "$out_bed"
}

ensure_fai() {
    local ref_file="$1"
    if [ ! -f "${ref_file}.fai" ]; then
        samtools faidx "$ref_file"
    fi
}

genome_has_mgs() {
    local ref_file="$1"
    local mgs_tsv="$2"
    [ -s "$mgs_tsv" ] || return 1
    [ -f "$ref_file" ] || return 1
    ensure_fai "$ref_file"
    awk -v m="$mask_bp" '
        NR==FNR {
            if ($2 > 2 * m) {
                ms[$1] = m + 1
                me[$1] = $2 - m
            }
            next
        }
        ($1 in ms) && NF>=3 {
            s = ($2 > ms[$1] ? $2 : ms[$1])
            e = ($3 < me[$1] ? $3 : me[$1])
            if (e >= s) { found=1; exit }
        }
        END { exit !found }
    ' "${ref_file}.fai" "$mgs_tsv"
}

filter_regions_by_fai() {
    local fai_file="$1"
    local region_file="$2"
    local out_file="$3"
    awk 'NR==FNR { c[$1]=1; next } (NF>=3 && ($1 in c)) { sub(/\r$/, "", $3); print $1"\t"$2"\t"$3 }' "$fai_file" "$region_file" > "$out_file"
}

LR_FILE=""
LR_TYPE=""
LR_PRESET=""

find_longread_for_base() {
    local base="$1"
    local dir="$2"
    local -a files=()
    local -a types=()
    local -a presets=()
    local suf pre ext cand

    for spec in "pbhifi:map-hifi" "pacbio:map-pb" "ont:map-ont"; do
        suf="${spec%%:*}"
        pre="${spec##*:}"
        for ext in fastq.gz fastq fq.gz fq; do
            cand="${dir}/${base}_${suf}.${ext}"
            if [ -f "$cand" ]; then
                files+=("$cand")
                types+=("$suf")
                presets+=("$pre")
                break
            fi
        done
    done

    if [ ${#files[@]} -gt 1 ]; then
        echo "Error: Sample ${base} has multiple long-read types: ${types[*]}"
        echo "       Only one of _pbhifi / _pacbio / _ont is allowed per sample."
        exit 1
    fi

    if [ ${#files[@]} -eq 1 ]; then
        LR_FILE="${files[0]}"
        LR_TYPE="${types[0]}"
        LR_PRESET="${presets[0]}"
    else
        LR_FILE=""
        LR_TYPE=""
        LR_PRESET=""
    fi
}

has_illumina_r1() {
    local base="$1"
    local dir="$2"
    [ -f "${dir}/${base}_1_clean.fastq.gz" ] || [ -f "${dir}/${base}_1_clean.fastq" ] \
    || [ -f "${dir}/${base}_1_clean.fq.gz" ] || [ -f "${dir}/${base}_1_clean.fq" ]
}

resolve_ref_str() {
    local base_name="$1"
    if [ -n "$index_table" ] && [ ${#sample_to_ref[@]} -gt 0 ]; then
        echo "${sample_to_ref[$base_name]:-}"
    else
        echo "$base_name"
    fi
}

normalize_ref_str() {
    local ref_str="$1"
    local -A seen=()
    local -a refs=()
    local r
    local IFS=','
    local -a tmp=()
    read -ra tmp <<< "$ref_str"
    for r in "${tmp[@]}"; do
        r="${r// /}"
        [ -z "$r" ] && continue
        [ -n "${seen[$r]+x}" ] && continue
        seen[$r]=1
        refs+=("$r")
    done
    local IFS=','
    echo "${refs[*]}"
}

run_bbmap_pe() {
    local ref_file="$1"
    local data1="$2"
    local data2="$3"
    local out_bam="$4"

    echo "        -> BBMap PE minid=0.95 ambiguous=random : $(basename "$data1") vs $(basename "$ref_file")"
    set -o pipefail
    bbmap.sh in1="$data1" in2="$data2" ref="$ref_file" out=stdout.sam threads=$threads \
        ambiguous=random minid=0.95 nodisk=t \
    | samtools view -Sb -@ "$threads" - \
    | samtools sort -o "$out_bam" -@ "$threads" -
    local rc=$?
    set +o pipefail
    if [ $rc -ne 0 ] || [ ! -s "$out_bam" ]; then
        echo "        Error: BBMap/samtools failed (exit code $rc)"
        rm -f "$out_bam"
        return 1
    fi
    return 0
}

run_minimap2() {
    local preset="$1"
    local ref_file="$2"
    local reads="$3"
    local out_bam="$4"

    echo "        -> Minimap2 -ax ${preset} --secondary=no : $(basename "$reads") vs $(basename "$ref_file")"
    set -o pipefail
    minimap2 -ax "$preset" --secondary=no -t "$threads" "$ref_file" "$reads" \
    | samtools view -Sb -@ "$threads" - \
    | samtools sort -o "$out_bam" -@ "$threads" -
    local rc=$?
    set +o pipefail
    if [ $rc -ne 0 ] || [ ! -s "$out_bam" ]; then
        echo "        Error: Minimap2/samtools failed (exit code $rc)"
        rm -f "$out_bam"
        return 1
    fi
    return 0
}

need_recompute() {
    local target="$1"
    local source="$2"
    local force_flag="${3:-false}"
    [ "$force_flag" = true ] || [ ! -s "$target" ] || [ "$source" -nt "$target" ]
}

merge_sr_lr_bam() {
    local bam_file="$1"
    local temp_sr_bam="$2"
    local temp_lr_bam="$3"
    local has_sr_bam="$4"
    local has_lr_bam="$5"
    local tag="$6"

    if [ "$has_sr_bam" = true ] && [ "$has_lr_bam" = true ]; then
        echo "        -> Merging short-read and long-read BAMs ($tag)..."
        samtools merge -@ "$threads" -f "$bam_file" "$temp_sr_bam" "$temp_lr_bam"
        rm -f "$temp_sr_bam" "$temp_lr_bam"
    elif [ "$has_sr_bam" = true ]; then
        mv "$temp_sr_bam" "$bam_file"
        rm -f "$temp_lr_bam"
    elif [ "$has_lr_bam" = true ]; then
        mv "$temp_lr_bam" "$bam_file"
        rm -f "$temp_sr_bam"
    else
        rm -f "$temp_sr_bam" "$temp_lr_bam"
        return 1
    fi
    samtools index "$bam_file"
}

extract_depth_for_bam() {
    local bam_file="$1"
    local ref_file="$2"
    local sample_dir="$3"
    local ref_name="$4"

    local depth_opts=( -@ "$threads" -a )
    if [ "$use_mapq" = true ]; then
        depth_opts+=( -Q 20 )
    fi

    ensure_fai "$ref_file"

    local work_dir="${sample_dir}/global/${ref_name}"
    mkdir -p "$work_dir"

    local phage_region_sub="$work_dir/phage_region.tsv"
    filter_regions_by_fai "${ref_file}.fai" "$phage_region" "$phage_region_sub"

    if [ ! -s "$phage_region_sub" ]; then
        echo "    [Skip] No matching prophage coordinates found for $ref_name"
        return 0
    fi

    local genome_depth_dir="${sample_dir}/depth_to_stat/${ref_name}"
    mkdir -p "$genome_depth_dir"

    local param_stamp_file="$genome_depth_dir/.depth_params.stamp"
    local current_stamp="mask=${mask_bp};mapq=${use_mapq}"
    local force_depth=false
    if [ ! -f "$param_stamp_file" ] || [ "$(cat "$param_stamp_file" 2>/dev/null)" != "$current_stamp" ]; then
        force_depth=true
        rm -f "$genome_depth_dir"/phage_depth_*.txt \
              "$genome_depth_dir"/*_mgs_depth.txt \
              "$genome_depth_dir"/*_host_nonphage_depth.txt
    fi

    local valid_bed="$genome_depth_dir/${ref_name}.valid.bed"
    make_valid_contig_bed "${ref_file}.fai" "$mask_bp" "$valid_bed"

    if [ ! -s "$valid_bed" ]; then
        echo "    [Skip] $ref_name has no valid contig sequence remaining after ${mask_bp} bp end masking"
        return 0
    fi

    local raw_phage_bed="$genome_depth_dir/${ref_name}.phage.raw.bed"
    local masked_phage_bed="$genome_depth_dir/${ref_name}.phage.masked.bed"
    make_bed_from_tsv "$phage_region_sub" "$raw_phage_bed"
    bedtools intersect -a "$raw_phage_bed" -b "$valid_bed" > "$masked_phage_bed"

    local -A contig_len=()
    while IFS=$'\t' read -r cname clen _; do
        [ -n "$cname" ] && contig_len["$cname"]="$clen"
    done < "${ref_file}.fai"

    local n_phage_total=0
    local n_phage_kept=0
    local n_phage_masked_out=0

    while IFS=$'\t' read -r contig start end; do
        [ -z "$contig" ] && continue
        n_phage_total=$((n_phage_total + 1))

        local clen="${contig_len[$contig]:-0}"
        local min_pos=$((mask_bp + 1))
        local max_pos=$((clen - mask_bp))

        local eff_start="$start"
        local eff_end="$end"
        [ "$eff_start" -lt 1 ] && eff_start=1
        [ "$eff_end" -gt "$clen" ] && eff_end="$clen"
        [ "$eff_start" -lt "$min_pos" ] && eff_start="$min_pos"
        [ "$eff_end" -gt "$max_pos" ] && eff_end="$max_pos"

        local out_d="$genome_depth_dir/phage_depth_${contig}_${start}_${end}.txt"

        if [ "$eff_start" -gt "$eff_end" ]; then
            n_phage_masked_out=$((n_phage_masked_out + 1))
            rm -f "$out_d"
            echo "    [Mask-Skip] Prophage ${contig}:${start}-${end} falls completely within ${mask_bp} bp masked ends"
            continue
        fi

        n_phage_kept=$((n_phage_kept + 1))
        if need_recompute "$out_d" "$bam_file" "$force_depth"; then
            samtools depth "${depth_opts[@]}" -r "${contig}:${eff_start}-${eff_end}" "$bam_file" \
            | awk -v OFS='\t' '{print $1,$2,$3}' > "$out_d"
        fi
    done < "$phage_region_sub"

    echo "    [Depth] $ref_name : ${n_phage_total} total prophages, ${n_phage_kept} retained (${n_phage_masked_out} masked out)"

    if [ "$n_phage_kept" -eq 0 ]; then
        echo "    [Skip] $ref_name has no valid prophage regions remaining for statistical testing"
        echo "$current_stamp" > "$param_stamp_file"
        return 0
    fi

    if [ "$mode" = "fast" ]; then
        local mgs_region_sub="$work_dir/mgs_region.tsv"
        filter_regions_by_fai "${ref_file}.fai" "$output/mgs_ordinate.tsv" "$mgs_region_sub"
        local out_mgs="$genome_depth_dir/${ref_name}_mgs_depth.txt"

        if [ -s "$mgs_region_sub" ]; then
            local raw_mgs_bed="$genome_depth_dir/${ref_name}.mgs.raw.bed"
            local masked_mgs_bed="$genome_depth_dir/${ref_name}.mgs.bed"
            make_bed_from_tsv "$mgs_region_sub" "$raw_mgs_bed"
            bedtools intersect -a "$raw_mgs_bed" -b "$valid_bed" > "$masked_mgs_bed"

            if [ -s "$masked_mgs_bed" ]; then
                if need_recompute "$out_mgs" "$bam_file" "$force_depth"; then
                    samtools depth "${depth_opts[@]}" -b "$masked_mgs_bed" "$bam_file" \
                    | awk -v OFS='\t' '{print $1,$2,$3}' > "$out_mgs"
                fi
            else
                rm -f "$out_mgs"
                echo "    [Warning] MGS regions for $ref_name are empty after ${mask_bp} bp end masking"
            fi
        else
            rm -f "$out_mgs"
            echo "    [Warning] No matching MGS regions for $ref_name"
        fi
    else
        local nonphage_bed="$genome_depth_dir/${ref_name}.nonphage.bed"
        local out_host="$genome_depth_dir/${ref_name}_host_nonphage_depth.txt"
        make_nonphage_bed "$valid_bed" "$raw_phage_bed" "$nonphage_bed"

        if [ -s "$nonphage_bed" ]; then
            if need_recompute "$out_host" "$bam_file" "$force_depth"; then
                samtools depth "${depth_opts[@]}" -b "$nonphage_bed" "$bam_file" \
                | awk -v OFS='\t' '{print $1,$2,$3}' > "$out_host"
            fi
        else
            rm -f "$out_host"
            echo "    [Warning] No host non-prophage background remaining for $ref_name after subtracting prophages and masked ends"
        fi
    fi

    echo "$current_stamp" > "$param_stamp_file"
}

process_one_sample_ref() {
    local base_name="$1"
    local ref_name="$2"
    local data1="$3"
    local data2="$4"
    local lr_file="$5"
    local lr_preset="$6"

    local ref_file="$in_hostgenome/${ref_name}.fasta"
    if [ ! -f "$ref_file" ]; then
        echo "    [Skip] Reference genome file not found: $ref_file"
        return 0
    fi

    if [ "$mode" = "fast" ]; then
        if ! genome_has_mgs "$ref_file" "$output/mgs_ordinate.tsv"; then
            echo "    [Skip] fast mode: $ref_name has no valid MGS background; skipping alignment and testing"
            return 0
        fi
    fi

    local sample_dir="${output_map_dir}/${base_name}"
    mkdir -p "$sample_dir"
    local bam_file="${sample_dir}/${base_name}_${ref_name}.bam"

    if [ -s "$bam_file" ]; then
        echo "    [Reuse] Existing BAM found: $bam_file"
    else
        echo "    [Mapping-per-bin] $base_name -> $ref_name"

        local temp_sr_bam="${sample_dir}/temp_sr_${base_name}_${ref_name}.bam"
        local temp_lr_bam="${sample_dir}/temp_lr_${base_name}_${ref_name}.bam"
        local has_sr_bam=false
        local has_lr_bam=false

        if [ -n "$data1" ] && [ -n "$data2" ]; then
            if run_bbmap_pe "$ref_file" "$data1" "$data2" "$temp_sr_bam"; then
                has_sr_bam=true
            fi
        fi

        if [ -n "$lr_file" ]; then
            if run_minimap2 "$lr_preset" "$ref_file" "$lr_file" "$temp_lr_bam"; then
                has_lr_bam=true
            fi
        fi

        if ! merge_sr_lr_bam "$bam_file" "$temp_sr_bam" "$temp_lr_bam" "$has_sr_bam" "$has_lr_bam" "per-bin"; then
            echo "    Error: Failed to produce valid BAM for $base_name vs $ref_name"
            return 1
        fi
    fi

    extract_depth_for_bam "$bam_file" "$ref_file" "$sample_dir" "$ref_name"

    if [ "$rm_bam" = true ]; then
        echo "    [--rmbam] Removing BAM: $bam_file"
        rm -f "$bam_file" "${bam_file}.bai"
    fi
}

process_competitive_sample() {
    local base_name="$1"
    local ref_str="$2"
    local data1="$3"
    local data2="$4"
    local lr_file="$5"
    local lr_preset="$6"

    local concat_dir="$output/competitive_refs"
    mkdir -p "$concat_dir"
    local concat_fa="$concat_dir/${base_name}.fasta"
    local concat_list="$concat_dir/${base_name}.bins.txt"

    local -a valid_refs=()
    local -a valid_fas=()
    local ref_name ref_file

    local IFS=','
    local -a REF_GENOMES=()
    read -ra REF_GENOMES <<< "$ref_str"
    unset IFS

    for ref_name in "${REF_GENOMES[@]}"; do
        ref_name="${ref_name// /}"
        [ -z "$ref_name" ] && continue
        ref_file="$in_hostgenome/${ref_name}.fasta"
        if [ ! -f "$ref_file" ]; then
            echo "    [Skip] Reference genome file not found: $ref_file"
            continue
        fi
        valid_refs+=("$ref_name")
        valid_fas+=("$ref_file")
    done

    if [ ${#valid_fas[@]} -eq 0 ]; then
        echo "    [Skip] Sample $base_name has no valid reference FASTA files"
        return 0
    fi

    local sample_dir="${output_map_dir}/${base_name}"
    mkdir -p "$sample_dir"
    local bam_file="${sample_dir}/${base_name}_competitive.bam"
    local temp_sr_bam="${sample_dir}/temp_sr_${base_name}_competitive.bam"
    local temp_lr_bam="${sample_dir}/temp_lr_${base_name}_competitive.bam"

    local need_rebuild=false
    if [ ! -s "$concat_fa" ]; then
        need_rebuild=true
    else
        local fa
        for fa in "${valid_fas[@]}"; do
            if [ "$fa" -nt "$concat_fa" ]; then
                need_rebuild=true
                break
            fi
        done
        local expected
        expected=$(printf '%s\n' "${valid_refs[@]}" | sort)
        if [ ! -f "$concat_list" ] || [ "$expected" != "$(sort "$concat_list")" ]; then
            need_rebuild=true
        fi
    fi

    if [ "$need_rebuild" = true ]; then
        echo "    [Concat-Ref] $base_name : concatenating ${#valid_fas[@]} bins -> $(basename "$concat_fa")"
        local tmp_fa="${concat_fa}.tmp"
        cat "${valid_fas[@]}" > "$tmp_fa"
        if ! awk '/^>/ {
                h=$1
                sub(/^>/, "", h)
                if (h in seen) {
                    printf("Error: Duplicate contig ID in concatenated reference: %s\n", h) > "/dev/stderr"
                    exit 1
                }
                seen[h]=1
            }' "$tmp_fa"; then
            rm -f "$tmp_fa"
            echo "    Error: Duplicate contig names across bins of $base_name; cannot run competitive mapping"
            return 1
        fi
        mv -f "$tmp_fa" "$concat_fa"
        printf '%s\n' "${valid_refs[@]}" | sort > "$concat_list"
        samtools faidx "$concat_fa"
        rm -f "$bam_file" "${bam_file}.bai" "$temp_sr_bam" "$temp_lr_bam"
    else
        echo "    [Concat-Ref] Reusing concatenated reference: $concat_fa (${#valid_fas[@]} bins)"
        ensure_fai "$concat_fa"
    fi

    if [ -s "$bam_file" ]; then
        echo "    [Reuse] Existing competitive BAM found: $bam_file"
    else
        echo "    [Mapping-competitive] $base_name -> ${#valid_fas[@]} bins"

        local has_sr_bam=false
        local has_lr_bam=false

        if [ -n "$data1" ] && [ -n "$data2" ]; then
            if run_bbmap_pe "$concat_fa" "$data1" "$data2" "$temp_sr_bam"; then
                has_sr_bam=true
            fi
        fi

        if [ -n "$lr_file" ]; then
            if run_minimap2 "$lr_preset" "$concat_fa" "$lr_file" "$temp_lr_bam"; then
                has_lr_bam=true
            fi
        fi

        if ! merge_sr_lr_bam "$bam_file" "$temp_sr_bam" "$temp_lr_bam" "$has_sr_bam" "$has_lr_bam" "competitive"; then
            echo "    Error: Competitive mapping failed to produce BAM for $base_name"
            return 1
        fi
    fi

    for ref_name in "${valid_refs[@]}"; do
        ref_file="$in_hostgenome/${ref_name}.fasta"
        if [ "$mode" = "fast" ]; then
            if ! genome_has_mgs "$ref_file" "$output/mgs_ordinate.tsv"; then
                echo "    [Skip-Stat] fast mode: $ref_name has no valid MGS (participated in competitive alignment only)"
                continue
            fi
        fi
        extract_depth_for_bam "$bam_file" "$ref_file" "$sample_dir" "$ref_name"
    done

    if [ "$rm_bam" = true ]; then
        echo "    [--rmbam] Removing competitive BAM and concatenated reference: $bam_file"
        rm -f "$bam_file" "${bam_file}.bai" "$concat_fa" "${concat_fa}.fai" "$concat_list"
    fi
}

output_host_faa="$output/prodigal_results"
fetchmgs_dir="$output/fetchmgs_results"
genomes_with_mgs_file="$output/genomes_with_mgs.txt"

if [ "$mode" = "fast" ]; then
    echo ">>> [Step 1] Predicting host CDS with Prodigal (-p ${prodigal_mode}, max jobs ${threads})..."
    mkdir -p "$output_host_faa"

    running=0
    prodigal_fail=0
    for hostgenome in "$in_hostgenome"/*.fasta; do
        [ -e "$hostgenome" ] || continue
        host_name=$(basename "$hostgenome" .fasta)
        outdir="$output_host_faa/$host_name"
        mkdir -p "$outdir"

        target_gff="$outdir/$host_name.gff"
        target_faa="$outdir/$host_name.faa"

        if [ -s "$target_gff" ] && [ -s "$target_faa" ]; then
            echo "    [Skip] Prodigal output exists: $host_name"
            continue
        fi

        (
            echo "    [Run] Prodigal -p ${prodigal_mode}: $host_name"
            prodigal -i "$hostgenome" -o "$target_gff" -f gff -a "$target_faa" -p "$prodigal_mode" -q
        ) &

        running=$((running + 1))
        if [ "$running" -ge "$threads" ]; then
            if ! wait -n; then
                prodigal_fail=1
            fi
            running=$((running - 1))
        fi
    done

    while [ "$running" -gt 0 ]; do
        if ! wait -n; then
            prodigal_fail=1
        fi
        running=$((running - 1))
    done

    if [ "$prodigal_fail" -ne 0 ]; then
        echo "Error: Prodigal failed on one or more reference genomes. Aborting."
        exit 1
    fi

    echo ">>> [Step 2] Extracting single-copy marker genes with FetchMGs..."
    mkdir -p "$fetchmgs_dir"

    find "$output_host_faa" -name "*.faa" | while read -r hostfaa; do
        host_base=$(basename "${hostfaa%.faa}")
        target_scores_dir="$fetchmgs_dir/$host_base"
        target_scores_file="$fetchmgs_dir/$host_base/${host_base}.faa.fetchMGs.scores"

        if [ -s "$fetchmgs_dir/$host_base.faa.fetchMGs.scores" ]; then
             echo "    [Skip] FetchMGs scores exist: $host_base"
        elif [ -d "$target_scores_dir" ] && [ -s "$target_scores_file" ]; then
             echo "    [Skip] FetchMGs directory exists: $host_base"
        else
             echo "    [Run] FetchMGs: $host_base"
             fetchMGs extraction "$hostfaa" gene "$fetchmgs_dir/$host_base" -t "$threads" >/dev/null 2>&1
        fi
    done

    find "$fetchmgs_dir" -mindepth 2 -name "*fetchMGs.scores" -exec mv {} "$fetchmgs_dir/" \;

    echo "    Filtering qualifying COGs (bit-score > ${mgs_score})..."
    cog_list_file="$output/target_cogs_list.txt"
    awk -v score="$mgs_score" \
    '!/^#/ && $2 > score {
        if (match($3, /[0-9]{4}/))
            print "COG" substr($3, RSTART, RLENGTH)
    }' "$fetchmgs_dir"/*.scores 2>/dev/null \
    | sort -u > "$cog_list_file"

    if [ ! -s "$cog_list_file" ]; then
        echo "Error: No COGs found with bit-score > ${mgs_score}. Fast mode cannot proceed."
        echo "       Check FetchMGs output, lower --score, or use -m for precise mode."
        exit 1
    fi

    cog_count=$(wc -l < "$cog_list_file")
    echo "    Qualifying COG count: $cog_count"

    echo "    [Run] Extracting MGS genomic coordinates..."
    PY_BIN="python3"
    command -v python3 >/dev/null 2>&1 || PY_BIN="python"

    "$PY_BIN" - \
        --input "$output_host_faa" \
        --scores "$fetchmgs_dir" \
        --targets "$cog_list_file" \
        --outdir "$output/mgs_ext_results" \
    <<'PYEOF'
import argparse, sys, csv, re
from pathlib import Path

def eprint(*a, **k): print(*a, file=sys.stderr, **k)

def base_stem(p):
    n = p.name
    if n.endswith(".gff") or n.endswith(".faa"):
        return n[:-4]
    return p.stem

def discover_pairs(base):
    gffs, faas = {}, {}
    for q in base.rglob("*"):
        if not q.is_file(): continue
        if q.name.endswith(".gff"):   gffs[base_stem(q)] = q
        elif q.name.endswith(".faa"): faas[base_stem(q)] = q
    pairs = {k: (gffs[k], faas[k]) for k in gffs.keys() & faas.keys()}
    eprint(f"[INFO] Paired (gff, faa) genomes: {len(pairs)}")
    return pairs

def parse_gff_by_id(gff):
    by_id = {}
    with gff.open("rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS": continue
            seqid, _, _, start, end, _, strand, _, attr = parts
            m = re.search(r"ID=([^;]+)", attr)
            if not m: continue
            gid = m.group(1)
            lt   = re.search(r"locus_tag=([^;]+)", attr)
            gene = re.search(r"(?:gene|Name)=([^;]+)", attr)
            prod = re.search(r"product=([^;]+)", attr)
            by_id[gid] = {
                "contig": seqid, "start": start, "end": end, "strand": strand,
                "locus_tag": (lt.group(1)   if lt   else ""),
                "gene":      (gene.group(1) if gene else ""),
                "product":   (prod.group(1) if prod else ""),
            }
    return by_id

def build_pid2gid(faa):
    m = {}
    with faa.open("rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"): continue
            pid = line[1:].split()[0]
            mid = re.search(r"\bID=([^;]+)", line)
            if mid and pid not in m:
                m[pid] = mid.group(1)
    return m

def load_scores_for_genome(scores_dir, genome_name):
    candidates = [
        scores_dir / f"{genome_name}.faa.fetchMGs.scores",
        scores_dir / genome_name / f"{genome_name}.faa.fetchMGs.scores",
    ]
    sp = next((c for c in candidates if c.is_file()), None)
    if sp is None:
        return {}
    m = {}
    with sp.open("rt", encoding="utf-8", errors="replace") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#") or len(row) < 3: continue
            pid = row[0].strip()
            try:    score = float(row[1])
            except: continue
            mm = re.search(r"(\d{4})", row[2].strip())
            if not mm: continue
            cog = f"COG{mm.group(1)}"
            if pid not in m or score > m[pid][1]:
                m[pid] = (cog, score)
    return m

def read_targets(s):
    if not s: return None
    p = Path(s)
    vals = []
    if p.exists():
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if line and not line.startswith("#"): vals.append(line)
    out = set()
    for v in vals:
        mm = re.search(r"(\d{4})", v)
        if mm: out.add(f"COG{mm.group(1)}")
    return out if out else None

ap = argparse.ArgumentParser()
ap.add_argument("--input",   required=True)
ap.add_argument("--scores",  required=True)
ap.add_argument("--targets", required=True)
ap.add_argument("--outdir",  required=True)
args = ap.parse_args()

pairs = discover_pairs(Path(args.input))
if not pairs:
    eprint("ERROR: No (gff, faa) pairs found."); sys.exit(2)

scores_dir = Path(args.scores)
targets = read_targets(args.targets)

merged = []
for name in sorted(pairs.keys()):
    gff, faa = pairs[name]
    sub_scores = load_scores_for_genome(scores_dir, name)
    if not sub_scores:
        continue
    by_id   = parse_gff_by_id(gff)
    pid2gid = build_pid2gid(faa)
    for pid, (cog, score) in sub_scores.items():
        if targets and cog not in targets: continue
        gid = pid2gid.get(pid)
        if not gid: continue
        rec = by_id.get(gid)
        if not rec: continue
        merged.append({
            "genome": name,
            "contig": rec["contig"], "start": rec["start"],
            "end": rec["end"], "strand": rec["strand"],
            "locus_tag": rec.get("locus_tag") or gid,
            "gene": rec.get("gene", ""), "product": rec.get("product", ""),
            "COG": cog, "bit_score": score,
        })

outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
out = outdir / "all_from_scores.tsv"
with out.open("wt", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["genome","contig","start","end","strand",
                "locus_tag","gene","product","COG","bit_score"])
    for r in merged:
        w.writerow([r["genome"], r["contig"], r["start"], r["end"],
                    r["strand"], r["locus_tag"], r["gene"], r["product"],
                    r["COG"], r["bit_score"]])
eprint(f"[DONE] MGS coordinates written: {out} (rows={len(merged)})")
PYEOF

    if [ -s "$output/mgs_ext_results/all_from_scores.tsv" ]; then
        awk -F'\t' -v ms="$mgs_score" 'NR>1 && ($10+0) > ms {print $2 "\t" $3 "\t" $4}' \
            "$output/mgs_ext_results/all_from_scores.tsv" > "$output/mgs_ordinate.tsv"
    else
        echo "Error: Extracted MGS table is empty. Fast mode cannot proceed."
        exit 1
    fi

    if [ ! -s "$output/mgs_ordinate.tsv" ]; then
        echo "Error: Filtered MGS coordinate file is empty. Fast mode cannot proceed."
        exit 1
    fi
    echo "    MGS coordinate file ready: $output/mgs_ordinate.tsv"

    mgs_missing_log="$output/mgs_missing_genomes.log"
    : > "$mgs_missing_log"
    : > "$genomes_with_mgs_file"
    echo ">>> Validating MGS coverage across reference genomes (--mask ${mask_bp} bp)..."
    missing_count=0
    for hostgenome in "$in_hostgenome"/*.fasta; do
        [ -e "$hostgenome" ] || continue
        host_name=$(basename "$hostgenome" .fasta)
        if genome_has_mgs "$hostgenome" "$output/mgs_ordinate.tsv"; then
            echo "$host_name" >> "$genomes_with_mgs_file"
        else
            echo "    [Missing] $host_name has no valid MGS; skipping statistical testing in fast mode"
            echo "$host_name" >> "$mgs_missing_log"
            missing_count=$((missing_count + 1))
        fi
    done

    if [ -s "$genomes_with_mgs_file" ]; then
        echo "    Reference genomes with valid MGS background: $(wc -l < "$genomes_with_mgs_file")"
    else
        echo "Error: None of the reference genomes contain valid MGS coordinates."
        exit 1
    fi
    if [ "$missing_count" -gt 0 ]; then
        echo "    $missing_count genome(s) lacked valid MGS (logged to $mgs_missing_log)."
        [ "$map_mode" = "competitive" ] && echo "    In competitive mode, these genomes still participate in read alignment."
    fi
else
    echo ">>> [Precise Mode] Skipping Prodigal and FetchMGs preprocessing."
    : > "$genomes_with_mgs_file"
fi

echo ">>> [Step 3] Running read alignment and depth extraction (map-mode=$map_mode, mask=${mask_bp}bp)..."

output_map_dir="$output/mapping_results/$map_mode"
mkdir -p "$output_map_dir"

for data1 in "$in_seq"/*_1_clean.fastq "$in_seq"/*_1_clean.fastq.gz \
             "$in_seq"/*_1_clean.fq "$in_seq"/*_1_clean.fq.gz; do
    [ -e "$data1" ] || continue

    case "$data1" in
        *_1_clean.fastq.gz) base_name=$(basename "$data1" _1_clean.fastq.gz); data2="$in_seq/${base_name}_2_clean.fastq.gz" ;;
        *_1_clean.fastq)    base_name=$(basename "$data1" _1_clean.fastq);    data2="$in_seq/${base_name}_2_clean.fastq" ;;
        *_1_clean.fq.gz)    base_name=$(basename "$data1" _1_clean.fq.gz);    data2="$in_seq/${base_name}_2_clean.fq.gz" ;;
        *_1_clean.fq)       base_name=$(basename "$data1" _1_clean.fq);       data2="$in_seq/${base_name}_2_clean.fq" ;;
        *) continue ;;
    esac

    if [ ! -f "$data2" ]; then
        echo "    Warning: Missing paired R2 file $data2; skipping sample $base_name"
        continue
    fi

    find_longread_for_base "$base_name" "$in_seq"
    if [ -n "$LR_FILE" ]; then
        echo "    [Hybrid] Sample $base_name paired with long reads (${LR_TYPE}, ${LR_PRESET}): $(basename "$LR_FILE")"
    else
        echo "    [Short-Read PE] Sample $base_name"
    fi

    ref_str=$(resolve_ref_str "$base_name")
    ref_str=$(normalize_ref_str "$ref_str")
    if [ -z "$ref_str" ]; then
        echo "    [Skip] Sample $base_name not found in index table"
        continue
    fi

    if [ "$map_mode" = "competitive" ]; then
        process_competitive_sample "$base_name" "$ref_str" "$data1" "$data2" "$LR_FILE" "$LR_PRESET"
    else
        IFS=',' read -ra REF_GENOMES <<< "$ref_str"
        for ref_name in "${REF_GENOMES[@]}"; do
            ref_name="${ref_name// /}"
            [ -z "$ref_name" ] && continue
            process_one_sample_ref "$base_name" "$ref_name" "$data1" "$data2" "$LR_FILE" "$LR_PRESET"
        done
    fi
done

shopt -s nullglob
lr_candidates=(
    "$in_seq"/*_pbhifi.fastq.gz "$in_seq"/*_pbhifi.fastq
    "$in_seq"/*_pacbio.fastq.gz "$in_seq"/*_pacbio.fastq
    "$in_seq"/*_ont.fastq.gz    "$in_seq"/*_ont.fastq
    "$in_seq"/*_pbhifi.fq.gz    "$in_seq"/*_pbhifi.fq
    "$in_seq"/*_pacbio.fq.gz    "$in_seq"/*_pacbio.fq
    "$in_seq"/*_ont.fq.gz       "$in_seq"/*_ont.fq
)
shopt -u nullglob

declare -A seen_lr_bases

for lr_path in "${lr_candidates[@]}"; do
    [ -e "$lr_path" ] || continue

    lr_name=$(basename "$lr_path")
    lr_stem="${lr_name%.gz}"
    lr_stem="${lr_stem%.fastq}"
    lr_stem="${lr_stem%.fq}"

    case "$lr_stem" in
        *_pbhifi) base_name="${lr_stem%_pbhifi}" ;;
        *_pacbio) base_name="${lr_stem%_pacbio}" ;;
        *_ont)    base_name="${lr_stem%_ont}" ;;
        *) continue ;;
    esac

    [ -n "${seen_lr_bases[$base_name]:-}" ] && continue
    seen_lr_bases["$base_name"]=1

    has_illumina_r1 "$base_name" "$in_seq" && continue

    find_longread_for_base "$base_name" "$in_seq"
    [ -z "$LR_FILE" ] && continue

    echo "    [Long-Read Only] Sample $base_name (${LR_TYPE}, ${LR_PRESET}): $(basename "$LR_FILE")"

    ref_str=$(resolve_ref_str "$base_name")
    ref_str=$(normalize_ref_str "$ref_str")
    if [ -z "$ref_str" ]; then
        echo "    [Skip] Sample $base_name not found in index table"
        continue
    fi

    if [ "$map_mode" = "competitive" ]; then
        process_competitive_sample "$base_name" "$ref_str" "" "" "$LR_FILE" "$LR_PRESET"
    else
        IFS=',' read -ra REF_GENOMES <<< "$ref_str"
        for ref_name in "${REF_GENOMES[@]}"; do
            ref_name="${ref_name// /}"
            [ -z "$ref_name" ] && continue
            process_one_sample_ref "$base_name" "$ref_name" "" "" "$LR_FILE" "$LR_PRESET"
        done
    fi
done

echo ">>> Depth extraction stage completed."

final_res_dir="$output/activity_results/$mode/$map_mode"
mkdir -p "$final_res_dir"
summary_tsv="$final_res_dir/${mode}_activity_summary.tsv"
basic_tsv="$final_res_dir/${mode}_basic_stats.tsv"

echo ">>> [Step 4] Running R statistical analysis (mode=$mode, map-mode=$map_mode)..."

if [ "$mode" = "fast" ]; then

Rscript - \
    --threads "$r_threads" \
    --mapdir "$output_map_dir" \
    --coverage "$coverage_threshold" \
    --cutoff "$cutoff" \
    --conservative-off "$conservative_off" \
    --out-summary "$summary_tsv" \
    --out-basic "$basic_tsv" \
<<'REOF_FAST'
suppressPackageStartupMessages({
  library(tseries); library(brunnermunzel)
  library(effsize); library(doParallel); library(foreach)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", flag, "="), "", hit[1]))
  idx <- which(args == flag)
  if (length(idx) > 0 && length(args) >= idx[1] + 1L) return(args[idx[1] + 1L])
  default
}

cores <- suppressWarnings(as.integer(get_arg("--threads", "20")))
if (is.na(cores) || cores < 1L) cores <- 20L

coverage_threshold <- suppressWarnings(as.numeric(get_arg("--coverage", "0.75")))
if (is.na(coverage_threshold)) coverage_threshold <- 0.75

cutoff <- suppressWarnings(as.numeric(get_arg("--cutoff", "0.8")))
if (is.na(cutoff)) cutoff <- 0.8

conservative_off <- tolower(get_arg("--conservative-off", "false")) %in% c("true", "1", "yes")
root_dir <- get_arg("--mapdir", "")
out_summary <- get_arg("--out-summary", "")
out_basic <- get_arg("--out-basic", "")

if (!nzchar(root_dir) || !dir.exists(root_dir)) stop("Invalid --mapdir", call. = FALSE)
if (!nzchar(out_summary) || !nzchar(out_basic)) stop("Missing output TSV paths", call. = FALSE)

cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
message(sprintf("[R-Fast] Cores: %d | Coverage >= %.4f | Cutoff >= %.4f | Conservative_Off: %s",
                cores, coverage_threshold, cutoff, conservative_off))

sample_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[!grepl("global$", basename(sample_dirs))]
alpha <- 0.01

all_results <- data.frame(
  Folder=character(), Genome=character(), Compare=character(),
  Phage_File=character(), Coverage=numeric(), Test_Type=character(),
  p_value=numeric(), Significant=character(), Effect_Size_Type=character(),
  Effect_Size_Value=numeric(), BM_Prob_Superiority=numeric(), Phage_Sites=numeric(),
  prediction=character(), stringsAsFactors=FALSE)

all_stat_results <- data.frame(
  Folder=character(), Genome=character(), Region=character(),
  File=character(), Mean=numeric(), Median=numeric(), Variance=numeric(), COVERAGE=numeric(),
  stringsAsFactors=FALSE)

for (sample_dir in sample_dirs) {
  folder_name <- basename(sample_dir)
  message(sprintf("======== Sample: %s ========", folder_name))
  stat_main_dir <- file.path(sample_dir, "depth_to_stat")
  if (!dir.exists(stat_main_dir)) { message("  No depth_to_stat; skipping."); next }
  target_dirs <- list.dirs(stat_main_dir, full.names = TRUE, recursive = FALSE)

  for (datadir in target_dirs) {
    genomename <- basename(datadir)
    mgs_file        <- list.files(datadir, pattern="_mgs_depth\\.txt$",     full.names=TRUE)
    phage_files_raw <- list.files(datadir, pattern="^phage_depth.*\\.txt$", full.names=TRUE)
    if (length(mgs_file)==0)        { message(sprintf("  [%s] Missing MGS depth; skipping.", genomename)); next }
    gi <- file.info(mgs_file[1])
    if (is.na(gi$size)||gi$size==0) { message(sprintf("  [%s] Empty MGS depth; skipping.", genomename)); next }
    if (length(phage_files_raw)==0) { message(sprintf("  [%s] Missing phage depth; skipping.", genomename)); next }
    pfi <- file.info(phage_files_raw)
    phage_files <- phage_files_raw[pfi$size > 0]
    if (length(phage_files)==0)     { message(sprintf("  [%s] All phage depth files empty; skipping.", genomename)); next }

    mgs_df <- tryCatch(read.table(mgs_file[1], header=FALSE), error=function(e) NULL)
    if (is.null(mgs_df)||ncol(mgs_df)<3) next
    colnames(mgs_df) <- c("id","pos","depth")
    mgs_depth <- as.numeric(mgs_df$depth)
    mgs_depth <- mgs_depth[is.finite(mgs_depth)]
    if (length(mgs_depth) < 2) next

    all_stat_results <- rbind(all_stat_results, data.frame(
      Folder=folder_name, Genome=genomename, Region="host_mgs",
      File=basename(mgs_file[1]), Mean=mean(mgs_depth),
      Median=median(mgs_depth), Variance=var(mgs_depth), COVERAGE=NA_real_,
      stringsAsFactors=FALSE))

    plist <- foreach(f=phage_files,
      .packages=c("tseries","brunnermunzel","effsize"), .errorhandling="pass") %dopar% {
      phage <- tryCatch(read.table(f, header=FALSE), error=function(e) NULL)
      if (is.null(phage)||ncol(phage)<3) return(NULL)
      colnames(phage) <- c("id","pos","depth")
      phage_depth <- as.numeric(phage$depth)
      phage_depth <- phage_depth[is.finite(phage_depth)]
      phage_sites <- length(phage_depth)
      if (phage_sites == 0) return(NULL)

      coverage <- sum(phage_depth > 0) / phage_sites
      stat_line <- data.frame(
        Folder=folder_name, Genome=genomename, Region="phage", File=basename(f),
        Mean=mean(phage_depth), Median=median(phage_depth),
        Variance=if (phage_sites>1) var(phage_depth) else NA_real_, COVERAGE=round(coverage,4),
        stringsAsFactors=FALSE)

      if (coverage < coverage_threshold || phage_sites < 2) {
        return(list(stat_line=stat_line, result_line=data.frame(
          Folder=folder_name, Genome=genomename, Compare="phage vs mgs",
          Phage_File=basename(f), Coverage=round(coverage,4),
          Test_Type=NA_character_, p_value=NA_real_, Significant="Coverage failed",
          Effect_Size_Type=NA_character_, Effect_Size_Value=NA_real_,
          BM_Prob_Superiority=NA_real_, Phage_Sites=phage_sites,
          prediction="low_coverage", stringsAsFactors=FALSE)))
      }

      norm1 <- FALSE
      norm2 <- FALSE
      if (conservative_off) {
        jb1   <- tryCatch(tseries::jarque.bera.test(phage_depth), error=function(e) NULL)
        jb2   <- tryCatch(tseries::jarque.bera.test(mgs_depth),   error=function(e) NULL)
        norm1 <- !is.null(jb1) && !is.na(jb1$p.value) && jb1$p.value > 0.05
        norm2 <- !is.null(jb2) && !is.na(jb2$p.value) && jb2$p.value > 0.05
      }

      test_type <- NA_character_; pval <- NA_real_; significant <- "No"
      effect_type <- NA_character_; effect_value <- NA_real_; bm_prob <- NA_real_
      pred <- "dormant"

      if (norm1 && norm2) {
        vtest <- tryCatch(var.test(phage_depth, mgs_depth), error=function(e) NULL)
        eq_var <- !is.null(vtest) && !is.na(vtest$p.value) && vtest$p.value > 0.05
        ttest <- t.test(phage_depth, mgs_depth, alternative="greater", var.equal=eq_var)
        test_type <- if (eq_var) "Student's t-test" else "Welch's t-test"
        pval <- ttest$p.value
        if (!is.na(pval) && pval < alpha) {
          significant <- "Yes"
          cd <- tryCatch(effsize::cohen.d(phage_depth, mgs_depth, hedges.correction=TRUE), error=function(e) NULL)
          if (!is.null(cd)) {
            effect_type <- "Cohen's D"
            effect_value <- as.numeric(cd$estimate)
            if (is.finite(effect_value) && effect_value >= 0.7) pred <- "active"
          }
        }
      } else {
        bmtest <- tryCatch(brunnermunzel::brunnermunzel.test(phage_depth, mgs_depth, alternative="greater"), error=function(e) NULL)
        test_type <- "Brunner-Munzel"
        if (!is.null(bmtest)) {
          pval    <- bmtest$p.value
          est_raw <- suppressWarnings(as.numeric(bmtest$estimate))
          if (is.finite(est_raw)) bm_prob <- 1 - est_raw
        }
        if (!is.na(pval) && pval < alpha) {
          significant <- "Yes"
          cd <- tryCatch(suppressWarnings(effsize::cliff.delta(phage_depth, mgs_depth)), error=function(e) NULL)
          if (!is.null(cd)) {
            effect_type <- "Cliff's Delta"
            effect_value <- as.numeric(cd$estimate)
            if (is.finite(effect_value) && effect_value >= cutoff) pred <- "active"
          }
        }
      }

      list(stat_line=stat_line, result_line=data.frame(
        Folder=folder_name, Genome=genomename, Compare="phage vs mgs",
        Phage_File=basename(f), Coverage=round(coverage,4),
        Test_Type=test_type, p_value=pval, Significant=significant,
        Effect_Size_Type=effect_type, Effect_Size_Value=effect_value,
        BM_Prob_Superiority=bm_prob, Phage_Sites=phage_sites,
        prediction=pred, stringsAsFactors=FALSE))
    }

    plist <- Filter(function(x) is.list(x) && all(c("stat_line","result_line") %in% names(x)), plist)
    if (length(plist) > 0) {
      all_stat_results <- rbind(all_stat_results, do.call(rbind, lapply(plist, `[[`, "stat_line")))
      all_results      <- rbind(all_results,      do.call(rbind, lapply(plist, `[[`, "result_line")))
    }
  }
}

write.table(all_results,      out_summary, sep="\t", row.names=FALSE, quote=FALSE)
write.table(all_stat_results, out_basic,   sep="\t", row.names=FALSE, quote=FALSE)
REOF_FAST

else

Rscript - \
    --threads "$r_threads" \
    --mapdir "$output_map_dir" \
    --coverage "$coverage_threshold" \
    --cutoff "$cutoff" \
    --conservative-off "$conservative_off" \
    --iter "$n_iter" \
    --out-summary "$summary_tsv" \
    --out-basic "$basic_tsv" \
<<'REOF_PRECISE'
suppressPackageStartupMessages({
  library(tseries); library(brunnermunzel)
  library(effsize); library(doParallel); library(foreach)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", flag, "="), "", hit[1]))
  idx <- which(args == flag)
  if (length(idx) > 0 && length(args) >= idx[1] + 1L) return(args[idx[1] + 1L])
  default
}

cores <- suppressWarnings(as.integer(get_arg("--threads", "20")))
if (is.na(cores) || cores < 1L) cores <- 20L

coverage_threshold <- suppressWarnings(as.numeric(get_arg("--coverage", "0.75")))
if (is.na(coverage_threshold)) coverage_threshold <- 0.75

cutoff <- suppressWarnings(as.numeric(get_arg("--cutoff", "0.8")))
if (is.na(cutoff)) cutoff <- 0.8

conservative_off <- tolower(get_arg("--conservative-off", "false")) %in% c("true", "1", "yes")

n_iter <- suppressWarnings(as.integer(get_arg("--iter", "500")))
if (is.na(n_iter) || n_iter < 1L) n_iter <- 500L

root_dir <- get_arg("--mapdir", "")
out_summary <- get_arg("--out-summary", "")
out_basic <- get_arg("--out-basic", "")

if (!nzchar(root_dir) || !dir.exists(root_dir)) stop("Invalid --mapdir", call. = FALSE)
if (!nzchar(out_summary) || !nzchar(out_basic)) stop("Missing output TSV paths", call. = FALSE)

alpha <- 0.01
cl <- parallel::makeCluster(cores)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
doParallel::registerDoParallel(cl)
message(sprintf("[R-Precise] Cores: %d | Coverage >= %.4f | Cutoff >= %.4f | Conservative_Off: %s | Iter: %d",
                cores, coverage_threshold, cutoff, conservative_off, n_iter))

sample_dirs <- list.dirs(root_dir, full.names=TRUE, recursive=FALSE)
sample_dirs <- sample_dirs[!grepl("global$", basename(sample_dirs))]

all_results <- data.frame(
  Folder=character(), Genome=character(), Compare=character(),
  Phage_File=character(), Coverage=numeric(), Test_Type=character(),
  p_value=numeric(), Significant=character(), BM_Prob_Superiority=numeric(),
  RR_Q95=numeric(), RR_Q99=numeric(),
  Delta_Mean=numeric(), Delta_CI_2.5=numeric(), Delta_CI_97.5=numeric(),
  prediction=character(), stringsAsFactors=FALSE)

all_stat_results <- data.frame(
  Folder=character(), Genome=character(), Region=character(),
  File=character(), Mean=numeric(), Median=numeric(), Variance=numeric(), COVERAGE=numeric(),
  stringsAsFactors=FALSE)

for (s_idx in seq_along(sample_dirs)) {
  sample_dir  <- sample_dirs[s_idx]
  folder_name <- basename(sample_dir)
  message(sprintf("======== [%d/%d] Sample: %s ========", s_idx, length(sample_dirs), folder_name))
  stat_main_dir <- file.path(sample_dir, "depth_to_stat")
  if (!dir.exists(stat_main_dir)) next
  target_dirs <- list.dirs(stat_main_dir, full.names=TRUE, recursive=FALSE)

  for (datadir in target_dirs) {
    genomename  <- basename(datadir)
    host_file   <- list.files(datadir, pattern="_host_nonphage_depth\\.txt$", full.names=TRUE)
    phage_files <- list.files(datadir, pattern="^phage_depth.*\\.txt$",        full.names=TRUE)
    if (length(host_file)==0 || length(phage_files)==0) next

    hi <- file.info(host_file[1])
    if (is.na(hi$size) || hi$size == 0) next
    host <- tryCatch(read.table(host_file[1], header=FALSE), error=function(e) NULL)
    if (is.null(host) || ncol(host) < 3) next
    colnames(host) <- c("id","pos","depth")
    host_depth <- as.numeric(host$depth); host_depth <- host_depth[is.finite(host_depth)]
    if (length(host_depth) < 2) next
    nonphage_sites <- length(host_depth)

    all_stat_results <- rbind(all_stat_results, data.frame(
      Folder=folder_name, Genome=genomename, Region="host_non_phage",
      File=basename(host_file[1]), Mean=mean(host_depth),
      Median=median(host_depth), Variance=stats::var(host_depth), COVERAGE=NA_real_,
      stringsAsFactors=FALSE))

    phage_all <- unlist(lapply(phage_files, function(f) {
      df <- tryCatch(read.table(f, header=FALSE), error=function(e) NULL)
      if (is.null(df) || ncol(df) < 3) return(numeric(0)); as.numeric(df[[3]])
    }))
    phage_all <- phage_all[is.finite(phage_all)]
    all_d <- c(phage_all, host_depth)
    if (length(all_d) == 0) next

    q95 <- as.numeric(quantile(all_d, 0.95, names=FALSE, type=7))
    q99 <- as.numeric(quantile(all_d, 0.99, names=FALSE, type=7))
    p0_95 <- sum(host_depth > q95) / nonphage_sites
    p0_99 <- sum(host_depth > q99) / nonphage_sites

    for (f in phage_files) {
      phage <- tryCatch(read.table(f, header=FALSE), error=function(e) NULL)
      if (is.null(phage) || ncol(phage) < 3) next
      colnames(phage) <- c("id","pos","depth")
      phage_depth <- as.numeric(phage$depth); phage_depth <- phage_depth[is.finite(phage_depth)]
      phage_sites <- length(phage_depth)
      if (phage_sites == 0) next

      coverage <- sum(phage_depth > 0) / phage_sites
      all_stat_results <- rbind(all_stat_results, data.frame(
        Folder=folder_name, Genome=genomename, Region="phage", File=basename(f),
        Mean=mean(phage_depth), Median=median(phage_depth),
        Variance=if (phage_sites>1) stats::var(phage_depth) else NA_real_, COVERAGE=round(coverage,4),
        stringsAsFactors=FALSE))

      if (coverage < coverage_threshold || phage_sites < 2) {
        all_results <- rbind(all_results, data.frame(
          Folder=folder_name, Genome=genomename, Compare="phage vs host_nonphage",
          Phage_File=basename(f), Coverage=round(coverage,4),
          Test_Type=NA_character_, p_value=NA_real_, Significant="Coverage failed",
          BM_Prob_Superiority=NA_real_, RR_Q95=NA_real_, RR_Q99=NA_real_,
          Delta_Mean=NA_real_, Delta_CI_2.5=NA_real_, Delta_CI_97.5=NA_real_,
          prediction="low_coverage", stringsAsFactors=FALSE))
        next
      }

      test_type <- NA_character_; pval <- NA_real_; significant <- "No"
      bm_prob <- NA_real_; rr_95 <- NA_real_; rr_99 <- NA_real_
      delta_mean <- NA_real_; delta_CI <- c(NA_real_, NA_real_)
      pred <- "dormant"

      norm1 <- FALSE
      norm2 <- FALSE
      if (conservative_off) {
        jb1   <- tryCatch(tseries::jarque.bera.test(phage_depth), error=function(e) NULL)
        jb2   <- tryCatch(tseries::jarque.bera.test(host_depth),  error=function(e) NULL)
        norm1 <- !is.null(jb1) && !is.na(jb1$p.value) && jb1$p.value > 0.05
        norm2 <- !is.null(jb2) && !is.na(jb2$p.value) && jb2$p.value > 0.05
      }

      if (norm1 && norm2) {
        vtest <- tryCatch(var.test(phage_depth, host_depth), error=function(e) NULL)
        eq_var <- !is.null(vtest) && !is.na(vtest$p.value) && vtest$p.value > 0.05
        ttest <- t.test(phage_depth, host_depth, alternative="greater", var.equal=eq_var)
        test_type <- if (eq_var) "Student's t-test" else "Welch's t-test"
        pval <- ttest$p.value
        if (!is.na(pval) && pval < alpha) {
          significant <- "Yes"
          cd <- tryCatch(effsize::cohen.d(phage_depth, host_depth, hedges.correction=TRUE), error=function(e) NULL)
          if (!is.null(cd)) {
            cohen_val <- suppressWarnings(as.numeric(cd$estimate))
            if (is.finite(cohen_val) && cohen_val >= 0.7) pred <- "active"
          }
        }
      } else {
        bmtest <- tryCatch(brunnermunzel::brunnermunzel.test(phage_depth, host_depth, alternative="greater"), error=function(e) NULL)
        test_type <- "Brunner-Munzel"
        if (!is.null(bmtest)) {
          pval    <- bmtest$p.value
          est_raw <- suppressWarnings(as.numeric(bmtest$estimate))
          if (is.finite(est_raw)) bm_prob <- 1 - est_raw
        }
        if (!is.na(pval) && pval < alpha) {
          significant    <- "Yes"
          phage_over_q95 <- sum(phage_depth > q95)
          phage_over_q99 <- sum(phage_depth > q99)
          rr_95 <- ifelse(p0_95==0, NA_real_, (phage_over_q95 / phage_sites) / p0_95)
          rr_99 <- ifelse(p0_99==0, NA_real_, (phage_over_q99 / phage_sites) / p0_99)

          set.seed(123)
          len_h <- min(10000L, length(host_depth))
          cliff_deltas <- foreach(i=1:n_iter, .combine=c, .packages="effsize") %dopar% {
            sub_h <- host_depth[sample.int(length(host_depth), len_h, replace=FALSE)]
            res   <- try(effsize::cliff.delta(phage_depth, sub_h), silent=TRUE)
            if (inherits(res, "try-error")) NA_real_ else suppressWarnings(as.numeric(res$estimate))
          }
          cliff_deltas <- cliff_deltas[is.finite(cliff_deltas)]
          if (length(cliff_deltas) > 0) {
            delta_mean <- mean(cliff_deltas, na.rm=TRUE)
            qs <- stats::quantile(cliff_deltas, c(0.025, 0.975), na.rm=TRUE, names=FALSE, type=7)
            delta_CI <- c(qs[1], qs[2])
            if (is.finite(delta_mean) && delta_mean >= cutoff) pred <- "active"
          }
        }
      }

      all_results <- rbind(all_results, data.frame(
        Folder=folder_name, Genome=genomename, Compare="phage vs host_nonphage",
        Phage_File=basename(f), Coverage=round(coverage,4),
        Test_Type=test_type, p_value=pval, Significant=significant,
        BM_Prob_Superiority=bm_prob, RR_Q95=rr_95, RR_Q99=rr_99,
        Delta_Mean=delta_mean, Delta_CI_2.5=delta_CI[1], Delta_CI_97.5=delta_CI[2],
        prediction=pred, stringsAsFactors=FALSE))
    }
  }
}

write.table(all_results,      out_summary, sep="\t", row.names=FALSE, quote=FALSE)
write.table(all_stat_results, out_basic,   sep="\t", row.names=FALSE, quote=FALSE)
REOF_PRECISE

fi

echo ">>> Statistical analysis completed:"
echo "    Summary table : $summary_tsv"
echo "    Basic stats   : $basic_tsv"

if [ "$do_plot" = true ]; then
    plot_out_dir="$final_res_dir/plots"
    echo ">>> [Step 4.5] Running visualization module (mode=$mode, plot-all=$plot_all, overview=$plot_overview)..."

    Rscript - \
        --map-dir "$output_map_dir" \
        --summary "$summary_tsv" \
        --out-dir "$plot_out_dir" \
        --genome-dir "$in_hostgenome" \
        --mode "$mode" \
        --plot-all "$plot_all" \
        --overview "$plot_overview" \
        --buffer "$plot_buffer" \
        --min-contig "$plot_min_contig" \
    <<'REOF_PLOT'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

FONT_FAMILY   <- "Arial"
BIN_SIZE      <- 100L
CONTIG_EXPAND <- 0.9

label_host    <- "Host Background"
label_mgs     <- "Host MGS"
label_active  <- "Active Phage"
label_dormant <- "Dormant Phage"
label_lowcov  <- "Low-coverage Phage"

color_contig_hl <- "#EFEFEF"
color_phage_hl  <- "#D0E1F9"

group_levels <- c(label_host, label_mgs, label_active, label_dormant, label_lowcov)
custom_colors <- setNames(
  c("#222222", "#377EB8", "#E41A1C", "#E69F00", "#888888"),
  group_levels
)
line_widths <- setNames(
  c(1.5, 0.55, 0.75, 0.65, 0.60),
  group_levels
)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", flag, "="), "", hit[1]))
  idx <- which(args == flag)
  if (length(idx) > 0 && length(args) >= idx[1] + 1L) return(args[idx[1] + 1L])
  default
}

map_dir    <- get_arg("--map-dir", "")
sum_file   <- get_arg("--summary", "")
out_dir    <- get_arg("--out-dir", "")
genome_dir <- get_arg("--genome-dir", "")
run_mode   <- get_arg("--mode", "fast")
plot_all   <- tolower(get_arg("--plot-all", "false")) %in% c("true", "1", "yes")
overview   <- tolower(get_arg("--overview", "false")) %in% c("true", "1", "yes") && identical(run_mode, "precise")
buffer_bp  <- suppressWarnings(as.integer(get_arg("--buffer", "30000")))
min_contig <- suppressWarnings(as.integer(get_arg("--min-contig", "8000")))

if (is.na(buffer_bp) || buffer_bp < 0L) buffer_bp <- 30000L
if (is.na(min_contig) || min_contig < 0L) min_contig <- 8000L
if (!nzchar(genome_dir)) genome_dir <- NULL

if (!file.exists(sum_file)) {
  message("[Plot] Summary table not found; skipping visualization.")
  quit(save = "no", status = 0)
}

sum_dt <- tryCatch(fread(sum_file, sep = "\t", header = TRUE), error = function(e) NULL)
if (is.null(sum_dt) || !nrow(sum_dt)) {
  message("[Plot] Summary table has 0 prophage rows; skipping visualization.")
  quit(save = "no", status = 0)
}

sum_dt[, Folder := as.character(Folder)]
sum_dt[, Genome := as.character(Genome)]
sum_dt[, Phage_File := as.character(Phage_File)]
sum_dt[, prediction := as.character(prediction)]

if (plot_all) {
  target_dt <- copy(sum_dt)
} else {
  target_dt <- sum_dt[prediction == "active"]
}

if (!nrow(target_dt)) {
  message("[Plot] No prophages matched plotting filter (prediction == 'active'); skipping visualization.")
  quit(save = "no", status = 0)
}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pred_to_label <- function(pred) {
  if (is.na(pred)) return(label_dormant)
  if (pred == "active") return(label_active)
  if (pred == "low_coverage") return(label_lowcov)
  label_dormant
}

fmt_int <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_stat <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !is.finite(as.numeric(x))) return("NA")
  v <- as.numeric(x)
  if (v != 0 && abs(v) < 1e-3) sprintf("%.2e", v) else sprintf("%.4f", v)
}
safe_token <- function(x) gsub("[^A-Za-z0-9._-]", "_", as.character(x))

find_existing <- function(...) {
  for (p in c(...)) if (!is.null(p) && file.exists(p)) return(p)
  NULL
}

read_depth <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  dt <- tryCatch(fread(path, header = FALSE, select = 1:3, col.names = c("contig", "position", "depth")),
                 error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  dt[, `:=`(position = as.numeric(position), depth = as.numeric(depth))]
  dt[is.finite(position) & is.finite(depth)][]
}

parse_phage_filename <- function(path) {
  bn <- sub("\\.gz$", "", basename(path))
  m <- str_match(bn, "^phage_depth_(.+)_([0-9]+)_([0-9]+)\\.txt$")
  if (anyNA(m[1, ])) return(NULL)
  list(contig = m[1, 2], start = as.integer(m[1, 3]), end = as.integer(m[1, 4]), bn = basename(path))
}

read_fai_lengths <- function(gdir, genome) {
  if (is.null(gdir)) return(NULL)
  fai <- find_existing(
    file.path(gdir, paste0(genome, ".fasta.fai")),
    file.path(gdir, paste0(genome, ".fa.fai")),
    file.path(gdir, paste0(genome, ".fna.fai")),
    file.path(gdir, paste0(genome, ".fai"))
  )
  if (is.null(fai)) return(NULL)
  dt <- fread(fai, header = FALSE, select = 1:2, col.names = c("contig", "length"))
  dt[, length := as.numeric(length)][]
}

read_bed_lengths <- function(data_dir, genome) {
  bed <- find_existing(
    file.path(data_dir, paste0(genome, ".valid.bed")),
    file.path(data_dir, paste0(genome, ".all.bed"))
  )
  if (is.null(bed)) return(NULL)
  dt <- fread(bed, header = FALSE, select = 1:3, col.names = c("contig", "start0", "end"))
  dt[, .(length = max(as.numeric(end))), by = contig]
}

infer_lengths <- function(dt_list) {
  keep <- dt_list[!vapply(dt_list, is.null, logical(1))]
  if (!length(keep)) return(NULL)
  all_dt <- rbindlist(keep, use.names = TRUE, fill = TRUE)
  all_dt[, .(length = max(position)), by = contig]
}

build_offsets <- function(len_dt) {
  len_dt <- copy(len_dt)[is.finite(length) & length > 0]
  if (!nrow(len_dt)) stop("No valid contig lengths")
  len_dt[, offset := c(0, cumsum(length)[-.N])]
  len_dt[, end_global := offset + length]
  len_dt[]
}

add_global_pos <- function(dt, off_dt) {
  if (is.null(dt) || !nrow(dt)) return(dt)
  out <- copy(dt)
  out[off_dt, on = "contig", global_pos := position + i.offset]
  out <- out[!is.na(global_pos)]
  if (nrow(out)) setorder(out, global_pos)
  out[]
}

filter_to_contigs <- function(dt, keep_contigs) {
  if (is.null(dt) || !nrow(dt)) return(dt)
  dt[contig %in% keep_contigs]
}

typical_step <- function(pos) {
  d <- diff(pos); d <- d[is.finite(d) & d > 0]
  if (!length(d)) return(1)
  as.numeric(stats::median(d))
}

mark_line_segments <- function(dt, gap = NULL) {
  if (is.null(dt) || !nrow(dt)) return(dt)
  out <- copy(dt)
  setorder(out, global_pos)
  step <- typical_step(out$global_pos)
  if (is.null(gap)) gap <- max(step * 8, 20)
  out[, seg := cumsum(c(0L, as.integer(diff(global_pos) > gap)))]
  out[]
}

downsample_depth <- function(dt, max_points = 80000L) {
  if (is.null(dt) || nrow(dt) <= max_points) return(copy(dt))
  out <- mark_line_segments(copy(dt))
  k <- max(1L, as.integer(ceiling(nrow(out) / max_points)))
  out[, `:=`(i = seq_len(.N), n = .N), by = seg]
  out <- out[i == 1L | i == n | (i %% k == 0L)]
  out[, `:=`(i = NULL, n = NULL, seg = NULL)][]
}

bin_mean_depth <- function(dt, bin_size) {
  if (is.null(dt) || !nrow(dt) || is.na(bin_size) || bin_size <= 1L) return(copy(dt))
  out <- copy(dt)
  out[, bin_id := as.integer(floor((pmax(position, 1) - 1) / bin_size))]
  if ("Group" %in% names(out)) {
    out <- out[, .(position = bin_id[1] * bin_size + (bin_size + 1) / 2,
                   depth = mean(depth, na.rm = TRUE)), by = .(contig, Group, bin_id)]
  } else {
    out <- out[, .(position = bin_id[1] * bin_size + (bin_size + 1) / 2,
                   depth = mean(depth, na.rm = TRUE)), by = .(contig, bin_id)]
  }
  out[, bin_id := NULL]
  setorder(out, contig, position)
  out[]
}

prepare_full_track <- function(dt, off_dt, bin_size) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  cols <- intersect(c("contig", "position", "depth", "Group"), names(dt))
  binned <- bin_mean_depth(dt[, ..cols], bin_size)
  add_global_pos(binned, off_dt)
}

safe_var_fix <- function(x) {
  if (length(x) >= 2 && (is.na(var(x)) || var(x) == 0)) {
    x[1] <- x[1] + 1e-6; x[2] <- x[2] - 1e-6
  }
  x
}

build_contig_groups <- function(phage_list, len_dt, expand_frac = CONTIG_EXPAND) {
  if (!length(phage_list) || is.null(len_dt) || !nrow(len_dt)) return(list())
  n <- nrow(len_dt)
  names_all <- len_dt$contig
  groups <- list()

  for (ph in phage_list) {
    idx <- match(ph$contig, names_all)
    if (is.na(idx)) next
    clen <- as.numeric(len_dt$length[idx])
    plen <- as.numeric(ph$end) - as.numeric(ph$start) + 1
    frac <- if (is.finite(clen) && clen > 0) plen / clen else NA_real_
    keep <- idx
    want_expand <- is.finite(frac) && frac >= expand_frac
    if (isTRUE(want_expand)) {
      if (idx > 1L) keep <- c(idx - 1L, keep)
      if (idx < n) keep <- c(keep, idx + 1L)
    }
    keep <- sort(unique(keep))
    contigs <- names_all[keep]
    key <- paste(contigs, collapse = "\001")

    if (is.null(groups[[key]])) {
      groups[[key]] <- list(contigs = contigs, center = ph$contig,
                            expanded = length(contigs) > 1L, phages = list())
    }
    groups[[key]]$phages[[length(groups[[key]]$phages) + 1L]] <- ph
  }
  unname(groups)
}

shared_theme <- function(font, legend = TRUE) {
  theme_classic(base_family = font) +
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 13, color = "black"),
      plot.subtitle     = element_text(hjust = 0.5, size = 9.5, color = "grey20"),
      axis.title        = element_text(size = 13, face = "bold", color = "black"),
      axis.text         = element_text(size = 11, color = "black"),
      axis.line         = element_line(color = "black", linewidth = 0.8),
      axis.ticks        = element_line(color = "black", linewidth = 0.8),
      legend.position   = if (legend) "top" else "none",
      legend.title      = element_blank(),
      legend.text       = element_text(size = 11),
      legend.background = element_blank()
    )
}

prep_line_dt <- function(dt, default_group = NULL, downsample = FALSE, gap = NULL) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  tmp <- copy(dt)
  if (!("Group" %in% names(tmp))) tmp[, Group := default_group]
  if (isTRUE(downsample)) tmp <- downsample_depth(tmp)
  out_list <- lapply(split(tmp, by = "Group"), function(sub_dt) {
    s <- mark_line_segments(sub_dt, gap = gap)
    s[, .(global_pos, depth, Group, seg)]
  })
  rbindlist(out_list)
}

make_profile <- function(dt_host, dt_phage, xmin, xmax, y_lab, title, subtitle,
                         contig_regions, phage_regions, hline_y,
                         font, show_labels = FALSE, legend = TRUE,
                         downsample_host = FALSE, line_gap = NULL) {
  host_dt  <- prep_line_dt(dt_host, label_host, downsample_host, gap = line_gap)
  phage_dt <- prep_line_dt(dt_phage, label_active, FALSE, gap = line_gap)

  p <- ggplot()

  if (!is.null(contig_regions) && nrow(contig_regions)) {
    c_sub <- contig_regions[end_g >= xmin & start_g <= xmax]
    if (nrow(c_sub)) {
      p <- p + geom_rect(
        data = c_sub,
        aes(xmin = pmax(start_g, xmin), xmax = pmin(end_g, xmax), ymin = -Inf, ymax = Inf),
        fill = color_contig_hl, color = "grey80", linetype = "dashed", linewidth = 0.4,
        inherit.aes = FALSE
      )
    }
  }

  if (!is.null(phage_regions) && nrow(phage_regions)) {
    p_sub <- phage_regions[end_g >= xmin & start_g <= xmax]
    if (nrow(p_sub)) {
      p <- p + geom_rect(
        data = p_sub,
        aes(xmin = start_g, xmax = end_g, ymin = -Inf, ymax = Inf),
        fill = color_phage_hl, color = NA, inherit.aes = FALSE
      )
    }
  }

  add_line <- function(plt, dt) {
    if (is.null(dt) || !nrow(dt)) return(plt)
    dt <- dt[global_pos >= xmin & global_pos <= xmax]
    if (!nrow(dt)) return(plt)
    dt[, Group := factor(Group, levels = group_levels)]
    plt + geom_line(
      data = dt,
      aes(x = global_pos, y = depth, color = Group,
          group = interaction(Group, seg), linewidth = Group),
      lineend = "round", na.rm = TRUE
    )
  }

  p <- add_line(p, phage_dt)
  p <- add_line(p, host_dt)

  if (is.finite(hline_y)) {
    p <- p + geom_hline(yintercept = hline_y, linetype = "dotted", color = "grey20", linewidth = 0.45)
  }

  if (isTRUE(show_labels) && !is.null(phage_regions) && nrow(phage_regions)) {
    p_sub <- phage_regions[end_g >= xmin & start_g <= xmax]
    if (nrow(p_sub)) {
      p <- p + geom_text(
        data = p_sub, aes(x = (start_g + end_g) / 2, y = Inf, label = label),
        vjust = 1.3, size = 3.1, fontface = "bold", family = font, inherit.aes = FALSE
      )
    }
  }

  p +
    scale_color_manual(values = custom_colors, drop = TRUE) +
    scale_linewidth_manual(values = line_widths, guide = "none") +
    labs(title = title, subtitle = subtitle, x = "Global position (bp)", y = y_lab) +
    shared_theme(font, legend = legend) +
    coord_cartesian(xlim = c(xmin, xmax), expand = FALSE)
}

make_density <- function(dt_host, dt_phage, dt_mgs, phage_group_label, x_lab, title, subtitle, font, legend = TRUE) {
  chunks <- list()
  if (!is.null(dt_host) && nrow(dt_host))   chunks[[length(chunks) + 1L]] <- data.table(depth = dt_host$depth, Group = label_host)
  if (!is.null(dt_mgs) && nrow(dt_mgs))     chunks[[length(chunks) + 1L]] <- data.table(depth = dt_mgs$depth, Group = label_mgs)
  if (!is.null(dt_phage) && nrow(dt_phage)) chunks[[length(chunks) + 1L]] <- data.table(depth = dt_phage$depth, Group = phage_group_label)
  if (!length(chunks)) return(ggplot() + theme_void() + ggtitle("Density Data Missing"))

  plot_dt <- rbindlist(chunks)[is.finite(depth)]
  present_levels <- intersect(group_levels, unique(plot_dt$Group))
  plot_dt[, Group := factor(Group, levels = present_levels)]
  plot_dt[, depth := safe_var_fix(depth), by = Group]
  plot_dt[, depth := pmax(depth, 1e-6)]

  ggplot(plot_dt, aes(x = depth, fill = Group)) +
    geom_density(color = "black", alpha = 0.5, linewidth = 0.45, na.rm = TRUE) +
    scale_x_log10(labels = label_scientific()) +
    scale_fill_manual(values = custom_colors, drop = TRUE) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = "Density") +
    shared_theme(font, legend = legend)
}

save_pdf <- function(path, plot, width, height, font) {
  ok <- FALSE
  if (isTRUE(capabilities("cairo"))) {
    ok <- tryCatch({
      ggsave(path, plot, width = width, height = height, device = cairo_pdf, family = font)
      TRUE
    }, error = function(e) FALSE)
  }
  if (!ok) ggsave(path, plot, width = width, height = height, device = "pdf")
}

save_overview <- function(h_src, p_src, off_dt, genome_end,
                          contig_regions, phage_regions, y_lab, title_main, subtitle,
                          hline, out_prefix) {
  h_bin <- prepare_full_track(h_src, off_dt, BIN_SIZE)
  p_bin <- prepare_full_track(p_src, off_dt, BIN_SIZE)
  line_gap <- 1.5 * BIN_SIZE
  bin_note <- sprintf("binned mean %sbp", fmt_int(BIN_SIZE))

  p_ov <- make_profile(
    h_bin, p_bin, 1, genome_end, y_lab,
    title_main, paste(c(subtitle, bin_note), collapse = "  |  "),
    contig_regions, phage_regions, hline, FONT_FAMILY,
    show_labels = TRUE, legend = TRUE, downsample_host = FALSE, line_gap = line_gap
  )
  out_file <- paste0(out_prefix, ".pdf")
  save_pdf(out_file, p_ov, width = 14, height = 4.8, font = FONT_FAMILY)
  message("  Saved: ", out_file)
}

target_pairs <- unique(target_dt[, .(Folder, Genome)])

for (pair_idx in seq_len(nrow(target_pairs))) {
  sample <- target_pairs$Folder[pair_idx]
  genome <- target_pairs$Genome[pair_idx]
  gd <- file.path(map_dir, sample, "depth_to_stat", genome)
  if (!dir.exists(gd)) next

  sub_all_meta <- sum_dt[Folder == sample & Genome == genome]
  target_files <- target_dt[Folder == sample & Genome == genome, Phage_File]

  phage_files <- list.files(gd, pattern = "^phage_depth_.*\\.txt(\\.gz)?$", full.names = TRUE)
  if (!length(phage_files)) next

  dt_host_raw <- read_depth(find_existing(
    file.path(gd, paste0(genome, "_host_nonphage_depth.txt")),
    file.path(gd, paste0(genome, "_host_nonphage_depth.txt.gz"))
  ))
  dt_mgs_raw <- read_depth(find_existing(
    file.path(gd, paste0(genome, "_mgs_depth.txt")),
    file.path(gd, paste0(genome, "_mgs_depth.txt.gz"))
  ))

  phage_list <- list()
  for (pf in phage_files) {
    meta <- parse_phage_filename(pf)
    if (is.null(meta)) next
    dt <- read_depth(pf)
    if (is.null(dt) || !nrow(dt)) next

    row_m <- sub_all_meta[Phage_File == meta$bn]
    pred_val <- if (nrow(row_m)) row_m$prediction[1] else "dormant"
    pval_val <- if (nrow(row_m)) row_m$p_value[1] else NA_real_
    eff_val  <- if (nrow(row_m)) {
      if ("Effect_Size_Value" %in% names(row_m)) row_m$Effect_Size_Value[1]
      else if ("Delta_Mean" %in% names(row_m)) row_m$Delta_Mean[1]
      else NA_real_
    } else NA_real_

    grp_lbl <- pred_to_label(pred_val)
    dt[, Group := grp_lbl]

    phage_list[[length(phage_list) + 1L]] <- list(
      file = pf, bn = meta$bn, contig = meta$contig, start = meta$start, end = meta$end,
      prediction = pred_val, p_value = pval_val, effect = eff_val,
      group_label = grp_lbl, is_target = meta$bn %in% target_files, dt = dt
    )
  }
  if (!length(phage_list)) next

  len_dt <- read_fai_lengths(genome_dir, genome)
  if (is.null(len_dt)) len_dt <- read_bed_lengths(gd, genome)
  if (is.null(len_dt)) len_dt <- infer_lengths(c(list(dt_host_raw, dt_mgs_raw), lapply(phage_list, `[[`, "dt")))

  off_all <- build_offsets(len_dt)
  genome_end_all <- max(off_all$end_global)

  for (i in seq_along(phage_list)) {
    phage_list[[i]]$dt <- add_global_pos(phage_list[[i]]$dt, off_all)
    off <- off_all[contig == phage_list[[i]]$contig, offset]
    phage_list[[i]]$start_g <- if (length(off)) off[1] + phage_list[[i]]$start else phage_list[[i]]$start
    phage_list[[i]]$end_g   <- if (length(off)) off[1] + phage_list[[i]]$end   else phage_list[[i]]$end
  }

  o <- order(vapply(phage_list, `[[`, numeric(1), "start_g"))
  phage_list <- phage_list[o]
  for (i in seq_along(phage_list)) {
    phage_list[[i]]$label  <- paste0("P", i)
    phage_list[[i]]$pretty <- sprintf("%s %s:%s-%s", phage_list[[i]]$label, phage_list[[i]]$contig, phage_list[[i]]$start, phage_list[[i]]$end)
  }

  target_phages <- Filter(function(ph) isTRUE(ph$is_target), phage_list)
  if (!length(target_phages)) next

  tag <- paste(sample, genome, sep = " | ")
  safe_stem <- safe_token(paste(sample, genome, sep = "_"))
  message(">>> [Plot] ", tag, " (plotting ", length(target_phages), " / ", length(phage_list), " prophage(s))")

  scale_dt <- function(dt, s_mode, med) {
    if (is.null(dt) || !nrow(dt)) return(dt)
    out <- copy(dt)
    if (s_mode == "norm" && is.finite(med) && med > 0) out[, depth := depth / med]
    out[]
  }

  if (identical(run_mode, "fast")) {
    if (is.null(dt_mgs_raw) || !nrow(dt_mgs_raw)) next
    ref_median <- stats::median(dt_mgs_raw$depth, na.rm = TRUE)

    for (s_mode in c("raw", "norm")) {
      if (s_mode == "norm" && (!is.finite(ref_median) || ref_median <= 0)) next
      m_sc   <- scale_dt(dt_mgs_raw, s_mode, ref_median)
      dens_x <- if (s_mode == "raw") "Depth (count + 1)" else "Depth / MGS median"
      plus1  <- if (s_mode == "raw") function(d) { if (is.null(d)) return(d); copy(d)[, depth := depth + 1] } else identity

      for (ph in target_phages) {
        sub_txt <- sprintf("Prediction: %s | p = %s | Effect = %s | MGS median = %.3f",
                           ph$prediction, fmt_stat(ph$p_value), fmt_stat(ph$effect), ref_median)
        p_den <- make_density(
          NULL, plus1(scale_dt(ph$dt, s_mode, ref_median)), plus1(m_sc),
          ph$group_label, dens_x,
          paste0("Phage density (", s_mode, "): ", tag, " | ", ph$pretty),
          sub_txt, FONT_FAMILY, legend = TRUE
        )
        phage_file <- file.path(out_dir, sprintf("%s_%s_%s_%s_%s_density_%s.pdf",
                                                 safe_stem, ph$label, safe_token(ph$contig), ph$start, ph$end, s_mode))
        save_pdf(phage_file, p_den, width = 7.5, height = 4.6, font = FONT_FAMILY)
        message("  Saved: ", phage_file)
      }
    }
  } else {
    if (is.null(dt_host_raw) || !nrow(dt_host_raw)) next
    dt_host <- add_global_pos(dt_host_raw, off_all)
    ref_median <- stats::median(dt_host$depth, na.rm = TRUE)
    ref_name <- "host_nonphage"

    regions_all <- rbindlist(lapply(phage_list, function(ph) {
      data.table(start_g = ph$start_g, end_g = ph$end_g, label = ph$label)
    }))
    dt_phage_all <- rbindlist(lapply(phage_list, `[[`, "dt"), use.names = TRUE, fill = TRUE)
    phage_contigs_unique <- unique(vapply(phage_list, `[[`, character(1), "contig"))
    contig_regions_all <- off_all[contig %in% phage_contigs_unique, .(start_g = offset, end_g = end_global, contig)]

    if (isTRUE(overview)) {
      len_keep <- len_dt[is.finite(length) & length >= min_contig]
      if (nrow(len_keep) > 0) {
        off_full <- build_offsets(len_keep)
        genome_end_full <- max(off_full$end_global)
        keep_names <- len_keep$contig

        dt_h_full <- filter_to_contigs(dt_host_raw, keep_names)
        phage_full <- Filter(function(ph) ph$contig %in% keep_names, phage_list)
        dt_p_full <- if (length(phage_full)) rbindlist(lapply(phage_full, function(ph) ph$dt[, .(contig, position, depth, Group)]), use.names = TRUE, fill = TRUE) else NULL

        regions_full <- if (length(phage_full)) {
          rbindlist(lapply(phage_full, function(ph) {
            off <- off_full[contig == ph$contig, offset]
            data.table(start_g = off + ph$start, end_g = off + ph$end, label = ph$label)
          }))
        } else data.table(start_g = numeric(), end_g = numeric(), label = character())

        contig_regs_full <- off_full[contig %in% phage_contigs_unique, .(start_g = offset, end_g = end_global, contig)]

        for (s_mode in c("raw", "norm")) {
          if (s_mode == "norm" && (!is.finite(ref_median) || ref_median <= 0)) next
          y_lab <- if (s_mode == "raw") "Depth" else paste0("Depth / ", ref_name, " median")
          hline <- if (s_mode == "raw") ref_median else 1
          save_overview(
            scale_dt(dt_h_full, s_mode, ref_median), scale_dt(dt_p_full, s_mode, ref_median),
            off_full, genome_end_full, contig_regs_full, regions_full,
            y_lab, paste0("Genome depth (", s_mode, "): ", tag),
            sprintf("%s median=%.3f | contigs>=%sbp", ref_name, ref_median, fmt_int(min_contig)),
            hline, file.path(out_dir, paste0(safe_stem, "_full_", s_mode))
          )
        }
      }

      groups <- build_contig_groups(target_phages, len_dt, CONTIG_EXPAND)
      for (grp in groups) {
        len_sub <- len_dt[contig %in% grp$contigs]
        if (!nrow(len_sub)) next
        off_sub <- build_offsets(len_sub)
        end_sub <- max(off_sub$end_global)

        phage_on <- Filter(function(ph) ph$contig %in% grp$contigs, phage_list)
        p_src_raw <- if (length(phage_on)) rbindlist(lapply(phage_on, function(ph) ph$dt[, .(contig, position, depth, Group)]), use.names = TRUE, fill = TRUE) else NULL

        regions_sub <- if (length(phage_on)) {
          rbindlist(lapply(phage_on, function(ph) {
            off <- off_sub[contig == ph$contig, offset]
            data.table(start_g = off + ph$start, end_g = off + ph$end, label = ph$label)
          }))
        } else data.table(start_g = numeric(), end_g = numeric(), label = character())

        contig_regs_sub <- off_sub[contig %in% grp$center, .(start_g = offset, end_g = end_global, contig)]

        for (s_mode in c("raw", "norm")) {
          if (s_mode == "norm" && (!is.finite(ref_median) || ref_median <= 0)) next
          y_lab <- if (s_mode == "raw") "Depth" else paste0("Depth / ", ref_name, " median")
          hline <- if (s_mode == "raw") ref_median else 1
          file_tag <- paste0(safe_stem, "_", safe_token(grp$center), if (isTRUE(grp$expanded)) "_plusflank" else "", "_bgcontig_", s_mode)
          title_main <- paste0("Contig depth (", s_mode, "): ", tag, " | ", paste(vapply(grp$phages, `[[`, character(1), "label"), collapse = ","))

          save_overview(
            scale_dt(filter_to_contigs(dt_host_raw, grp$contigs), s_mode, ref_median),
            scale_dt(p_src_raw, s_mode, ref_median),
            off_sub, end_sub, contig_regs_sub, regions_sub,
            y_lab, title_main,
            sprintf("%s median=%.3f | Host contig: %s", ref_name, ref_median, grp$center),
            hline, file.path(out_dir, file_tag)
          )
        }
      }
    }

    for (s_mode in c("raw", "norm")) {
      if (s_mode == "norm" && (!is.finite(ref_median) || ref_median <= 0)) next
      h_scaled <- scale_dt(dt_host, s_mode, ref_median)
      p_all_sc <- scale_dt(dt_phage_all, s_mode, ref_median)

      y_lab  <- if (s_mode == "raw") "Depth" else paste0("Depth / ", ref_name, " median")
      dens_x <- if (s_mode == "raw") "Depth (count + 1)" else y_lab
      hline  <- if (s_mode == "raw") ref_median else 1
      plus1  <- if (s_mode == "raw") function(d) { if (is.null(d)) return(d); copy(d)[, depth := depth + 1] } else identity

      for (ph in target_phages) {
        win_min <- max(1, ph$start_g - buffer_bp)
        win_max <- min(genome_end_all, ph$end_g + buffer_bp)
        sub_txt <- sprintf("Prediction: %s | p = %s | Delta = %s | %s median = %.3f",
                           ph$prediction, fmt_stat(ph$p_value), fmt_stat(ph$effect), ref_name, ref_median)

        p_win <- make_profile(
          h_scaled, p_all_sc, win_min, win_max, y_lab,
          paste0(ph$pretty, "  window"), sub_txt,
          contig_regions_all, regions_all, hline, FONT_FAMILY,
          show_labels = FALSE, legend = TRUE, downsample_host = FALSE
        )

        p_den <- make_density(
          plus1(h_scaled), plus1(scale_dt(ph$dt, s_mode, ref_median)), NULL,
          ph$group_label, dens_x, paste0(ph$label, "  density (log10)"), NULL, FONT_FAMILY, legend = TRUE
        )

        combined <- (p_win | p_den) +
          plot_annotation(
            title = paste0("Phage depth (", s_mode, "): ", tag, " | ", ph$pretty),
            theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5, family = FONT_FAMILY))
          )

        phage_file <- file.path(out_dir, sprintf("%s_%s_%s_%s_%s_%s.pdf",
                                                 safe_stem, ph$label, safe_token(ph$contig), ph$start, ph$end, s_mode))
        save_pdf(phage_file, combined, width = 14, height = 4.6, font = FONT_FAMILY)
        message("  Saved: ", phage_file)
      }
    }
  }
}

message("[Plot] Visualization stage finished. Output directory: ", normalizePath(out_dir, mustWork = FALSE))
REOF_PLOT
fi

if [ "$rm_depth" = true ] && [ "$rm_all" = false ]; then
    echo ">>> [--rmdepth] Removing per-base depth (.txt) files..."
    find "$output_map_dir" -name "phage_depth_*.txt" -delete
    find "$output_map_dir" -name "*_mgs_depth.txt" -delete
    find "$output_map_dir" -name "*_host_nonphage_depth.txt" -delete
    find "$output_map_dir" -name ".depth_params.stamp" -delete
fi

if [ "$rm_prodigal" = true ] && [ "$mode" = "fast" ]; then
    echo ">>> [--rmprodigal] Removing Prodigal directory: $output_host_faa"
    rm -rf "$output_host_faa"
fi

if [ "$rm_mgs" = true ] && [ "$mode" = "fast" ]; then
    echo ">>> [--rmmgs] Removing FetchMGs and MGS coordinate intermediate files..."
    rm -rf "$fetchmgs_dir" "$output/mgs_ext_results"
    rm -f  "$output/target_cogs_list.txt" \
           "$output/mgs_missing_genomes.log" \
           "$output/mgs_ordinate.tsv" \
           "$genomes_with_mgs_file"
fi

if [ "$rm_all" = true ]; then
    echo ">>> [--rmall] Removing all alignment and intermediate directories..."
    rm -rf "$output/mapping_results" \
           "$output/competitive_refs" \
           "$output_host_faa" \
           "$fetchmgs_dir" \
           "$output/mgs_ext_results"
    rm -f  "$output/target_cogs_list.txt" \
           "$output/mgs_missing_genomes.log" \
           "$output/mgs_ordinate.tsv" \
           "$genomes_with_mgs_file"
    echo ">>> [--rmall] Cleanup complete. Retained final output directory: $output/activity_results"
fi

echo ">>> Pipeline execution completed successfully!"
