#!/bin/bash

# ==========================================
# 0. Core setup: locate script directory
# ==========================================
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# ==========================================
# 1. Help function and argument parsing
# ==========================================
show_help() {
cat <<'EOF'
ACTIVE - Activity Comparison Test for Induced Viral Elements

Tests whether prophage regions show significantly higher sequencing depth
compared to the bacterial host genome background. Supports hybrid mapping
with Illumina short reads (BBMap) and Nanopore long reads (Minimap2).

Usage:
  bash active.sh -g <host_genome_dir> -s <sequencing_data_dir> \
                         -p <phage_region_file> -o <output_dir>
                         [-t <threads>] [-r <r_threads>] [-i <index_table>]
                         [-m] [--score <fetchmgs_score>] [--cov <coverage>]
                         [--iter <n>]
                         [--rmbam] [--rmdepth] [--rmprodigal] [--rmmgs]
                         [--clean] [--rmall]

Arguments:
  -g   Directory containing host genome(s); files must use .fasta extension
  -s   Directory containing sequencing data samples
         Short reads : *_1_clean.fastq[.gz]  and  *_2_clean.fastq[.gz]
         Long reads  : *_ont.fastq[.gz]  (optional)
  -p   Phage coordinate file (tab-separated: contig  start  end)
  -o   Output directory
  -t   Threads for mapping steps (default: 20)
  -r   Threads for R statistical analysis (default: 20)
  -i   Index table mapping sample names to reference genome names
         (tab-separated, first row = header)
  -m   Run in precise mode (default: fast mode)
  --score
       FetchMGs bit-score threshold; fast mode only (default: 300)
  --cov
       Minimum coverage fraction required for a phage region to be tested
       (default: 0.75; passed to the R statistics step)
  --iter
       Number of permutation iterations for Cliff's delta in precise mode
       (default: 500)
  --rmbam
       Delete BAM files after depth extraction to save disk space
       (default: keep)
  --rmdepth
       Delete per-region depth .txt files after R analysis
       (default: keep)
  --rmprodigal
       Delete Prodigal output directory after pipeline completes
       (fast mode only; default: keep)
  --rmmgs
       Delete all FetchMGs-related files and directories after pipeline
       completes (fast mode only; default: keep)
  --clean
       Shorthand for --rmbam --rmdepth: delete all BAM and depth files
  --rmall
       Delete ALL intermediate files and directories after the pipeline
       completes. Only the final activity_results/<mode>/ folder is kept.
       Equivalent to --rmbam --rmdepth --rmprodigal --rmmgs, and also
       removes the entire mapping_results/ tree (except the CSV summaries
       which have already been moved to activity_results/).
       (default: keep everything)
  -h, --help
       Show this help message and exit

Modes:
  fast      Runs Prodigal + FetchMGs to identify single-copy marker gene
            (MGS) regions as the host background for statistical comparison.
  precise   Skips Prodigal/FetchMGs; uses all non-phage genomic positions
            as the host background.

Output:
  <output>/activity_results/<mode>/   Final per-sample CSV summary files
  <output>/mapping_results/           BAM files and per-region depth tables
  <output>/prodigal_results/          Prodigal GFF/FAA outputs (fast only)
  <output>/fetchmgs_results/          FetchMGs scores (fast only)
  <output>/mgs_ext_results/           MGS coordinate TSV (fast only)

EOF
}

# ---- Default values ----
threads=20
r_threads=20
index_table=""
skip_all_steps=false
mgs_score=300
coverage_threshold=0.75
n_iter=500
rm_bam=false
rm_depth=false
rm_prodigal=false
rm_mgs=false
rm_all=false

# ---- Pre-process long options before getopts ----
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cov)         coverage_threshold="$2"; shift 2 ;;
        --cov=*)       coverage_threshold="${1#*=}"; shift ;;
        --score)       mgs_score="$2"; shift 2 ;;
        --score=*)     mgs_score="${1#*=}"; shift ;;
        --iter)        n_iter="$2"; shift 2 ;;
        --iter=*)      n_iter="${1#*=}"; shift ;;
        --rmbam)       rm_bam=true; shift ;;
        --rmdepth)     rm_depth=true; shift ;;
        --rmprodigal)  rm_prodigal=true; shift ;;
        --rmmgs)       rm_mgs=true; shift ;;
        --clean)       rm_bam=true; rm_depth=true; shift ;;
        --rmall)       rm_all=true; shift ;;
        --help)        show_help; exit 0 ;;
        *)             ARGS+=("$1"); shift ;;
    esac
done
set -- "${ARGS[@]}"

while getopts "hg:s:p:o:t:r:i:m" opt; do
    case $opt in
        h) show_help; exit 0 ;;
        g) in_hostgenome="$OPTARG" ;;
        s) in_seq="$OPTARG" ;;
        p) phage_region="$OPTARG" ;;
        o) output="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        r) r_threads="$OPTARG" ;;
        i) index_table="$OPTARG" ;;
        m) skip_all_steps=true ;;
        *) echo "Unknown argument. Use -h to see help."; exit 1 ;;
    esac
done

# ---- Validate required arguments ----
if [ -z "$in_hostgenome" ] || [ -z "$in_seq" ] || [ -z "$phage_region" ] || [ -z "$output" ]; then
    echo "Error: Missing required arguments. Use -h to see help."
    exit 1
fi

if [ "$skip_all_steps" = true ]; then
    mode="precise"
else
    mode="fast"
fi

# --rmall implies all individual cleanup flags
if [ "$rm_all" = true ]; then
    rm_bam=true
    rm_depth=true
    rm_prodigal=true
    rm_mgs=true
fi

echo "Config: threads=$threads, r_threads=$r_threads, mode=$mode, cov=$coverage_threshold, iter=$n_iter"
[ "$mode"        == "fast" ] && echo ">>> Fast mode: FetchMGs score threshold = $mgs_score"
[ "$rm_all"      = true ]   && echo ">>> --rmall: all intermediate files will be deleted after pipeline completes."
[ "$rm_bam"      = true ] && [ "$rm_all" = false ] && echo ">>> --rmbam: BAM files will be deleted after depth extraction."
[ "$rm_depth"    = true ] && [ "$rm_all" = false ] && echo ">>> --rmdepth: Depth files will be deleted after R analysis."
[ "$rm_prodigal" = true ] && [ "$rm_all" = false ] && echo ">>> --rmprodigal: Prodigal output will be deleted at pipeline end."
[ "$rm_mgs"      = true ] && [ "$rm_all" = false ] && echo ">>> --rmmgs: FetchMGs-related files will be deleted at pipeline end."

mkdir -p "$output"

# ==========================================
# 2. Read index table (if provided)
# ==========================================
declare -A sample_to_ref
if [ -n "$index_table" ] && [ -f "$index_table" ]; then
    echo ">>> Reading index table: $index_table"
    while IFS=$'\t' read -r sample ref; do
        if [ -z "${sample_to_ref["$sample"]}" ]; then
            sample_to_ref["$sample"]="$ref"
        else
            sample_to_ref["$sample"]="${sample_to_ref["$sample"]},$ref"
        fi
    done < <(tail -n +2 "$index_table")
else
    echo ">>> No index table provided; sample name is assumed to match reference genome name."
fi

# ==========================================
# 3. Host genome preprocessing (Prodigal & FetchMGs)  [fast mode only]
# ==========================================
output_host_faa="$output/prodigal_results"
fetchmgs_dir="$output/fetchmgs_results"

if [ "$mode" == "fast" ]; then
    echo ">>> [Step 1] Running host CDS prediction (Prodigal)..."
    mkdir -p "$output_host_faa"

    for hostgenome in "$in_hostgenome"/*.fasta; do
        [ -e "$hostgenome" ] || continue
        host_name=$(basename "$hostgenome" .fasta)
        outdir="$output_host_faa/$host_name"
        mkdir -p "$outdir"
        target_gff="$outdir/$host_name.gff"
        target_faa="$outdir/$host_name.faa"
        if [ -s "$target_gff" ] && [ -s "$target_faa" ]; then
            echo "    [SKIP] Prodigal output already exists: $host_name"
        else
            echo "    [RUN]  Prodigal: $host_name"
            prodigal -i "$hostgenome" -o "$target_gff" -f gff -a "$target_faa" -p meta -q >/dev/null 2>&1
        fi
    done

    echo ">>> [Step 2] Extracting single-copy marker genes (FetchMGs)..."
    mkdir -p "$fetchmgs_dir"
    FETCHMGS_THREADS=28   # Fixed at 28; not exposed as a user parameter

    find "$output_host_faa" -name "*.faa" | while read -r hostfaa; do
        host_base=$(basename "${hostfaa%.faa}")
        target_scores_dir="$fetchmgs_dir/$host_base"
        target_scores_file="$fetchmgs_dir/$host_base/${host_base}.faa.fetchMGs.scores"
        if [ -s "$fetchmgs_dir/$host_base.faa.fetchMGs.scores" ]; then
            echo "    [SKIP] FetchMGs output already exists: $host_base"
        elif [ -d "$target_scores_dir" ] && [ -s "$target_scores_file" ]; then
            echo "    [SKIP] FetchMGs directory already exists: $host_base"
        else
            echo "    [RUN]  FetchMGs: $host_base (threads=$FETCHMGS_THREADS)"
            fetchMGs extraction "$hostfaa" gene "$fetchmgs_dir/$host_base" -t "$FETCHMGS_THREADS" >/dev/null 2>&1
        fi
    done

    # Consolidate scores files to the top-level fetchmgs directory
    find "$fetchmgs_dir" -mindepth 2 -name "*fetchMGs.scores" -exec mv {} "$fetchmgs_dir/" \;

    echo "    Generating COG list..."
    cog_list_file="$output/target_cogs_list.txt"
    awk -v score="$mgs_score" \
    '!/^#/ && $2 > score {
        if (match($3, /[0-9]{4}/))
            print "COG" substr($3, RSTART, RLENGTH)
    }' "$fetchmgs_dir"/*.scores 2>/dev/null \
    | sort -u > "$cog_list_file"

    if [ ! -s "$cog_list_file" ]; then
        echo "    Warning: No COGs found with score > $mgs_score."
    else
        cog_count=$(wc -l < "$cog_list_file")
        echo "    MGS COG list length: $cog_count"
        > "$output/mgs_missing_genomes.log"

        # ------------------------------------------------------------------
        # MGS coordinate extraction: Python code runs inline via stdin.
        # No .py file is written to disk.
        # ------------------------------------------------------------------
        echo "    [RUN]  Python (inline): extracting MGS coordinates..."

        python3 - \
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
    if n.endswith(".gff"): return n[:-4]
    if n.endswith(".faa"): return n[:-4]
    return p.stem

def discover_pairs(base):
    gffs, faas = {}, {}
    for q in base.rglob("*"):
        if not q.is_file(): continue
        if q.name.endswith(".gff"):   gffs[base_stem(q)] = q
        elif q.name.endswith(".faa"): faas[base_stem(q)] = q
    pairs = {k:(gffs[k], faas[k]) for k in gffs.keys() & faas.keys()}
    eprint(f"[INFO] Paired samples: {len(pairs)}")
    for s in sorted(set(faas)-set(gffs))[:5]: eprint(f"[WARN] Missing GFF for: {s}")
    for s in sorted(set(gffs)-set(faas))[:5]: eprint(f"[WARN] Missing FAA for: {s}")
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

def load_scores(p):
    files = [p] if p.is_file() else list(p.rglob("*.scores")) if p.is_dir() else []
    m = {}
    for sp in files:
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
    else:
        for tok in re.split(r"[,\s;]+", s.strip()):
            if tok: vals.append(tok)
    out = set()
    for v in vals:
        mm = re.search(r"(\d{4})", v)
        if mm: out.add(f"COG{mm.group(1)}")
    return out if out else None

def join_namefree(by_id, pid2gid, scoremap, targets):
    out = []
    for pid, (cog, score) in scoremap.items():
        gid = pid2gid.get(pid)
        if not gid: continue
        rec = by_id.get(gid)
        if not rec: continue
        if targets and cog not in targets: continue
        out.append({
            "contig": rec["contig"], "start": rec["start"],
            "end": rec["end"], "strand": rec["strand"],
            "locus_tag": rec.get("locus_tag") or gid,
            "gene": rec.get("gene", ""), "product": rec.get("product", ""),
            "COG": cog, "bit_score": score,
        })
    return out

ap = argparse.ArgumentParser()
ap.add_argument("--input",   required=True)
ap.add_argument("--scores",  required=True)
ap.add_argument("--targets", default=None)
ap.add_argument("--outdir",  default="cog_coords_min_out")
args = ap.parse_args()

base     = Path(args.input)
pairs    = discover_pairs(base)
if not pairs: eprint("ERROR: No (gff,faa) pairs found."); sys.exit(2)

scoremap = load_scores(Path(args.scores))
if not scoremap: eprint("ERROR: No scores loaded."); sys.exit(3)
eprint(f"[INFO] Unique protein IDs in scores: {len(scoremap)}")

targets = read_targets(args.targets)
if targets: eprint(f"[INFO] Target COG count: {len(targets)}")

merged = []
for name, (gff, faa) in pairs.items():
    by_id   = parse_gff_by_id(gff)
    pid2gid = build_pid2gid(faa)
    # Restrict to this sample's proteins only to prevent cross-sample contamination
    sub     = {p: scoremap[p] for p in pid2gid if p in scoremap}
    rows    = join_namefree(by_id, pid2gid, sub, targets)
    for r in rows: r["genome"] = name
    eprint(f"[INFO] {name}: FAA={len(pid2gid)}, intersect={len(sub)}, hits={len(rows)}")
    merged.extend(rows)

outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
out = outdir / "all_from_scores.tsv"
with out.open("wt", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["genome","contig","start","end","strand",
                "locus_tag","gene","product","COG","bit_score"])
    for r in merged:
        w.writerow([r.get("genome",""), r["contig"], r["start"], r["end"],
                    r["strand"], r["locus_tag"], r["gene"], r["product"],
                    r["COG"], r["bit_score"]])
eprint(f"[DONE] Written: {out}  rows={len(merged)}")
PYEOF

        if [ -s "$output/mgs_ext_results/all_from_scores.tsv" ]; then
            awk -F'\t' 'NR>1 {print $2 "\t" $3 "\t" $4}' \
                "$output/mgs_ext_results/all_from_scores.tsv" > "$output/mgs_ordinate.tsv"
            echo "    MGS coordinate file created: $output/mgs_ordinate.tsv"
        else
            echo "    Error: MGS result file is empty."
            exit 1
        fi
    fi
else
    echo ">>> [Precise Mode] Skipping Prodigal and FetchMGs steps."
fi

# ==========================================
# 4. Mapping and depth extraction (hybrid mapping)
# ==========================================
echo ">>> [Step 3] Starting hybrid mapping and depth extraction..."

make_bed_from_tsv() {
  awk -v OFS='\t' 'NF>=3{ s=$2; e=$3; if(s<1)s=1; if(e>=s) print $1, s-1, e }' "$1" \
  | sort -k1,1 -k2,2n | bedtools merge -i - > "$2"
}
make_nonphage_bed() { bedtools subtract -a "$1" -b "$2" > "$3"; }

output_map_dir="$output/mapping_results"
mkdir -p "$output_map_dir"

for data1 in "$in_seq"/*_1_clean.fastq*; do
  [ -e "$data1" ] || continue

  if [[ "$data1" == *.gz ]]; then
    base_name=$(basename "$data1" _1_clean.fastq.gz)
    data2="$in_seq/${base_name}_2_clean.fastq.gz"
  else
    base_name=$(basename "$data1" _1_clean.fastq)
    data2="$in_seq/${base_name}_2_clean.fastq"
  fi

  if [ ! -f "$data2" ]; then
    echo "    Warning: Missing paired-end file $data2, skipping $base_name"; continue
  fi

  # Detect optional Nanopore long-read data
  ont_data=""
  [ -f "$in_seq/${base_name}_ont.fastq.gz" ] && ont_data="$in_seq/${base_name}_ont.fastq.gz"
  [ -z "$ont_data" ] && [ -f "$in_seq/${base_name}_ont.fastq" ] && ont_data="$in_seq/${base_name}_ont.fastq"

  if [ -n "$ont_data" ]; then
      echo "    [Hybrid mode]     Sample $base_name: ONT data detected: $(basename "$ont_data")"
  else
      echo "    [Short-read mode] Sample $base_name: no ONT data found."
  fi

  if [ -n "$index_table" ] && [ ${#sample_to_ref[@]} -gt 0 ]; then
    ref_str="${sample_to_ref[$base_name]}"
    [ -z "$ref_str" ] && continue
  else
    ref_str="$base_name"
  fi

  sample_dir="${output_map_dir}/${base_name}"
  mkdir -p "$sample_dir"

  IFS=',' read -ra REF_GENOMES <<< "$ref_str"
  for ref_name in "${REF_GENOMES[@]}"; do
      ref_file="$in_hostgenome/${ref_name}.fasta"
      [ ! -f "$ref_file" ] && continue

      bam_file="${sample_dir}/${base_name}_${ref_name}.bam"

      if [ -s "$bam_file" ]; then
          echo "    [SKIP] BAM already exists: $bam_file"
      else
          echo "    [Mapping] $base_name -> $ref_name"
          temp_sr_bam="${sample_dir}/temp_sr_${base_name}_${ref_name}.bam"
          temp_lr_bam="${sample_dir}/temp_lr_${base_name}_${ref_name}.bam"

          # [A] Short-read alignment (BBMap): ambiguous=random, minid=0.95
          bbmap.sh in1="$data1" in2="$data2" ref="$ref_file" out=stdout.sam threads=$threads \
              ambiguous=random minid=0.95 nodisk=t \
          | samtools view -Sb -@ "$threads" - \
          | samtools sort -o "$temp_sr_bam" -@ "$threads" -

          # [B] Long-read alignment (Minimap2 ONT preset) — only when ONT data exists
          has_lr_bam=false
          if [ -n "$ont_data" ]; then
              echo "        -> Running Minimap2 (ONT mode)..."
              minimap2 -ax map-ont -t "$threads" "$ref_file" "$ont_data" \
              | samtools view -Sb -@ "$threads" - \
              | samtools sort -o "$temp_lr_bam" -@ "$threads" -
              [ -s "$temp_lr_bam" ] && has_lr_bam=true
          fi

          # [C] Merge BAMs or rename
          if [ "$has_lr_bam" = true ]; then
              echo "        -> Merging short-read and long-read BAMs..."
              samtools merge -@ "$threads" -f "$bam_file" "$temp_sr_bam" "$temp_lr_bam"
              rm -f "$temp_sr_bam" "$temp_lr_bam"
          else
              mv "$temp_sr_bam" "$bam_file"
          fi
          samtools index "$bam_file"
      fi

      [ -f "${ref_file}.fai" ] || samtools faidx "$ref_file"

      # ---- Depth extraction ----
      global_dir="$sample_dir/global"
      mkdir -p "$global_dir"
      grep ">" "$ref_file" | sed 's/>//' | cut -d " " -f1 > "$global_dir/contigs.txt"
      cut -d "_" -f1 "$global_dir/contigs.txt" | sort | uniq > "$global_dir/genomes_prefix.txt"

      while read -r prefix; do
        grep "^${prefix}_" "$phage_region" > "$global_dir/${prefix}_phage_region" 2>/dev/null
        if [ "$mode" == "fast" ]; then
             grep "^${prefix}_" "$output/mgs_ordinate.tsv" > "$global_dir/${prefix}_mgs_region" 2>/dev/null
        fi
      done < "$global_dir/genomes_prefix.txt"

      while read -r prefix; do
        genome_depth_dir="${sample_dir}/depth_to_stat/${prefix}"
        mkdir -p "$genome_depth_dir"
        phage_region_sub="$global_dir/${prefix}_phage_region"

        if [ -s "$phage_region_sub" ]; then
            phage_bed="$genome_depth_dir/${prefix}.phage.bed"
            make_bed_from_tsv "$phage_region_sub" "$phage_bed"

            # 1. Per-region phage depth
            while IFS=$'\t' read -r contig start end; do
                [ -z "$contig" ] && continue
                out_d="$genome_depth_dir/phage_depth_${contig}_${start}_${end}.txt"
                [ -s "$out_d" ] && continue
                samtools depth -@ "$threads" -a -r "${contig}:${start}-${end}" "$bam_file" \
                | awk -v OFS='\t' '{print $1,$2,$3}' > "$out_d"
            done < "$phage_region_sub"

            # 2. Host background depth
            if [ "$mode" == "fast" ]; then
                mgs_region_sub="$global_dir/${prefix}_mgs_region"
                out_mgs="$genome_depth_dir/${prefix}_mgs_depth.txt"
                if [ -s "$mgs_region_sub" ] && [ ! -s "$out_mgs" ]; then
                    mgs_bed="$genome_depth_dir/${prefix}.mgs.bed"
                    make_bed_from_tsv "$mgs_region_sub" "$mgs_bed"
                    samtools depth -@ "$threads" -a -b "$mgs_bed" "$bam_file" \
                    | awk -v OFS='\t' '{print $1,$2,$3}' > "$out_mgs"
                fi
            else
                # Precise mode: all non-phage positions as background
                all_bed="$genome_depth_dir/${prefix}.all.bed"
                nonphage_bed="$genome_depth_dir/${prefix}.nonphage.bed"
                out_host="$genome_depth_dir/${prefix}_host_nonphage_depth.txt"
                if [ ! -s "$out_host" ]; then
                    awk -v OFS='\t' '{print $1,0,$2}' "${ref_file}.fai" \
                    | grep "^${prefix}_" | sort -k1,1 -k2,2n > "$all_bed"
                    make_nonphage_bed "$all_bed" "$phage_bed" "$nonphage_bed"
                    samtools depth -@ "$threads" -a -b "$nonphage_bed" "$bam_file" \
                    | awk -v OFS='\t' '{print $1,$2,$3}' > "$out_host"
                fi
            fi
        fi
      done < "$global_dir/genomes_prefix.txt"

      # ---- Optional: delete BAM after depth extraction ----
      if [ "$rm_bam" = true ]; then
          echo "    [--rmbam] Removing BAM: $bam_file"
          rm -f "$bam_file" "${bam_file}.bai"
      fi
  done
done
echo ">>> Depth extraction stage complete."

# ==========================================
# 5. R statistical analysis (inline via Rscript -)
# ==========================================
echo ">>> [Step 4] Running R statistical analysis (mode: $mode)..."

if [ "$mode" == "fast" ]; then

# ---------- fast mode: R code runs inline via stdin ----------
Rscript - \
    --threads "$r_threads" \
    --mapdir  "$output_map_dir" \
    --coverage "$coverage_threshold" \
<<'REOF_FAST'
suppressPackageStartupMessages({
  library(tseries); library(brunnermunzel)
  library(effsize);  library(doParallel); library(foreach)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse --threads
cores <- {
  hit <- grep("^--threads=", args, value = TRUE)
  if (length(hit) > 0) as.integer(sub("^--threads=", "", hit[1]))
  else { idx <- which(args == "--threads")
    if (length(idx) > 0 && length(args) >= idx+1) as.integer(args[idx+1]) else 22L }
}

# Parse --coverage (minimum phage region coverage fraction; default 0.75)
coverage_threshold <- {
  hit <- grep("^--coverage=", args, value = TRUE)
  if (length(hit) > 0) {
    val <- suppressWarnings(as.numeric(sub("^--coverage=", "", hit[1])))
    if (!is.na(val)) val else 0.75
  } else {
    idx <- which(args == "--coverage")
    if (length(idx) > 0 && length(args) >= idx+1) {
      val <- suppressWarnings(as.numeric(args[idx+1]))
      if (!is.na(val)) val else 0.75
    } else 0.75
  }
}

# Parse --mapdir (required)
root_dir <- {
  hit <- grep("^--mapdir=", args, value = TRUE)
  if (length(hit) > 0) sub("^--mapdir=", "", hit[1])
  else { idx <- which(args == "--mapdir")
    if (length(idx) > 0 && length(args) >= idx+1) args[idx+1] else "" }
}
if (!nzchar(root_dir))     stop("Missing --mapdir", call. = FALSE)
if (!dir.exists(root_dir)) stop(sprintf("mapdir not found: %s", root_dir), call. = FALSE)

cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
message(sprintf("Parallel backend: %d cores | coverage threshold: %.4f", cores, coverage_threshold))

sample_dirs <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[!grepl("global$", basename(sample_dirs))]
message(sprintf("Sample directories: %d", length(sample_dirs)))
alpha <- 0.01

for (sample_dir in sample_dirs) {
  folder_name <- basename(sample_dir)
  message(sprintf("======== %s ========", folder_name))
  stat_main_dir <- file.path(sample_dir, "depth_to_stat")
  if (!dir.exists(stat_main_dir)) { message("  No depth_to_stat; skipping."); next }
  target_dirs <- list.dirs(stat_main_dir, full.names = TRUE, recursive = FALSE)

  results <- data.frame(
    Folder=character(), Genome=character(), Compare=character(),
    Phage_File=character(), Coverage=numeric(), Test_Type=character(),
    p_value=numeric(), Significant=character(), Effect_Size_Type=character(),
    Effect_Size_Value=numeric(), BM_Prob_Superiority=numeric(), Phage_Sites=numeric(),
    stringsAsFactors=FALSE)
  stat_results <- data.frame(
    Folder=character(), Genome=character(), Region=character(),
    File=character(), Mean=numeric(), Median=numeric(), Variance=numeric(), COVERAGE=numeric(),
    stringsAsFactors=FALSE)

  for (datadir in target_dirs) {
    genomename <- basename(datadir)
    message(sprintf("---- Genome: %s", genomename))
    gyrB_file       <- list.files(datadir, pattern="_mgs_depth\\.txt$",    full.names=TRUE)
    phage_files_raw <- list.files(datadir, pattern="^phage_depth.*\\.txt$", full.names=TRUE)
    if (length(gyrB_file)==0)       { message("  Missing MGS file; skipping.");       next }
    gi <- file.info(gyrB_file[1])
    if (is.na(gi$size)||gi$size==0) { message("  MGS file empty; skipping.");         next }
    if (length(phage_files_raw)==0) { message("  Missing phage files; skipping.");    next }
    pfi <- file.info(phage_files_raw)
    phage_files <- phage_files_raw[pfi$size > 0]
    if (length(phage_files)==0)     { message("  All phage files empty; skipping.");  next }

    gyrB <- tryCatch(read.table(gyrB_file[1], header=FALSE), error=function(e) NULL)
    if (is.null(gyrB)||ncol(gyrB)<3) { message("  MGS format error; skipping."); next }
    colnames(gyrB) <- c("id","pos","depth")
    gyrB_depth <- gyrB$depth
    stat_results <- rbind(stat_results, data.frame(
      Folder=folder_name, Genome=genomename, Region="host_mgs",
      File=basename(gyrB_file[1]), Mean=mean(gyrB_depth),
      Median=median(gyrB_depth), Variance=var(gyrB_depth), COVERAGE=NA))

    # Parallel processing: one task per phage region file
    plist <- foreach(f=phage_files,
      .packages=c("tseries","brunnermunzel","effsize"), .errorhandling="pass") %dopar% {
      phage <- tryCatch(read.table(f, header=FALSE), error=function(e) NULL)
      if (is.null(phage)||ncol(phage)<3) return(NULL)
      colnames(phage) <- c("id","pos","depth")
      phage_depth <- phage$depth
      phage_sites <- length(phage_depth)
      coverage <- if (phage_sites>0) sum(phage_depth>0)/phage_sites else 0
      stat_line <- data.frame(
        Folder=folder_name, Genome=genomename, Region="phage", File=basename(f),
        Mean=mean(phage_depth), Median=median(phage_depth),
        Variance=var(phage_depth), COVERAGE=round(coverage,4))
      if (coverage < coverage_threshold) {
        return(list(stat_line=stat_line, result_line=data.frame(
          Folder=folder_name, Genome=genomename, Compare="phage vs mgs",
          Phage_File=basename(f), Coverage=round(coverage,4),
          Test_Type=NA, p_value=NA, Significant="Coverage failed",
          Effect_Size_Type=NA, Effect_Size_Value=NA,
          BM_Prob_Superiority=NA, Phage_Sites=phage_sites)))
      }
      jb1  <- tryCatch(tseries::jarque.bera.test(phage_depth), error=function(e) NULL)
      jb2  <- tryCatch(tseries::jarque.bera.test(gyrB_depth),  error=function(e) NULL)
      norm1 <- !is.null(jb1)&&!is.na(jb1$p.value)&&jb1$p.value>0.05
      norm2 <- !is.null(jb2)&&!is.na(jb2$p.value)&&jb2$p.value>0.05
      test_type<-NA; pval<-NA; significant<-"No"
      effect_type<-NA; effect_value<-NA; bm_prob<-NA
      if (norm1 & norm2) {
        vtest <- var.test(phage_depth, gyrB_depth)
        if (vtest$p.value>0.05) {
          ttest <- t.test(phage_depth,gyrB_depth,alternative="greater",var.equal=TRUE)
          test_type <- "Student's t-test"
        } else {
          ttest <- t.test(phage_depth,gyrB_depth,alternative="greater",var.equal=FALSE)
          test_type <- "Welch's t-test"
        }
        pval <- ttest$p.value
        if (!is.na(pval)&&pval<alpha) {
          significant  <- "Yes"
          cd <- effsize::cohen.d(phage_depth,gyrB_depth,hedges.correction=TRUE)
          effect_type  <- "Cohen's D"; effect_value <- cd$estimate
        }
      } else {
        bmtest    <- brunnermunzel::brunnermunzel.test(phage_depth,gyrB_depth,alternative="greater")
        test_type <- "Brunner-Munzel"; pval <- bmtest$p.value
        bm_prob   <- suppressWarnings(as.numeric(bmtest$estimate))
        if (!is.na(pval)&&pval<alpha) {
          significant <- "Yes"
          cd <- tryCatch(suppressWarnings(effsize::cliff.delta(phage_depth,gyrB_depth)),
                         error=function(e) NULL)
          if (!is.null(cd)) { effect_type <- "Cliff's Delta"; effect_value <- as.numeric(cd$estimate) }
        }
      }
      list(stat_line=stat_line, result_line=data.frame(
        Folder=folder_name, Genome=genomename, Compare="phage vs mgs",
        Phage_File=basename(f), Coverage=round(coverage,4),
        Test_Type=test_type, p_value=pval, Significant=significant,
        Effect_Size_Type=effect_type, Effect_Size_Value=effect_value,
        BM_Prob_Superiority=bm_prob, Phage_Sites=phage_sites))
    }
    plist <- Filter(function(x) is.list(x)&&all(c("stat_line","result_line")%in%names(x)), plist)
    if (length(plist)>0) {
      stat_results <- rbind(stat_results, do.call(rbind, lapply(plist,`[[`,"stat_line")))
      results      <- rbind(results,      do.call(rbind, lapply(plist,`[[`,"result_line")))
    }
  }
  out1 <- file.path(sample_dir, paste0(folder_name,"_phage_depth_stats_summary_final.csv"))
  out2 <- file.path(sample_dir, paste0(folder_name,"_phage_basic_stats_summary.csv"))
  write.csv(results,      out1, row.names=FALSE)
  write.csv(stat_results, out2, row.names=FALSE)
  message(sprintf("Written: %s | %s", out1, out2))
}
message("All samples processed.")
REOF_FAST

else

# ---------- precise mode: R code runs inline via stdin ----------
Rscript - \
    --threads "$r_threads" \
    --mapdir  "$output_map_dir" \
    --coverage "$coverage_threshold" \
    --iter "$n_iter" \
<<'REOF_PRECISE'
suppressPackageStartupMessages({
  library(tseries); library(brunnermunzel)
  library(effsize);  library(doParallel); library(foreach)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse --threads
cores <- {
  hit <- grep("^--threads=", args, value = TRUE)
  if (length(hit) > 0) as.integer(sub("^--threads=", "", hit[1]))
  else { idx <- which(args == "--threads")
    if (length(idx) > 0 && length(args) >= idx+1) as.integer(args[idx+1]) else 22L }
}

# Parse --coverage (minimum phage region coverage fraction; default 0.75)
coverage_threshold <- {
  hit <- grep("^--coverage=", args, value = TRUE)
  if (length(hit) > 0) {
    val <- suppressWarnings(as.numeric(sub("^--coverage=", "", hit[1])))
    if (!is.na(val)) val else 0.75
  } else {
    idx <- which(args == "--coverage")
    if (length(idx) > 0 && length(args) >= idx+1) {
      val <- suppressWarnings(as.numeric(args[idx+1]))
      if (!is.na(val)) val else 0.75
    } else 0.75
  }
}

# Parse --iter (number of permutation iterations for Cliff's delta; default 500)
n_iter <- {
  hit <- grep("^--iter=", args, value = TRUE)
  if (length(hit) > 0) {
    val <- suppressWarnings(as.integer(sub("^--iter=", "", hit[1])))
    if (!is.na(val) && val > 0L) val else 500L
  } else {
    idx <- which(args == "--iter")
    if (length(idx) > 0 && length(args) >= idx+1) {
      val <- suppressWarnings(as.integer(args[idx+1]))
      if (!is.na(val) && val > 0L) val else 500L
    } else 500L
  }
}

# Parse --mapdir (required)
root_dir <- {
  hit <- grep("^--mapdir=", args, value = TRUE)
  if (length(hit) > 0) sub("^--mapdir=", "", hit[1])
  else { idx <- which(args == "--mapdir")
    if (length(idx) > 0 && length(args) >= idx+1) args[idx+1] else "" }
}
if (!nzchar(root_dir))     stop("Missing --mapdir", call. = FALSE)
if (!dir.exists(root_dir)) stop(sprintf("mapdir not found: %s", root_dir), call. = FALSE)

alpha <- 0.01
cl    <- makeCluster(cores)
on.exit(stopCluster(cl), add = TRUE)
registerDoParallel(cl)
message(sprintf("Parallel threads: %d | coverage threshold: %.4f | iter: %d",
                cores, coverage_threshold, n_iter))

sample_dirs <- list.dirs(root_dir, full.names=TRUE, recursive=FALSE)
sample_dirs <- sample_dirs[!grepl("global$", basename(sample_dirs))]
message(sprintf("Sample directories: %d", length(sample_dirs)))

for (s_idx in seq_along(sample_dirs)) {
  sample_dir  <- sample_dirs[s_idx]
  folder_name <- basename(sample_dir)
  cat(sprintf("\n======== [%d/%d] %s ========\n", s_idx, length(sample_dirs), folder_name))
  stat_main_dir <- file.path(sample_dir, "depth_to_stat")
  if (!dir.exists(stat_main_dir)) { cat("  No depth_to_stat; skipping.\n"); next }
  target_dirs <- list.dirs(stat_main_dir, full.names=TRUE, recursive=FALSE)
  cat(sprintf("  Genome dirs: %d\n", length(target_dirs)))

  results <- data.frame(
    Folder=character(), Genome=character(), Compare=character(),
    Phage_File=character(), Coverage=numeric(), Test_Type=character(),
    p_value=numeric(), Significant=character(), BM_Prob_Superiority=numeric(),
    Q95=numeric(), Q99=numeric(), Phage_Sites=numeric(),
    Phage_OverQ95=numeric(), Phage_OverQ99=numeric(),
    RR_Q95=numeric(), RR_Q99=numeric(),
    Delta_Mean=numeric(), Delta_CI_2.5=numeric(), Delta_CI_97.5=numeric(),
    stringsAsFactors=FALSE)
  stat_results <- data.frame(
    Folder=character(), Genome=character(), Region=character(),
    File=character(), Mean=numeric(), Median=numeric(), Variance=numeric(), COVERAGE=numeric(),
    stringsAsFactors=FALSE)

  for (g_idx in seq_along(target_dirs)) {
    datadir    <- target_dirs[g_idx]
    genomename <- basename(datadir)
    cat(sprintf("\n---- [%d/%d] Genome: %s ----\n", g_idx, length(target_dirs), genomename))
    host_file   <- list.files(datadir, pattern="_host_nonphage_depth\\.txt$", full.names=TRUE)
    phage_files <- list.files(datadir, pattern="^phage_depth.*\\.txt$",        full.names=TRUE)
    cat(sprintf("  Host files: %d | Phage region files: %d\n", length(host_file), length(phage_files)))
    if (length(host_file)==0||length(phage_files)==0) {
      cat("  Missing files; skipping genome.\n"); next }

    host <- tryCatch(read.table(host_file[1], header=FALSE), error=function(e) NULL)
    if (is.null(host)||ncol(host)<3) { cat("  Host file format error; skipping.\n"); next }
    colnames(host) <- c("id","pos","depth")
    host_depth <- as.numeric(host$depth); host_depth <- host_depth[is.finite(host_depth)]
    if (length(host_depth)==0) { cat("  No valid host depth; skipping.\n"); next }
    nonphage_sites <- length(host_depth)
    stat_results <- rbind(stat_results, data.frame(
      Folder=folder_name, Genome=genomename, Region="host_non_phage",
      File=basename(host_file[1]), Mean=mean(host_depth),
      Median=median(host_depth), Variance=stats::var(host_depth), COVERAGE=NA_real_))

    # Global Q95/Q99 thresholds from all depth values (phage + host combined)
    phage_all <- unlist(lapply(phage_files, function(f) {
      df <- tryCatch(read.table(f, header=FALSE), error=function(e) NULL)
      if (is.null(df)||ncol(df)<3) return(numeric(0)); as.numeric(df[[3]])
    }))
    phage_all <- phage_all[is.finite(phage_all)]
    all_d <- c(phage_all, host_depth)
    if (length(all_d)==0) { cat("  No valid depth; skipping.\n"); next }
    q95 <- as.numeric(quantile(all_d, 0.95, names=FALSE, type=7))
    q99 <- as.numeric(quantile(all_d, 0.99, names=FALSE, type=7))
    p0_95 <- sum(host_depth > q95) / nonphage_sites
    p0_99 <- sum(host_depth > q99) / nonphage_sites

    for (f in phage_files) {
      phage <- tryCatch(read.table(f, header=FALSE), error=function(e) NULL)
      if (is.null(phage)||ncol(phage)<3) next
      colnames(phage) <- c("id","pos","depth")
      phage_depth <- as.numeric(phage$depth); phage_depth <- phage_depth[is.finite(phage_depth)]
      phage_sites <- length(phage_depth)
      if (phage_sites==0) next
      coverage <- sum(phage_depth>0) / phage_sites
      stat_results <- rbind(stat_results, data.frame(
        Folder=folder_name, Genome=genomename, Region="phage", File=basename(f),
        Mean=mean(phage_depth), Median=median(phage_depth),
        Variance=stats::var(phage_depth), COVERAGE=round(coverage,4)))

      if (coverage < coverage_threshold) {
        cat(sprintf("  [COVERAGE FAIL] %s: %.4f < %.4f\n", basename(f), coverage, coverage_threshold))
        results <- rbind(results, data.frame(
          Folder=folder_name, Genome=genomename, Compare="phage vs host_nonphage",
          Phage_File=basename(f), Coverage=round(coverage,4),
          Test_Type=NA_character_, p_value=NA_real_, Significant="Coverage failed",
          BM_Prob_Superiority=NA_real_, Q95=NA_real_, Q99=NA_real_, Phage_Sites=phage_sites,
          Phage_OverQ95=NA_real_, Phage_OverQ99=NA_real_, RR_Q95=NA_real_, RR_Q99=NA_real_,
          Delta_Mean=NA_real_, Delta_CI_2.5=NA_real_, Delta_CI_97.5=NA_real_))
        next
      }

      test_type<-NA_character_; pval<-NA_real_; significant<-"No"
      bm_prob<-NA_real_; phage_over_q95<-NA_real_; phage_over_q99<-NA_real_
      rr_95<-NA_real_; rr_99<-NA_real_
      delta_mean<-NA_real_; delta_CI<-c(NA_real_,NA_real_)

      jb1  <- tryCatch(jarque.bera.test(phage_depth), error=function(e) NULL)
      jb2  <- tryCatch(jarque.bera.test(host_depth),  error=function(e) NULL)
      norm1 <- !is.null(jb1)&&!is.na(jb1$p.value)&&jb1$p.value>0.05
      norm2 <- !is.null(jb2)&&!is.na(jb2$p.value)&&jb2$p.value>0.05

      if (norm1 & norm2) {
        vtest <- var.test(phage_depth, host_depth)
        if (vtest$p.value>0.05) {
          ttest <- t.test(phage_depth,host_depth,alternative="greater",var.equal=TRUE)
          test_type <- "Student's t-test"
        } else {
          ttest <- t.test(phage_depth,host_depth,alternative="greater",var.equal=FALSE)
          test_type <- "Welch's t-test"
        }
        pval <- ttest$p.value
        if (!is.na(pval)&&pval<alpha) significant <- "Yes"
      } else {
        bmtest    <- brunnermunzel::brunnermunzel.test(phage_depth,host_depth,alternative="greater")
        test_type <- "Brunner-Munzel"; pval <- bmtest$p.value
        bm_prob   <- suppressWarnings(as.numeric(bmtest$estimate))
        if (!is.na(pval)&&pval<alpha) {
          significant    <- "Yes"
          phage_over_q95 <- sum(phage_depth>q95); phage_over_q99 <- sum(phage_depth>q99)
          rr_95 <- ifelse(p0_95==0, NA_real_, (phage_over_q95/phage_sites)/p0_95)
          rr_99 <- ifelse(p0_99==0, NA_real_, (phage_over_q99/phage_sites)/p0_99)
          # Permutation-based Cliff's delta (parallel); iteration count from --iter
          set.seed(123)
          len_h <- min(10000L, length(host_depth))
          cliff_deltas <- foreach(i=1:n_iter, .combine=c, .packages="effsize") %dopar% {
            sub_h <- host_depth[sample.int(length(host_depth), len_h, replace=FALSE)]
            res   <- try(effsize::cliff.delta(phage_depth, sub_h), silent=TRUE)
            if (inherits(res,"try-error")) NA_real_ else suppressWarnings(as.numeric(res$estimate))
          }
          cliff_deltas <- cliff_deltas[is.finite(cliff_deltas)]
          if (length(cliff_deltas)>0) {
            delta_mean <- mean(cliff_deltas, na.rm=TRUE)
            qs <- stats::quantile(cliff_deltas, c(0.025,0.975), na.rm=TRUE, names=FALSE, type=7)
            delta_CI <- c(qs[1], qs[2])
          }
        }
      }
      results <- rbind(results, data.frame(
        Folder=folder_name, Genome=genomename, Compare="phage vs host_nonphage",
        Phage_File=basename(f), Coverage=round(coverage,4),
        Test_Type=test_type, p_value=pval, Significant=significant,
        BM_Prob_Superiority=bm_prob, Q95=q95, Q99=q99, Phage_Sites=phage_sites,
        Phage_OverQ95=phage_over_q95, Phage_OverQ99=phage_over_q99,
        RR_Q95=rr_95, RR_Q99=rr_99,
        Delta_Mean=delta_mean, Delta_CI_2.5=delta_CI[1], Delta_CI_97.5=delta_CI[2]))
    }
  }
  out1 <- file.path(sample_dir, sprintf("%s_phage_depth_stats_summary_final.csv", folder_name))
  out2 <- file.path(sample_dir, sprintf("%s_basic_stats_summary.csv", folder_name))
  write.csv(results,      out1, row.names=FALSE)
  write.csv(stat_results, out2, row.names=FALSE)
  cat(sprintf("  Output: %s | %s\n", out1, out2))
}
cat("\nAll samples processed.\n")
REOF_PRECISE

fi   # end mode branch

# Collect final result CSVs into activity_results directory
final_res_dir="$output/activity_results/$mode"
mkdir -p "$final_res_dir"
find "$output_map_dir" -name "*_phage_depth_stats_summary_final.csv" -print0 \
| while IFS= read -r -d '' file; do
    mv "$file" "$final_res_dir/"
done
echo "    Analysis complete. Results saved to: $final_res_dir"

# ==========================================
# 6. Optional cleanup after R analysis
# ==========================================

# Delete depth files (--rmdepth or --rmall)
if [ "$rm_depth" = true ]; then
    echo ">>> [--rmdepth] Removing per-region depth .txt files..."
    find "$output_map_dir" -name "phage_depth_*.txt" -delete
    find "$output_map_dir" -name "*_mgs_depth.txt" -delete
    find "$output_map_dir" -name "*_host_nonphage_depth.txt" -delete
    echo "    Depth files removed."
fi

# Delete Prodigal output (--rmprodigal or --rmall, fast mode only)
if [ "$rm_prodigal" = true ] && [ "$mode" == "fast" ]; then
    echo ">>> [--rmprodigal] Removing Prodigal output directory: $output_host_faa"
    rm -rf "$output_host_faa"
    echo "    Prodigal output removed."
fi

# Delete FetchMGs-related files (--rmmgs or --rmall, fast mode only)
if [ "$rm_mgs" = true ] && [ "$mode" == "fast" ]; then
    echo ">>> [--rmmgs] Removing FetchMGs-related files..."
    rm -rf "$fetchmgs_dir"
    rm -f  "$output/target_cogs_list.txt"
    rm -f  "$output/mgs_missing_genomes.log"
    rm -f  "$output/mgs_ordinate.tsv"
    rm -rf "$output/mgs_ext_results"
    echo "    FetchMGs-related files removed."
fi

# Delete entire mapping_results tree (--rmall only)
# BAMs were already removed per-sample during depth extraction if rm_bam=true.
# Here we remove the whole directory (depth subdirs, BED files, etc.).
if [ "$rm_all" = true ]; then
    echo ">>> [--rmall] Removing mapping_results directory: $output_map_dir"
    rm -rf "$output_map_dir"
    echo "    mapping_results removed."
    echo ">>> [--rmall] Cleanup complete. Only activity_results/ is retained."
fi

echo ">>> Pipeline complete."
