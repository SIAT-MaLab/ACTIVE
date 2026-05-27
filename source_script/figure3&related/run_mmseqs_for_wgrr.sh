#!/bin/bash
SENSITIVITY=7.5
EVALUE=1e-4
COVERAGE=0.5
PIDENT=0.35
THREADS=16

usage() {
    echo "Usage: $0 -i <input.faa> -o <output_prefix> [options]"
    echo "Options:"
    echo "  -i  Input protein FASTA file (required)"
    echo "  -o  Output file prefix (required, generates <prefix>_all_vs_all.tsv)"
    echo "  -t  Number of threads (default: $THREADS)"
    echo "  -s  MMseqs2 sensitivity (default: $SENSITIVITY)"
    echo "  -e  E-value threshold (default: $EVALUE)"
    echo "  -c  Coverage threshold (default: $COVERAGE)"
    echo "  -p  Sequence identity threshold (default: $PIDENT)"
    echo "  -h  Display this help message"
    exit 1
}

while getopts "i:o:t:s:e:c:p:h" opt; do
    case "$opt" in
        i) INPUT_FAA=$OPTARG ;;
        o) PREFIX=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        s) SENSITIVITY=$OPTARG ;;
        e) EVALUE=$OPTARG ;;
        c) COVERAGE=$OPTARG ;;
        p) PIDENT=$OPTARG ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "$INPUT_FAA" ] || [ -z "$PREFIX" ]; then
    echo "Error: Missing required parameters -i or -o"
    usage
fi

if [ ! -f "$INPUT_FAA" ]; then
    echo "Error: Cannot find input file $INPUT_FAA"
    exit 1
fi

echo "=================================================="
echo "Starting MMseqs2 All-vs-All alignment"
echo "Input file: $INPUT_FAA"
echo "Output prefix: $PREFIX"
echo "Parameters: -s $SENSITIVITY | -e $EVALUE | -c $COVERAGE | --min-seq-id $PIDENT | Threads $THREADS"
echo "=================================================="

DB_NAME="${PREFIX}_DB"
RES_NAME="${PREFIX}_ResultDB"
TMP_DIR="${PREFIX}_tmp"
OUT_TSV="${PREFIX}_all_vs_all.tsv"

mkdir -p "$TMP_DIR"

echo "1/3 Creating MMseqs2 database..."
mmseqs createdb "$INPUT_FAA" "$DB_NAME" > /dev/null

echo "2/3 Running all-vs-all alignment (this may take some time)..."
mmseqs search "$DB_NAME" "$DB_NAME" "$RES_NAME" "$TMP_DIR" -a \
    -s "$SENSITIVITY" \
    -e "$EVALUE" \
    -c "$COVERAGE" \
    --min-seq-id "$PIDENT" \
    --threads "$THREADS" > /dev/null

echo "3/3 Generating BLAST m8 format table..."
mmseqs convertalis "$DB_NAME" "$DB_NAME" "$RES_NAME" "$OUT_TSV" \
    --format-output "query,target,pident,evalue,qcov,tcov,qlen,tlen,alnlen" > /dev/null

echo "Alignment complete! Results saved to $OUT_TSV"

echo "Cleaning up temporary files..."
rm -rf "$TMP_DIR" "${DB_NAME}"* "${RES_NAME}"*

echo "Process finished successfully!"