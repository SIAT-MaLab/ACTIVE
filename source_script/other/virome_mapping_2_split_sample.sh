#!/bin/bash

# ===========================================
FASTA_FILE="virome_pp.fna"
TEMP_SPLIT_DIR="temp_split_seqs"
# ===========================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

if ! command -v seqkit &> /dev/null; then
    echo -e "${RED}Error: seqkit not found, please make sure it is installed.${NC}"
    exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then
    echo -e "${RED}Error: File $FASTA_FILE does not exist.${NC}"
    exit 1
fi

echo -e "${BLUE}=== Starting Dry Run ===${NC}"
echo -e "Analyzing sample IDs in $FASTA_FILE..."

grep "^>" "$FASTA_FILE" | awk -F'-' '{print substr($1,2)}' | sort | uniq > existing_fasta_ids.txt

declare -A fasta_id_map
while read id; do fasta_id_map["$id"]=1; done < existing_fasta_ids.txt
rm existing_fasta_ids.txt

ls *_virome_*.fastq.gz 2>/dev/null | sed 's/_virome.*//' | sort | uniq > sample_list.txt

if [ ! -s sample_list.txt ]; then
    echo -e "${RED}No sequencing files matching *_virome_*.fastq.gz found.${NC}"
    exit 1
fi

echo -e "${YELLOW}----------------------------------------${NC}"
echo -e "${YELLOW}            DRY RUN PLAN REPORT          ${NC}"
echo -e "${YELLOW}----------------------------------------${NC}"

declare -a plan_details

while read sample; do
    fq1="${sample}_virome_1.fastq.gz"
    fq2="${sample}_virome_2.fastq.gz"
    
    msg="${GREEN}[Create Dir]${NC} $sample/"
    
    if [[ -f "$fq1" && -f "$fq2" ]]; then
        msg="$msg\n    |-- ${BLUE}[Move File]${NC} $fq1, $fq2 -> $sample/"
    elif [[ -f "$fq1" ]]; then
         msg="$msg\n    |-- ${BLUE}[Move File]${NC} $fq1 -> $sample/ (Missing read2)"
    else
         msg="$msg\n    |-- ${RED}[Missing File]${NC} Sequencing files not found"
    fi

    if [ "${fasta_id_map[$sample]}" == "1" ]; then
        msg="$msg\n    |-- ${BLUE}[Assign Seq]${NC} Split $sample.fna from $FASTA_FILE -> $sample/"
    else
        msg="$msg\n    |-- ${YELLOW}[No Seq]${NC} No record for sample $sample in $FASTA_FILE"
    fi
    
    echo -e "$msg"
    
done < sample_list.txt

echo -e "${YELLOW}----------------------------------------${NC}"
echo -e "Note: After execution, fastq.gz files in the current directory will be moved to subdirectories."

read -p "Confirm to execute the above operations? [y/N] " response
if [[ ! "$response" =~ ^[yY]$ ]]; then
    echo "Operation cancelled."
    rm sample_list.txt
    exit 0
fi

echo -e "\n${BLUE}=== Execution Started ===${NC}"

echo -e "Splitting $FASTA_FILE (This may take a while)..."
mkdir -p "$TEMP_SPLIT_DIR"

seqkit split "$FASTA_FILE" --by-id --id-regexp "^([^-]+)-" -O "$TEMP_SPLIT_DIR" --force

if [ $? -ne 0 ]; then
    echo -e "${RED}SeqKit split failed, script terminated.${NC}"
    exit 1
fi

echo -e "Moving and organizing files..."

while read sample; do
    mkdir -p "$sample"
    
    mv "${sample}_virome_"*.fastq.gz "$sample/" 2>/dev/null
    
    split_file=$(find "$TEMP_SPLIT_DIR" -name "*id_${sample}.*" | head -n 1)
    
    if [ -z "$split_file" ]; then
        split_file=$(find "$TEMP_SPLIT_DIR" -name "*${sample}.*" | head -n 1)
    fi
    
    if [ -n "$split_file" ]; then
        mv "$split_file" "$sample/${sample}_virome.fna"
    fi
    
done < sample_list.txt

echo -e "Cleaning up temporary files..."
rm sample_list.txt

count_remaining=$(ls "$TEMP_SPLIT_DIR" 2>/dev/null | wc -l)
if [ "$count_remaining" -gt 0 ]; then
    echo -e "${YELLOW}Notice: $count_remaining files remain in $TEMP_SPLIT_DIR (possibly orphan sequences without matching fastq).${NC}"
else
    rmdir "$TEMP_SPLIT_DIR" 2>/dev/null
fi

echo -e "${GREEN}All operations completed successfully!${NC}"