import pandas as pd
from Bio import SeqIO 
import os
import sys
MAP_FILE = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/genome_to_vOTU.tsv"

ANNOT_FILE = "allpharokka_annot.tsv"

FASTA_FILE = "allterL.faa" 

OUTPUT_FASTA = "vOTU_terL_for_tree.faa"
OUTPUT_TSV = "vOTU_selection_summary.tsv"

def main():
    df_map = pd.read_csv(MAP_FILE, sep='\t', header=None, names=['Genome', 'vOTU', 'Rep'])
    votu_dict = {}
    for _, row in df_map.iterrows():
        votu = row['vOTU']
        genome = row['Genome']
        rep = row['Rep']
        
        if votu not in votu_dict:
            votu_dict[votu] = {'rep': rep, 'members': []}

        if genome != rep:
            votu_dict[votu]['members'].append(genome)

    df_annot = pd.read_csv(ANNOT_FILE, sep='\t')

    term_pattern = "terminase large"
    mask = df_annot['annot'].str.contains(term_pattern, case=False, na=False)
    df_terl = df_annot[mask][['gene', 'contig', 'annot']]

    genome_terl_map = {}
    valid_genes = set()     
    for _, row in df_terl.iterrows():
        genome = row['contig']
        gene_id = row['gene']
        
        if genome not in genome_terl_map:
            genome_terl_map[genome] = []
        genome_terl_map[genome].append(gene_id)
        valid_genes.add(gene_id)
       
    seq_dict = {}
    try:
        for record in SeqIO.parse(FASTA_FILE, "fasta"):
            if record.id in valid_genes:
                seq_dict[record.id] = str(record.seq)

    summary_records = []
    fasta_records = []
    
    for votu, info in votu_dict.items():
        rep_genome = info['rep']
        members = info['members']
        
        selected_gene = None
        selected_genome = None
        selection_type = None # Representative or Member
        note = ""
        rep_terls = genome_terl_map.get(rep_genome, [])
        
        if len(rep_terls) == 1:
            selected_gene = rep_terls[0]
            selected_genome = rep_genome
            selection_type = "Representative"
            note = "Pass: Representative has unique TerL"
        else:
            if len(rep_terls) == 0:
                failure_reason = "Rep has 0 TerL"
            else:
                failure_reason = f"Rep has {len(rep_terls)} TerLs"
            
            found_rescue = False
            for member in members:
                mem_terls = genome_terl_map.get(member, [])
                if len(mem_terls) == 1:
                    selected_gene = mem_terls[0]
                    selected_genome = member
                    selection_type = "Member_Rescue"
                    note = f"Pass: {failure_reason}, found valid member"
                    found_rescue = True
                    break 
            
            if not found_rescue:
                selection_type = "Failed"
                note = f"Fail: {failure_reason}, and no member has unique TerL"

        if selected_gene:
            if selected_gene in seq_dict:
                header = f"{votu} original_id={selected_gene} source={selection_type}"
                sequence = seq_dict[selected_gene]
                fasta_records.append(f">{header}\n{sequence}\n")
                
                summary_records.append({
                    'vOTU': votu,
                    'Status': 'Selected',
                    'Selection_Type': selection_type,
                    'Selected_Genome': selected_genome,
                    'Selected_Gene': selected_gene,
                    'Note': note
                })
            else:
                summary_records.append({
                    'vOTU': votu,
                    'Status': 'Error',
                    'Selection_Type': selection_type,
                    'Selected_Genome': selected_genome,
                    'Selected_Gene': selected_gene,
                    'Note': "Gene found in annot but missing in FASTA file"
                })
        else:
            summary_records.append({
                'vOTU': votu,
                'Status': 'Not_Found',
                'Selection_Type': 'None',
                'Selected_Genome': 'None',
                'Selected_Gene': 'None',
                'Note': note
            })

    
    df_summary = pd.DataFrame(summary_records)
    cols = ['vOTU', 'Status', 'Selection_Type', 'Selected_Genome', 'Selected_Gene', 'Note']
    df_summary = df_summary[cols]
    df_summary.to_csv(OUTPUT_TSV, sep='\t', index=False)
    with open(OUTPUT_FASTA, 'w') as f:
        f.writelines(fasta_records)

if __name__ == "__main__":
    main()
