import pandas as pd
metadata_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_checkv_quality/hqphage_metadata.tsv"
fasta_file = "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_annot/allpp_proteins.faa"
output_fasta = "hq_proteins.faa"
df = pd.read_csv(metadata_file, sep='\t')

hq_phages = set(df['Phage_Name'].dropna().astype(str))
extracted_count = 0
write_flag = False

with open(fasta_file, 'r') as fin, open(output_fasta, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
           header_id = line[1:].strip().split()[0]
              if '_CDS_' in header_id:
                phage_id = header_id.rsplit('_CDS_', 1)[0]
            else:
                phage_id = header_id 
            if phage_id in hq_phages:
                write_flag = True
                fout.write(line)
                extracted_count += 1
            else:
                write_flag = False
        else:

            if write_flag:
                fout.write(line)

