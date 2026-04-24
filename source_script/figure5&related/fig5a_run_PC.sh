#All proteins were clustered at 30% AAI and 70% alignment coverage using MMseqs2 ref:Metagenomic compendium of 189,680 DNA viruses from the human gut microbiome
mmseqs easy-cluster "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_checkv_quality/hq_proteins.faa" mmseqs_out tmp_dir --min-seq-id 0.3 -c 0.7 --cov-mode 0 --threads 16
cut -f1 mmseqs_out_cluster.tsv | sort | uniq | \
awk '{printf("%s\tPC_%05d\n", $0, NR)}' > temp_rep_to_id.tsv
awk 'NR==FNR{map[$1]=$2; next} {print $2"\t"map[$1]"\t"$1}' \
temp_rep_to_id.tsv mmseqs_out_cluster.tsv > allpp_PC_map.tsv
sed -i '1i protein_id\tpc_id\trep_seq_id' allpp_PC_map.tsv

metacerberus.py --protein mmseqs_out_rep_seq.faa --hmm COG, PFAM, KOFam_prokaryote, PHROG, VOG --dir_out mmseqs_repseq_annot --cpus 24
