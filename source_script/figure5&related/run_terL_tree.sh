NR==FNR {
    if(NR>1) {
        clean_id = $5;
        sub(/_CDS_0+/, "_", clean_id);
        map[clean_id] = $1; 
    }
    next
}

{
    if ($0 ~ /^>/) {
        full_id = $1;
        sub(/^>/, "", full_id);
        if (full_id in map) {
            sub(/^>[^ ]+/, ">"map[full_id], $0);
            print $0;
        } else {
            print $0;
        }
    } else {
        print $0;
    }
}' TerL_All_vOTU_Table.tsv vFAM1_vOTU_terL.faa > tmp_terL.faa
mafft --thread 20 --auto tmp_terL.faa > terL_aligned.fasta
echo "!!!!!!!!!!mafft done !!!!!!!!!"
fasttree -wag terL_aligned.fasta > terL_fasttree.tree
echo "!!!!!!!!!!fasttree done !!!!!!!!!"
mv terL_fasttree.tree vFAM1_votu_terL.tree
