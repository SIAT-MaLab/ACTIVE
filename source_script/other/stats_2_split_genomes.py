import sys
import os

def split_fasta(list_file, fasta_file, out_dir="split_genomes"):
    target_genomes = set()
    with open(list_file, 'r') as f:
        for line in f:
            if line.strip():
                target_genomes.add(line.strip())
    
    print(f" {len(target_genomes)}")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)   
    found_genomes = set()
    current_file_handle = None
    current_genome_id = None

    line_count = 0
    
    with open(fasta_file, 'r') as f_in:
        for line in f_in:
            line_count += 1
            if line_count % 100000 == 0:
                print(f" {line_count}", end='\r')

            if line.startswith('>'):
                header = line.strip().lstrip('>')  

                if '_' in header:
                    genome_id = header.rsplit('_', 1)[0]
                else:
                    genome_id = header 
 
                if genome_id in target_genomes:
                    found_genomes.add(genome_id)
                    if genome_id != current_genome_id:
                        if current_file_handle:
                            current_file_handle.close()
                        out_path = os.path.join(out_dir, f"{genome_id}.fasta")
                        current_file_handle = open(out_path, 'a')
                        current_genome_id = genome_id
                    current_file_handle.write(line)
                else:
                    current_genome_id = None
                    if current_file_handle:
                        current_file_handle.close()
                        current_file_handle = None
            else:
                if current_file_handle:
                    current_file_handle.write(line)

    if current_file_handle:
        current_file_handle.close()
    
    missing_genomes = target_genomes - found_genomes
    
    print("-" * 30)
    print(f"target: {len(target_genomes)}")
    print(f"actual {len(found_genomes)}")
    
    if len(missing_genomes) == 0:
        print("✅ all genome in ref.fasta 中。")
    else:
        print(f"⚠️ caution：there are {len(missing_genomes)} missing in ref.fasta")
        with open('missing_genomes.log', 'w') as f_log:
            for miss in missing_genomes:
                f_log.write(miss + '\n')

if __name__ == "__main__":
    genome_list = "genome.list"
    ref_fasta = "ref.fasta"
    split_fasta(genome_list, ref_fasta)

