import csv
import sys

def calculate_stats():
    file_genome_list = "genome.list"
    file_missing_log = "all_mgs_missing_genomes.log"
    file_results_csv = "all_final_sensitive_results.csv"
    output_active_phages = "active_phages_list.txt"
    
    ACTIVITY_THRESHOLD = 0.7
    expected_genomes = set()
    try:
        with open(file_genome_list, 'r') as f:
            for line in f:
                if line.strip():
                    expected_genomes.add(line.strip())
    except FileNotFoundError:
        print(f"missing {file_genome_list}")
        return

    print(f"1. E: {len(expected_genomes)}")
    missing_genomes = set()
    try:
        with open(file_missing_log, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    missing_genomes.add(parts[1])
    except FileNotFoundError:
        print(f"misiing {file_missing_log}")

    print(f"2. No MGS: {len(missing_genomes)}")

    valid_genomes = expected_genomes - missing_genomes
    n_total_calculated = len(valid_genomes)
    
    print(f"3. total_calculated_genome : {n_total_calculated}")
    if n_total_calculated == 0:
        print("no genome")
        return

    lysogen_genomes = set()          
    active_lysogen_genomes = set()       
    total_phages_count = 0           
    active_phages_count = 0             
    active_phage_names = []         
    try:
        with open(file_results_csv, 'r') as csvfile:
            reader = csv.reader(csvfile)
            
            for row in reader:
                if len(row) < 5: 
                    continue
                genome_id = row[1].strip()
                phage_id = row[3].strip()
                score_str = row[-3].strip()
                if genome_id in valid_genomes:
                    lysogen_genomes.add(genome_id)
                    total_phages_count += 1
                    try:
                        if score_str.upper() == "NA" or score_str == "":
                            score = 0.0
                        else:
                            score = float(score_str)
                    except ValueError:
                        score = 0.0 
                    if score >= ACTIVITY_THRESHOLD:
                        active_phages_count += 1
                        active_lysogen_genomes.add(genome_id)
                        active_phage_names.append(phage_id)
                        
    except FileNotFoundError:
        print(f"missing {file_results_csv}")
        return
    n_lysogens = len(lysogen_genomes)
    n_active_lysogens = len(active_lysogen_genomes)
    lysogeny_rate = n_lysogens / n_total_calculated if n_total_calculated > 0 else 0
    genome_induction_rate = n_active_lysogens / n_lysogens if n_lysogens > 0 else 0
    phage_induction_rate = active_phages_count / total_phages_count if total_phages_count > 0 else 0


    print(f"A. n_total_calculated: \t{n_total_calculated}")
    print(f"B. n_lysogens: \t{n_lysogens}")
    print(f"C. n_active_lysogens: \t{n_active_lysogens}")
    print("-" * 40)
    print(f"D. total_phages_count: \t\t{total_phages_count}")
    print(f"E. active_phages_count: \t\t{active_phages_count}")
    print("=" * 40)
    
    print(f"1. lysogeny_rate (B / A): \t\t{lysogeny_rate:.2%}")
    print(f"2. genome_induction_rate (C / B): \t\t{genome_induction_rate:.2%}")
    print(f"3. phage_induction_rate (E / D): \t\t{phage_induction_rate:.2%}")
    print("=" * 40)

    with open(output_active_phages, 'w') as f_out:
        for name in active_phage_names:
            f_out.write(name + "\n")


if __name__ == "__main__":
    calculate_stats()
