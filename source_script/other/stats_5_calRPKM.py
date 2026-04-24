import csv

def calculate_virome_rpkm():
    file_virome_txt = "all_virome_coverage.txt" 
    file_activity_csv = "all_final_sensitive_results.csv"
    output_file = "virome_abundance_RPKM.tsv"

    ACTIVITY_THRESHOLD = 0.7 

    COL_ID = 0
    COL_AVG_FOLD = 1
    COL_LENGTH = 2
    COL_COV_PCT = 4
    COL_PLUS_READS = 6  
    COL_MINUS_READS = 7 

    activity_map = {}
    try:
        with open(file_activity_csv, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) < 5: continue

                raw_name = row[3].strip()
                clean_id = raw_name.replace("phage_depth_", "").replace(".txt", "")

                try:
                    score_str = row[-3].strip()
                    if score_str.upper() == "NA" or score_str == "":
                        score = 0.0
                    else:
                        score = float(score_str)
                except:
                    score = 0.0
                activity_map[clean_id] = score
    except FileNotFoundError:
        print(f"missing {file_activity_csv}")
        return


    total_mapped_reads = 0
    line_count = 0
    
    try:
        with open(file_virome_txt, 'r') as f:
            for line in f:
                if line.startswith("#"): continue 
                if not line.strip(): continue
                
                parts = line.split()
                if len(parts) <= COL_MINUS_READS: continue
                
                try:
                    p_reads = int(parts[COL_PLUS_READS])
                    m_reads = int(parts[COL_MINUS_READS])
                    total_mapped_reads += (p_reads + m_reads)
                    line_count += 1
                except ValueError:
                    continue 
                    
    except FileNotFoundError:
        print(f"missing {file_virome_txt}")
        return

    count_active = 0
    count_inactive = 0
    
    with open(file_virome_txt, 'r') as f_in, open(output_file, 'w') as f_out:
        f_out.write("Phage_ID\tGroup\tActivity_Score\tTotal_Reads\tLength\tRPKM\tAvg_Fold\tCovered_Percent\n")
        
        for line in f_in:
            if line.startswith("#"): continue
            if not line.strip(): continue
            parts = line.split()
            if len(parts) <= COL_MINUS_READS: continue
            
            p_id = parts[COL_ID].strip()
            score = activity_map.get(p_id, None)
              if score is None: continue

            try:
                length = float(parts[COL_LENGTH])
                avg_fold = float(parts[COL_AVG_FOLD])
                covered_pct = float(parts[COL_COV_PCT])
                reads_count = int(parts[COL_PLUS_READS]) + int(parts[COL_MINUS_READS])
            except ValueError:
                continue

            group = "Active" if score >= ACTIVITY_THRESHOLD else "Inactive"
            if group == "Active": count_active += 1
            else: count_inactive += 1

            if length > 0:
                rpkm = (reads_count * 1_000_000_000) / (length * total_mapped_reads)
            else:
                rpkm = 0
            
            f_out.write(f"{p_id}\t{group}\t{score}\t{reads_count}\t{length}\t{rpkm:.4f}\t{avg_fold}\t{covered_pct}\n")


if __name__ == "__main__":
    calculate_virome_rpkm()
