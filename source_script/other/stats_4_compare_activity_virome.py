import csv

def analyze_virome_subsample():
    file_results_csv = "all_final_sensitive_results.csv"  
    file_virome_txt = "all_virome_coverage.txt"           
    output_detail_file = "phage_category_details.tsv"     
    ACTIVITY_THRESHOLD = 0.7        
    VIROME_COVERAGE_THRESHOLD = 70 
    activity_map = {}
    
    try:
        with open(file_results_csv, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if len(row) < 5: continue

                raw_name = row[3].strip()
                clean_id = raw_name.replace("phage_depth_", "").replace(".txt", "")

                try:
                    score_str = row[-3].strip() 
                except IndexError:
                    score_str = ""

                try:
                    if score_str.upper() == "NA" or score_str == "":
                        score = 0.0
                    else:
                        score = float(score_str)
                except ValueError:
                    score = 0.0
                
                activity_map[clean_id] = score
                
    except FileNotFoundError:
        print(f"missing {file_results_csv}")
        return

    stats = {
        "Active_Detected": 0,       
        "Active_NotDetected": 0,    
        "Inactive_Detected": 0,     
        "Inactive_NotDetected": 0,  
        "Missing_Data": 0
    }
    
    try:
        with open(file_virome_txt, 'r') as f_in, open(output_detail_file, 'w') as f_out:

            header_out = ["Phage_ID", "Activity_Score", "Predicted_Status", 
                          "Virome_Coverage", "Virome_Status", "Category"]
            f_out.write("\t".join(header_out) + "\n")
            header_in = f_in.readline() 
            
            for line in f_in:
                if not line.strip(): continue
                parts = line.split()
                if len(parts) < 5: continue

                p_id = parts[0].strip()
                try:
                    coverage = float(parts[4]) 
                except ValueError:
                    continue

                score = activity_map.get(p_id, None)
                
                if score is None:
                    stats["Missing_Data"] += 1
                    f_out.write(f"{p_id}\tNA\tUnknown\t{coverage}\tUnknown\tMissing_Score_Data\n")
                    continue

                is_active = score >= ACTIVITY_THRESHOLD
                is_detected = coverage >= VIROME_COVERAGE_THRESHOLD
                
                predicted_str = "Active" if is_active else "Inactive"
                virome_str = "Detected" if is_detected else "Not_Detected"

                category = ""
                if is_active and is_detected:
                    category = "Active_Detected"
                    stats["Active_Detected"] += 1
                elif is_active and not is_detected:
                    category = "Active_NotDetected"
                    stats["Active_NotDetected"] += 1
                elif not is_active and is_detected:
                    category = "Inactive_Detected"
                    stats["Inactive_Detected"] += 1
                elif not is_active and not is_detected:
                    category = "Inactive_NotDetected"
                    stats["Inactive_NotDetected"] += 1

                out_line = f"{p_id}\t{score}\t{predicted_str}\t{coverage}\t{virome_str}\t{category}\n"
                f_out.write(out_line)

    except FileNotFoundError:
        print(f"missing {file_virome_txt}")
        return


    print("\n" + "="*60)
    print("           Classification Summary")
    print("="*60)
    
    total_processed = sum(stats.values()) - stats["Missing_Data"]
    
    print(f"total_processed: {total_processed}")
    print("-" * 60)
    print(f"1. Active & Detected: \t{stats['Active_Detected']}")
    print(f"2.(Inactive & Not Det): \t{stats['Inactive_NotDetected']}")
    print("-" * 60)
    print(f"3.  (Active but Not Det): \t{stats['Active_NotDetected']}")
    print(f"4. (Inactive but Det): \t{stats['Inactive_Detected']}")
    print("-" * 60)


if __name__ == "__main__":
    analyze_virome_subsample()
