import pandas as pd

def analyze_sample_metadata(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
    except FileNotFoundError:
        print(f"error {file_path}")
        return

    expected_count = df['External ID'].count()
    actual_meta_count = df['calc_metagenomic'].count()
    actual_mag_count = df['used_bins'].count()
    print(f"1.  (External ID): {expected_count}")
    print(f"2.  (calc_metagenomic): {actual_meta_count}")
    print(f"3.  (used_bins): {actual_mag_count}")
    print("-" * 30)


    virome_set = set(df['calced_virome'].dropna().astype(str).str.strip())
    meta_set = set(df['calc_metagenomic'].dropna().astype(str).str.strip())
    mag_set = set(df['used_bins'].dropna().astype(str).str.strip())
    print(f"used_virome: {len(virome_set)}")
    
    missing_in_meta = virome_set - meta_set
    if len(missing_in_meta) == 0:
        print("✅ all virme included in metagenomic。")
    else:
        print(f"❌ there are {len(missing_in_meta)} virome missing")



    missing_in_mag = virome_set - mag_set
    if len(missing_in_mag) == 0:
        print("✅all virme included in metagenomic")
    else:
        print(f"❌ there are {len(missing_in_meta)} virome missing, indicating no high-quality MAG")


if __name__ == "__main__":

    file_name = 'sample_used_metadata.tsv'
    analyze_sample_metadata(file_name)
