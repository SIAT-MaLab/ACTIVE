import pandas as pd
import re
INPUT_FILE = "vfc_full_flow_breakdown_glmm.tsv"
OUTPUT_SANKEY = "sankey_links_with_pct.tsv"

def main():
    df = pd.read_csv(INPUT_FILE, sep='\t')

    links_dict = {}

    for _, row in df.iterrows():
        vfc_source = row['VFC_Source']
        if pd.isna(vfc_source):
            continue
        vfc_node = f"VFC_{int(float(vfc_source))}" 
        target_node = f"New_{row['Target_vFAM']}"
        vfc_count = int(row['Flow_Count'])
        links_dict[(vfc_node, target_node)] = links_dict.get((vfc_node, target_node), 0) + vfc_count
        old_fams_str = str(row['User_Original_Families'])
        if old_fams_str != "None" and old_fams_str != "nan":
            items = old_fams_str.split(';')
            for item in items:
                match = re.match(r'\s*(.+?)\((\d+)\)', item)
                if match:
                    old_fam_name = match.group(1).strip()
                    user_count = int(match.group(2))
                    old_fam_node = f"Old_{old_fam_name}"
                    links_dict[(target_node, old_fam_node)] = user_count

    sankey_data = []
    for (src, tgt), val in links_dict.items():
        sankey_data.append({
            'source': src,
            'target': tgt,
            'value': val
        })

    df_sankey = pd.DataFrame(sankey_data)

    source_sums = df_sankey.groupby('source')['value'].sum().reset_index(name='sumofsource')
    df_sankey = pd.merge(df_sankey, source_sums, on='source', how='left')
    df_sankey['%'] = df_sankey['value'] / df_sankey['sumofsource']
    df_sankey = df_sankey[['source', 'target', 'value', 'sumofsource', '%']]
    df_sankey = df_sankey.sort_values(by=['source', '%'], ascending=[True, False])
    df_sankey.to_csv(OUTPUT_SANKEY, sep='\t', index=False)

if __name__ == "__main__":
    main()
