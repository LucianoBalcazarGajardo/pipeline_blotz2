#!/usr/bin/env python3
import os
import re
import pandas as pd
import subprocess
import argparse
import matplotlib.pyplot as plt
from Bio import SeqIO

# ==========================================
# MASTER FASTA FILE CONFIGURATION
# Remember to edit this path for your file
# ==========================================
MASTER_FASTA_PATH = "specific file path for your RNA.fasta"

def generate_iptm_weblogo(folder, protein, percentage):
    print(f"=== STARTING WEBLOGO & HISTOGRAM ANALYSIS BY ipTM (Top {percentage}%) ===")
    
    # 1. Read the general index (all models before filtering)
    csv_path = os.path.join(folder, f"01_indice_general_modelos_{protein}.csv")
    
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found. Run the main blitz.py pipeline first.")
        return

    # 2. Sort by ipTM and isolate the top %
    df = pd.read_csv(csv_path)
    df_sorted = df.sort_values(by='ipTM', ascending=False).reset_index(drop=True)
    top_count = max(1, int(len(df_sorted) * (percentage / 100.0)))
    df_top = df_sorted.head(top_count)
    
    print(f"1. Selected {top_count} models based on ipTM score.")
    
    # Clean file naming (removed "_best_models")
    base_name = f"top{percentage}pct_ipTM_{protein}"
    csv_out_path = os.path.join(folder, f"{base_name}.csv")
    fasta_out_path = os.path.join(folder, f"{base_name}_sequences.fasta")
    logo_path = os.path.join(folder, f"{base_name}_weblogo.png")
    hist_path = os.path.join(folder, f"{base_name}_histogram.png")
    
    df_top.to_csv(csv_out_path, index=False)

    # 3. Generate ipTM Histogram for this Top subset
    print("2. Generating ipTM histogram")
    plt.figure(figsize=(8, 5))
    plt.hist(df_top['ipTM'], bins=15, color='#2CA02C', edgecolor='black')
    plt.title(f'ipTM Distribution (Top {percentage}% - {protein})', weight='bold')
    plt.xlabel('ipTM Score')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(hist_path, dpi=300)
    plt.close()

    # 4. Load master FASTA and cross-reference data
    print("3. Cross-referencing RNA sequences")
    rna_dict = {}
    try:
        for record in SeqIO.parse(MASTER_FASTA_PATH, "fasta"):
            rna_dict[record.id] = str(record.seq)
    except FileNotFoundError:
        print(f"Error: Master FASTA file not found at {MASTER_FASTA_PATH}")
        return
        
    matches_found = 0
    with open(fasta_out_path, 'w') as f_out:
        for idx, row in df_top.iterrows():
            # Looks for "ARN" followed by numbers
            match = re.search(r'(ARN\d+)', str(row['Carpeta']) + str(row['Nombre_CIF']))
            if match:
                arn_id = match.group(1)
                if arn_id in rna_dict:
                    f_out.write(f">{arn_id}\n{rna_dict[arn_id]}\n")
                    matches_found += 1

    if matches_found == 0:
        print("No RNA IDs were matched with the Master FASTA.")
        return

    # 5. Run WebLogo
    print("4. Executing WebLogo3")
    try:
        command = [
            "weblogo", 
            "-f", fasta_out_path, 
            "-o", logo_path, 
            "-F", "png", 
            "--resolution", "300", 
            "--title", f"Top {percentage}% ipTM - {protein}",
            "--title-fontsize", "10",
            "--size", "large",
            "--units", "probability",    
            "--sequence-type", "rna"     
        ]
        subprocess.run(command, check=True)
        print(f"=== SUCCESS ===")
        print(f"Files created:")
        print(f" - {csv_out_path}")
        print(f" - {hist_path}")
        print(f" - {logo_path}")
    except Exception as e:
        print(f"Error executing WebLogo: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WebLogo and Histogram Analysis based on ipTM")
    parser.add_argument("-i", "--input", required=True, help="Input directory")
    parser.add_argument("-p", "--protein", required=True, help="Protein name")
    parser.add_argument("-pct", "--percentage", type=float, required=True, help="Top percentage to filter")
    
    args = parser.parse_args()
    generate_iptm_weblogo(args.input, args.protein, args.percentage)
