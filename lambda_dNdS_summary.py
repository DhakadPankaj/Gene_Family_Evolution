# This script makes a tsv file with the dN/dS and Lambda values for each HOG (Immune and Control)

import re
import json
import pandas as pd
import argparse
import os
from Bio import SeqIO
import numpy as np
from concurrent.futures import ThreadPoolExecutor

# Set up argument parser
parser = argparse.ArgumentParser(description='Extract dN/dS and p-values from a HyPhy BUSTED JSON output file. Lambda values are also extracted from the CAFE output files.')
parser.add_argument('control_desc', type=str, help='Path to the "control_genes.tsv" file')
parser.add_argument('immune_desc', type=str, help='Path to the "HOGs_immune_genes.tsv" file')
parser.add_argument('json_folder', type=str, help='Path to the folder containing JSON files') # HyPhy_BUSTED
parser.add_argument('Lambda', type=str, help='Path to the file containing lambda of all HOGs') # CAFE/immune_k2/Gamma_lambda_per_family.txt
#parser.add_argument('Lambda_control', type=str, help='Path to the file containing lambda of Control HOGs') # CAFE/control_k2/Gamma_lambda_per_family.txt
parser.add_argument('fasta_folder', type=str, help='Path to the folder containing multi-FASTA files') # ~/fly_annotation/complement_annotations/Results_MSA/CDS_HOG/

# Example usage
# python lambda_dNdS_summary.py control_genes.tsv HOGs_immune_genes.tsv HyPhy_BUSTED CAFE/immune_k2/Gamma_lambda_per_family.txt CAFE/control_k2/Gamma_lambda_per_family.txt ~/fly_annotation/complement_annotations/Results_MSA/CDS_HOG

# Parse arguments
args = parser.parse_args()

# Read the TSV files
immune_control_df = pd.read_csv(args.control_desc, sep='\t')
immune_HOGs_df = pd.read_csv(args.immune_desc, sep='\t', encoding='ISO-8859-1')
#lambda_immune_df = pd.read_csv(args.Lambda_immune, sep='\s+', engine='python')
#lambda_control_df = pd.read_csv(args.Lambda_control, sep='\s+', engine='python')
lambda_df = pd.read_csv(args.Lambda, sep='\s+', engine='python')

# take the mean of the lambda values for each HOG
lambda_df = lambda_df.groupby('HOG', as_index=False).mean()

# Function to extract dN/dS and p-values from a HyPhy output file
def extract_dn_ds_pvalues(output_file):
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Extract dN/dS ratio from the full codon model section
        dn_ds_match = re.search(r'### Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model.*?non-synonymous/synonymous rate ratio for \*test\* =\s+([\d.]+)', content, re.DOTALL)
        dn_ds = dn_ds_match.group(1) if dn_ds_match else 'NA'
        
        # Extract p-value for evidence of diversifying selection
        p_value_match = re.search(r'Likelihood ratio test for episodic diversifying positive selection, \*\*p =\s+([\d.]+)', content)
        p_value = p_value_match.group(1) if p_value_match else 'NA'

        # Extract proportion of sites under diversifying selection
        prop_match = re.search(r'Diversifying selection\s+\|\s+[\d.]+\s+\|\s+([\d.]+)', content)
        if prop_match:
            proportion = float(prop_match.group(1))
            prop = proportion
        else:
            prop = 'NA'

    except (FileNotFoundError, AttributeError):
        dn_ds, p_value = 'NA', 'NA'
    
    return dn_ds, p_value, prop

# Function to calculate the median length of sequences in a multi-FASTA file
def calculate_median_length(fasta_file):
    lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))
    return np.median(lengths) if lengths else 'NA'

# Initialize a list to store the results
#results = []
new_table = []

def process_row(row):
    results = []
    gene_hog = row['Immune_OG']
    lambda_val = lambda_df.loc[lambda_df['HOG'] == gene_hog]
    if not lambda_val.empty:
        lambda_val = lambda_val['Lambda'].values[0]
    else:
        lambda_val = 'NA'
    if pd.isna(gene_hog) or gene_hog == 'NA':
        return results

    # Get Immune Class/pathways and all other information
    immune_info = immune_HOGs_df[immune_HOGs_df['HOG'] == gene_hog]
    if not immune_info.empty:
        for col in immune_info.columns:
            if col != 'HOG':
                row[col] = immune_info.iloc[0][col]
                # Extract all other columns from immune_HOGs_df using matching gene_hog and store them in a variable
                immune_info_dict = immune_info.iloc[0].to_dict()
                immune_info_dict.pop('HOG', None)  # Remove the 'HOG' key as it's not needed

    json_file = os.path.join(args.json_folder, f'{gene_hog}_output.txt')
    dn_ds, p_value, prop = extract_dn_ds_pvalues(json_file)

    fasta_file = os.path.join(args.fasta_folder, f'{gene_hog}_CDS.fa')
    if not os.path.exists(fasta_file):
        median_length = "NA"
    else:
        median_length = calculate_median_length(fasta_file)

    results.append({
        'HOG': gene_hog,
        'Type': 'Immune',
        'median_length': median_length,
        'dNdS': dn_ds,
        'pvalue': p_value,
        'proportion': prop,
        'lambda': lambda_val,
        'group': gene_hog,
        **immune_info_dict
    })

    control_hogs = {}
    if isinstance(row['Control_OGs'], str) and row['Control_OGs'] != 'NA':
        hogs = row['Control_OGs'].split(',')
        lambdas = []
        for hog in hogs:
            lambda_val = lambda_df.loc[lambda_df['HOG'] == hog]
            if not lambda_val.empty:
                lambdas.append(lambda_val['Lambda'].values[0])
            else:
                lambdas.append('NA')
        control_hogs = dict(zip(hogs, lambdas))

    for control_hog, control_lambda in control_hogs.items():
        control_json_file = os.path.join(args.json_folder, f'{control_hog}_output.txt')
        control_dn_ds, control_p_value, control_prop = extract_dn_ds_pvalues(control_json_file)

        control_fasta_file = os.path.join(args.fasta_folder, f'{control_hog}_CDS.fa')
        if not os.path.exists(control_fasta_file):
            control_median_length = "NA"
        else:
            control_median_length = calculate_median_length(control_fasta_file)

        # Add "NA" for the missing values (i.e. length of immune_info_dict)
        control_info_dict = {k: 'NA' for k in immune_info_dict.keys()}
        results.append({
            'HOG': control_hog,
            'Type': 'Control',
            'median_length': control_median_length,
            'dNdS': control_dn_ds,
            'pvalue': control_p_value,
            'proportion': control_prop,
            'lambda': control_lambda,
            'group': gene_hog,
            **control_info_dict
        })
    return results

with ThreadPoolExecutor(max_workers=40) as executor:
    futures = [executor.submit(process_row, row) for index, row in immune_control_df.iterrows()]
    for future in futures:
        new_table.extend(future.result())
    


# Convert new_table to a DataFrame and save to a new TSV file
new_table_df = pd.DataFrame(new_table)
new_table_df.to_csv('lambda_dNdS_summary.tsv', sep='\t', index=False)