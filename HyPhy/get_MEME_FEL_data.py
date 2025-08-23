#! /usr/bin/env python3
# Read Json files from MEME and FEL directories and stores site-level data in separate files
import json
import csv
import sys
import argparse

# example usage:
# python get_MEME_FEL_data.py input.json output.tsv --method MEME
# python get_MEME_FEL_data.py input.json output.tsv --method FEL
def parse_fel_json(json_file, output_file):
    """Parse FEL JSON and convert to TSV format"""
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract site data from MLE content
    site_data = data["MLE"]["content"]
    
    # Define column headers based on FEL output
    headers = [
        "Site",
        "alpha",
        "beta", 
        "alpha_equals_beta",
        "LRT",
        "p_value",
        "total_branch_length"
    ]
    
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(headers)
        
        # Get the number of sites from the first parameter array
        num_sites = len(site_data["0"])
        
        # Write data for each site
        for site_idx in range(num_sites):
            site_values = [site_idx + 1] + site_data["0"][site_idx]
            writer.writerow(site_values)

def parse_meme_json(json_file, output_file):
    """Parse MEME JSON and convert to TSV format"""
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract site data from MLE content
    site_data = data["MLE"]["content"]
    
    # Define column headers based on MEME output
    headers = [
        "Site",
        "alpha",
        "beta_negative", 
        "p_negative",
        "beta_positive",
        "p_positive", 
        "LRT",
        "p_value",
        "branches_under_selection",
        "total_branch_length",
        "MEME_LogL",
        "FEL_LogL"
    ]
    
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(headers)
        
        # Get the number of sites from the first parameter array
        num_sites = len(site_data["0"])
        
        # Write data for each site
        for site_idx in range(num_sites):
            site_values = [site_idx + 1] + site_data["0"][site_idx]
            writer.writerow(site_values)

def main():
    parser = argparse.ArgumentParser(description='Convert HyPhy MEME/FEL JSON to TSV')
    parser.add_argument('input_json', help='Input JSON file')
    parser.add_argument('output_tsv', help='Output TSV file')
    parser.add_argument('--method', choices=['MEME', 'FEL'], required=True,
                        help='Analysis method (MEME or FEL)')
    
    args = parser.parse_args()
    
    if args.method == 'MEME':
        parse_meme_json(args.input_json, args.output_tsv)
        #print(f"MEME data converted to {args.output_tsv}")
    elif args.method == 'FEL':
        parse_fel_json(args.input_json, args.output_tsv)
        #print(f"FEL data converted to {args.output_tsv}")

if __name__ == "__main__":
    main()