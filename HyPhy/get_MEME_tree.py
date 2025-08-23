import json
import argparse
import re
import os
import csv
import pandas as pd

def count_ebf_sites(branch_data, threshold, selected_sites):
    """
    Count EBF sites >= threshold for a branch
    """
    count = 0
    # match EBF sites in branch data. ex. "EBF site {site_selected} (partition 1)"
    for site_selected in selected_sites:
        site_key = f"EBF site {site_selected} (partition 1)"
        if site_key in branch_data:
            value = branch_data[site_key]
            if value >= threshold:
                count += 1
    #for key, value in branch_data.items():
    #    if key.startswith("EBF site") and value >= threshold:
    #        count += 1
    return count

def get_branch_length(branch_data):
    """
    Extract Global MG94xREV branch length
    """
    return branch_data.get("Global MG94xREV", 0.0)

def get_original_name(branch_data):
    """
    Extract original branch name
    """
    return branch_data.get("original name", "")

def annotate_tree_with_sites(tree_string, branch_info):
    """
    Add site counts to tree in format: taxon:length[sites=count]
    """
    annotated_tree = tree_string
    
    for name, info in branch_info.items():
        site_count = info['sites']
        branch_length = info['branch_length']
        
        # Find pattern: name:length and replace with name:length[sites=count]
        pattern = f'({re.escape(name)}):([0-9.e-]+)'
        replacement = f'\\1:{branch_length:.6f}[sites={site_count}]'
        annotated_tree = re.sub(pattern, replacement, annotated_tree)
    
    return annotated_tree

def parse_meme_json(json_file, output_file, threshold):
    """
    Parse MEME JSON file and create annotated tree
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Get the tree from "input" section
    tree = data["input"]["trees"]["0"]
    
    # Get branch attributes
    branch_attributes = data["branch attributes"]["0"]
    
    # Extract branch information
    branch_info = {}

    # sites under selection (p <= 0.05 and beta_positive > alpha) 
    filename = json_file.replace('HyPhy_MEME/', 'MEME_FEL_out/').replace('_MEME.json', '_MEME.out')
    if not os.path.isfile(filename):
        print(f"Warning: Site data file not found: {filename}. No sites will be selected.")
        selected_sites = set()
    else:
        site_data = pd.read_csv(filename, sep='\t')
        selected_sites = site_data[(site_data['p_value'] <= 0.05) & (site_data['beta_positive'] > site_data['alpha'])]['Site'].tolist()
        selected_sites = set(selected_sites)  # Use a set for faster lookups
    print(f"Selected sites: {selected_sites}")    
    for branch_name, branch_data in branch_attributes.items():
        #print(branch_name)
        site_count = count_ebf_sites(branch_data, threshold, selected_sites)
        branch_length = get_branch_length(branch_data)
        #original_name = get_original_name(branch_data)
        
        # Use original name if available, otherwise use branch name
        name_key = branch_name
        
        branch_info[name_key] = {
            'sites': site_count,
            'branch_length': branch_length
        }

    #print(branch_info)

    # Add branch lengths on the tree
    for name, info in branch_info.items():
        # Use word boundary to avoid partial matches (e.g., "Node96" vs "Node6")
        branch_length = info['branch_length']
        pattern = rf'({re.escape(name)}):[0-9.e-]+'
        replacement = f'\\1:{branch_length:.6f}'
        tree = re.sub(pattern, replacement, tree)
    
    print(tree)
    # Annotate tree with site counts
    annotated_tree = annotate_tree_with_sites(tree, branch_info)

    ## Write iTOL compatible file for number of sites under selection 
    itol_file = output_file.replace('.nw', '.itol')
    with open(itol_file, 'w') as f:
        # header
        f.write("METADATA\n\nSEPARATOR TAB\n\n")
        f.write("FIELD_LABELS\tSites\n\nDATA\n")
        for name, info in branch_info.items():
            f.write(f"{name}\t{info['sites']}\n")
    

    
    # Write output
    with open(output_file, 'w') as f:
        f.write(annotated_tree)
    
    print(f"Processed {json_file} -> {output_file}")
    print(f"Total branches with sites >= {threshold}: {sum(1 for info in branch_info.values() if info['sites'] > 0)}")
    print(f"Total sites under selection: {sum(info['sites'] for info in branch_info.values())}")

def get_output_filename(input_filename):
    # Remove _MEME.json suffix
    base_name = input_filename.replace('_MEME.json', '')
    
    # Add .nw extension
    output_name = base_name + '.nw'
    
    return output_name

def main():
    parser = argparse.ArgumentParser(description='Get annotated tree from HyPhy MEME JSON')
    parser.add_argument('input_json', nargs='+', help='Input JSON file(s)')
    parser.add_argument('--method', choices=['MEME'], required=True,
                        help='Analysis method (MEME)')
    parser.add_argument('--threshold', type=float, default=100,
                        help='EBF threshold for sites under selection (default: 100)')
    parser.add_argument('--output-dir', help='Output directory (default: same as input)')
    
    args = parser.parse_args()
    
    if args.method == 'MEME':
        for json_file in args.input_json:
            # Generate output filename
            if args.output_dir:
                base_name = os.path.basename(json_file)
                output_file = os.path.join(args.output_dir, get_output_filename(base_name))
            else:
                output_file = get_output_filename(json_file)
            
            try:
                                
                parse_meme_json(json_file, output_file, threshold=args.threshold)
                
            except Exception as e:
                print(f"Error processing {json_file}: {e}")
                continue

if __name__ == "__main__":
    main()
