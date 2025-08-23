from Bio import AlignIO
import os
import sys
import pandas as pd

# get HOG id from command line argument
if len(sys.argv) < 2:
    print("Usage: python get_sites_proteins.py <HOG_id>")
    sys.exit(1)

hog = sys.argv[1]

# Load alignments
masked_aln = AlignIO.read(f"og_alignments/{hog}/{hog}_final_align_AA.aln", "fasta")
unmasked_aln = AlignIO.read(f"og_alignments/{hog}/{hog}_final_unmask_align_AA.aln", "fasta")

# get Dmel sequence name from alignments (starting with 'rna-')
target_sequences = [rec.id for rec in masked_aln if rec.id.startswith('rna-')]

#target_sequences = ["rna-NM_078782.3", "rna-NM_078854.2"]

# Selected sites (MEME)
meme_dt = pd.read_csv(f"MEME_FEL_out/{hog}_MEME.out", sep="\t")
# sites that have a p-value <= 0.05 & beta_positive > alpha
selected_sites = meme_dt[(meme_dt['p_value'] <= 0.05) & (meme_dt['beta_positive'] > meme_dt['alpha'])]['Site'].tolist()

#selected_sites = [1, 5, 10, 30, 31, 39, 1240]

def get_flanking_sequence(seq_with_gaps, site_pos, upstream=5, downstream=5):
    """
    Get flanking sequence around a site position in an alignment
    Returns (flanking_sequence, relative_position_of_site)
    """
    site_idx = site_pos - 1  # Convert to 0-based
    
    if site_idx >= len(seq_with_gaps) or seq_with_gaps[site_idx] == '-':
        return None, None
    
    # Find start and end positions for flanking region
    start_pos = max(0, site_idx - upstream)
    end_pos = min(len(seq_with_gaps), site_idx + downstream + 1)
    
    # Extract flanking sequence (with gaps)
    flanking_with_gaps = seq_with_gaps[start_pos:end_pos]
    
    # Remove gaps to get ungapped flanking sequence
    flanking_ungapped = flanking_with_gaps.replace('-', '')
    
    # Find the relative position of our site in the ungapped flanking sequence
    chars_before_site = 0
    for i in range(start_pos, site_idx):
        if seq_with_gaps[i] != '-':
            chars_before_site += 1
    
    if len(flanking_ungapped) == 0:
        return None, None
    
    return flanking_ungapped, chars_before_site

def find_flanking_in_unmasked(flanking_seq, unmasked_ungapped, relative_site_pos):
    """
    Find the flanking sequence in the unmasked ungapped sequence
    Returns the absolute position of the site in the unmasked sequence
    """
    if not flanking_seq or relative_site_pos is None:
        return None
    
    # Find all occurrences of the flanking sequence
    positions = []
    start = 0
    while True:
        pos = unmasked_ungapped.find(flanking_seq, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    
    if len(positions) == 0:
        print(f"    Flanking sequence not found: {flanking_seq}")
        return None
    elif len(positions) == 1:
        # Unique match found
        site_pos_in_unmasked = positions[0] + relative_site_pos + 1  # Convert to 1-based
        print(f"    Unique match found for flanking: {flanking_seq}")
        return site_pos_in_unmasked
    else:
        # Multiple matches - return the first one but warn
        site_pos_in_unmasked = positions[0] + relative_site_pos + 1  # Convert to 1-based
        print(f"    Multiple matches found for flanking: {flanking_seq} (using first match)")
        print(f"    Match positions: {positions}")
        return site_pos_in_unmasked


def find_flanking_in_unmasked_fuzzy(flanking_seq, unmasked_ungapped, relative_site_pos, max_mismatches=1):
    """
    Find the flanking sequence in the unmasked ungapped sequence with fuzzy matching
    """
    if not flanking_seq or relative_site_pos is None:
        return None
    
    # First try exact match
    positions = []
    start = 0
    while True:
        pos = unmasked_ungapped.find(flanking_seq, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    
    if len(positions) == 1:
        site_pos_in_unmasked = positions[0] + relative_site_pos + 1
        print(f"    Exact match found")
        return site_pos_in_unmasked
    elif len(positions) > 1:
        site_pos_in_unmasked = positions[0] + relative_site_pos + 1
        print(f"    Multiple exact matches found (using first)")
        return site_pos_in_unmasked
    
    # If no exact match, try fuzzy matching
    print(f"    No exact match, trying fuzzy matching (max {max_mismatches} mismatches)")
    
    best_matches = []
    
    for i in range(len(unmasked_ungapped) - len(flanking_seq) + 1):
        window = unmasked_ungapped[i:i + len(flanking_seq)]
        mismatches = sum(1 for a, b in zip(flanking_seq, window) if a != b)
        
        if mismatches <= max_mismatches:
            best_matches.append((i, mismatches, window))
    
    if best_matches:
        # Sort by number of mismatches
        best_matches.sort(key=lambda x: x[1])
        best_pos, best_mismatches, best_window = best_matches[0]
        
        site_pos_in_unmasked = best_pos + relative_site_pos + 1
        print(f"    Fuzzy match found with {best_mismatches} mismatches:")
        print(f"    Expected: {flanking_seq}")
        print(f"    Found:    {best_window}")
        
        if len(best_matches) > 1 and best_matches[1][1] == best_mismatches:
            print(f"    Warning: Multiple matches with same mismatch count")
        
        return site_pos_in_unmasked
    
    return None

def try_different_flanking_sizes(seq_with_gaps, site_pos, unmasked_ungapped):
    """
    Try different flanking sizes to find a unique match
    """
    for upstream, downstream in [(2,2),(3, 3), (5, 5), (7, 7), (10, 10), (15, 15)]:
        print(f"    Trying flanking size: {upstream} upstream, {downstream} downstream")
        
        flanking_seq, relative_pos = get_flanking_sequence(
            seq_with_gaps, site_pos, upstream, downstream
        )
        
        if flanking_seq is None:
            continue
        
        print(f"    Flanking sequence: {flanking_seq} (site at position {relative_pos})")
        
        # Count occurrences
        count = unmasked_ungapped.count(flanking_seq)
        
        if count == 1:
            # Unique match found
            match_pos = find_flanking_in_unmasked_fuzzy(flanking_seq, unmasked_ungapped, relative_pos)
            print(f"    SUCCESS: Unique match with flanking size {upstream},{downstream}")
            return match_pos
        elif count == 0:
            print(f"    No match found")
        else:
            print(f"    {count} matches found - trying larger flanking region")
    
    print(f"    FAILED: Could not find unique match with any flanking size")
    return None



def try_different_flanking_sizes_improved(seq_with_gaps, site_pos, unmasked_ungapped):
    """
    Try different flanking sizes with improved edge case handling
    """
    # Get the position in the ungapped masked sequence
    masked_ungapped = seq_with_gaps.replace('-', '')
    site_idx = site_pos - 1
    
    # Count non-gap characters before this position
    ungapped_pos = 0
    for i in range(site_idx):
        if seq_with_gaps[i] != '-':
            ungapped_pos += 1
    
    ungapped_pos += 1  # Convert to 1-based
    
    print(f"    Site {site_pos} corresponds to position {ungapped_pos} in ungapped masked sequence")
    print(f"    Masked sequence length: {len(masked_ungapped)}")
    print(f"    Unmasked sequence length: {len(unmasked_ungapped)}")
    
    # Show the character at this position
    if ungapped_pos <= len(masked_ungapped):
        site_char = masked_ungapped[ungapped_pos - 1]
        print(f"    Character at site: '{site_char}'")
    
    # Try different flanking sizes, adapting to sequence edges
    max_upstream = min(15, ungapped_pos - 1)
    max_downstream = min(15, len(masked_ungapped) - ungapped_pos)
    
    print(f"    Available: {max_upstream} upstream, {max_downstream} downstream")
    
    for base_size in [3, 5, 7, 10, 15]:
        upstream = min(base_size, max_upstream)
        downstream = min(base_size, max_downstream)
        
        if upstream + downstream < 4:  # Need minimum flanking context
            continue
        
        print(f"    Trying flanking size: {upstream} upstream, {downstream} downstream")
        
        flanking_seq, relative_pos = get_flanking_sequence(
            seq_with_gaps, site_pos, upstream, downstream
        )
        
        if flanking_seq is None:
            print(f"    Could not extract flanking sequence")
            continue
        
        print(f"    Flanking sequence: {flanking_seq} (site at position {relative_pos})")
        
        # Try exact match first, then fuzzy
        match_pos = find_flanking_in_unmasked_fuzzy(
            flanking_seq, unmasked_ungapped, relative_pos, max_mismatches=1
        )
        
        if match_pos:
            print(f"    SUCCESS: Match found with flanking size {upstream},{downstream}")
            return match_pos
    
    # If still no match, try just the site with more context
    print(f"    Trying single character with maximum context")
    
    upstream = max_upstream
    downstream = max_downstream
    
    if upstream + downstream >= 6:  # Need reasonable context
        flanking_seq, relative_pos = get_flanking_sequence(
            seq_with_gaps, site_pos, upstream, downstream
        )
        
        if flanking_seq:
            print(f"    Final attempt - Flanking: {flanking_seq} (site at {relative_pos})")
            match_pos = find_flanking_in_unmasked_fuzzy(
                flanking_seq, unmasked_ungapped, relative_pos, max_mismatches=2
            )
            
            if match_pos:
                print(f"    SUCCESS: Match found with relaxed criteria")
                return match_pos
    
    print(f"    FAILED: Could not find match with any approach")
    return None

# Main processing
results = {}

for target_seq_id in target_sequences:
    print(f"\nProcessing {target_seq_id}...")
    
    # Find the sequence in both alignments
    masked_seq = None
    unmasked_seq = None
    
    for rec in masked_aln:
        if rec.id == target_seq_id:
            masked_seq = str(rec.seq)
            break
    
    for rec in unmasked_aln:
        if rec.id == target_seq_id:
            unmasked_seq = str(rec.seq)
            break
    
    if masked_seq is None or unmasked_seq is None:
        print(f"  Sequence {target_seq_id} not found in one or both alignments!")
        continue
    
    # Get ungapped unmasked sequence
    unmasked_ungapped = unmasked_seq.replace('-', '')
    
    print(f"  Masked length: {len(masked_seq)}")
    print(f"  Unmasked ungapped length: {len(unmasked_ungapped)}")
    #print(unmasked_ungapped)  
    
    # Process each selected site
    site_mappings = {}
    
    for site in selected_sites:
        print(f"\n  Processing site {site}...")
        
        # Check if site exists and is not a gap
        if site > len(masked_seq) or masked_seq[site - 1] == '-':
            print(f"    Site {site} is gap or out of bounds in masked alignment")
            site_mappings[site] = None
            continue
        
        # Try to find unique match using flanking sequences
        unmasked_pos = try_different_flanking_sizes(masked_seq, site, unmasked_ungapped)
        site_mappings[site] = unmasked_pos
        
        if unmasked_pos:
            print(f"    Site {site} maps to position {unmasked_pos} in unmasked protein")
        else:
            print(f"    Site {site} could not be mapped")
    
    results[target_seq_id] = site_mappings

# Summary output
print("\n" + "="*60)
print("FINAL RESULTS SUMMARY")
print("="*60)

# Store the selected sites and their mappings in a file
write_file = f"sites_mappings_MEME/{hog}_site_mappings.tsv"
with open(write_file, 'w') as f:
    f.write("Sequence_ID\tMasked_Site\tUnmasked_Position\tp_value\tbranches\n")
    for seq_id, mappings in results.items():
        print(f"\n{seq_id}:")
        # create a chimera attribute file 
        write_file = f"sites_mappings_MEME/{hog}_{seq_id}.defattr"
        with open(write_file, 'w') as f2:
            f2.write("attribute: branch_count\nrecipient: residues\n")
            for site in sorted(selected_sites):
                p_value = meme_dt[meme_dt['Site'] == site]['p_value'].values[0] if site in meme_dt['Site'].values else 'N/A'
                branches = meme_dt[meme_dt['Site'] == site]['branches_under_selection'].values[0] if site in meme_dt['Site'].values else 'N/A'
                if site in mappings:
                    pos = mappings[site]
                    if pos is not None:
                        print(f"  Masked alignment site {site} => Unmasked protein position {pos}")
                        #with open(write_file, 'a') as f:
                        f.write(f"{seq_id}\t{site}\t{pos}\t{p_value}\t{branches}\n")
                        f2.write(f"\t:{pos}\t{branches}\n")
                    else:
                        print(f"  Masked alignment site {site} => Could not map")
                        #with open(write_file, 'a') as f:
                        f.write(f"{seq_id}\t{site}\tN/A\t{p_value}\t{branches}\n")

# Verification (show the actual characters at mapped positions)
print("\n" + "="*60)
print("VERIFICATION - Characters at mapped positions")
print("="*60)

for seq_id, mappings in results.items():
    print(f"\n{seq_id}:")
    
    # Get sequences
    masked_seq = None
    unmasked_seq = None
    
    for rec in masked_aln:
        if rec.id == seq_id:
            masked_seq = str(rec.seq)
            break
    
    for rec in unmasked_aln:
        if rec.id == seq_id:
            unmasked_seq = str(rec.seq)
            break
    
    if masked_seq and unmasked_seq:
        unmasked_ungapped = unmasked_seq.replace('-', '')
        
        for site in sorted(selected_sites):
            if site in mappings and mappings[site] is not None:
                pos = mappings[site]
                
                # Get characters
                masked_char = masked_seq[site - 1] if site <= len(masked_seq) else 'N/A'
                unmasked_char = unmasked_ungapped[pos - 1] if pos <= len(unmasked_ungapped) else 'N/A'
                
                match = "✓" if masked_char == unmasked_char else "✗"
                print(f"  Site {site}: '{masked_char}' => '{unmasked_char}' (pos {pos}) {match}")
                
                # Show context
                start_ctx = max(0, pos - 6)
                end_ctx = min(len(unmasked_ungapped), pos + 5)
                context = unmasked_ungapped[start_ctx:end_ctx]
                marker = " " * (pos - start_ctx - 1) + "^"
                print(f"    Context: {context}")
                print(f"             {marker}")