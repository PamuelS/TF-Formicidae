import os
import re
import csv
from pathlib import Path

# Path to your OG folder
OG_DIR = "../../02_orthology/02_DISCO_OG"
OUTPUT_FILE = "orthogroups_top_annotations.tsv"

def score_annotation(header_text):
    """
    Scores how 'annotated' a sequence header is.
    Higher score = better annotation.
    """
    # Exclude generic or uninformative descriptions
    uninformative_patterns = [
        r"hypothetical", r"uncharacterized", r"predicted protein", 
        r"unknown", r"unnamed", r"putative protein"
    ]
    
    # Check if header contains generic keywords
    for pattern in uninformative_patterns:
        if re.search(pattern, header_text, re.IGNORECASE):
            return 1  # Low score for uncharacterized/hypothetical
            
    # Otherwise, score based on description length or presence of functional terms
    return len(header_text)

def parse_og_file(filepath):
    """
    Parses a single OG FASTA file and returns:
    (best_protein_id, best_species_name)
    """
    best_score = -1
    best_protein = "N/A"
    best_species = "N/A"

    with open(filepath, 'r') as f:
        current_header = ""
        for line in f:
            if line.startswith(">"):
                current_header = line.strip()[1:]  # Remove leading '>'
                
                # --- PARSING HEADER (Adjust split logic based on your headers) ---
                # Example assumption: >ProteinID SpeciesName description...
                parts = current_header.split(maxsplit=2)
                
                protein_id = parts[0] if len(parts) > 0 else "Unknown"
                
                # Extract species (Assumes species is 2nd word or encoded in ID)
                species_name = parts[1] if len(parts) > 1 else "Unknown"
                
                # Evaluate annotation quality
                score = score_annotation(current_header)
                
                if score > best_score:
                    best_score = score
                    best_protein = protein_id
                    best_species = species_name

    return best_protein, best_species

def main():
    og_path = Path(OG_DIR)
    
    if not og_path.exists():
        print(f"Error: Directory '{OG_DIR}' not found.")
        return

    # Find all FASTA files (.fa, .fasta, .txt, or files without extension)
    og_files = sorted([f for f in og_path.iterdir() if f.is_file() and not f.name.startswith(".")])

    print(f"Found {len(og_files)} OG files in {OG_DIR}. Processing...")

    with open(OUTPUT_FILE, 'w', newline='') as out_tsv:
        writer = csv.writer(out_tsv, delimiter='\t')
        # Write Header
        writer.writerow(["Orthogroup", "Most_Annotated_Protein", "Species"])

        for filepath in og_files:
            og_name = filepath.stem  # Extracts 'OG0000001' from 'OG0000001.fa'
            protein, species = parse_og_file(filepath)
            writer.writerow([og_name, protein, species])

    print(f"Done! Results written to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
