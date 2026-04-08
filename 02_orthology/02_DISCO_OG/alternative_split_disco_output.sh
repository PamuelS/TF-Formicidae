#!/bin/bash

echo "Splitting DISCO outputs ..."

mkdir split

# Split the multiline output in multiple outputs
for disco in *.nwk; do
	OGname=$(basename -s .nwk "$disco")
	split -d -a 2 -l 1 --additional-suffix=.nwk "$disco" split/"$OGname"_
done

mkdir ../disco_OG

# Define original sequence folder with $1
OG_folder=$(realpath $1)

echo "Recreating orthogroups ..."

cd split

for tree in *_*.nwk; do
    name="${tree%.nwk}"
    OG="${name%%_*}"  # Extract the part before the number
    output="../../disco_OG/$name.faa"

    # Extract unique sequence names from tree - MODIFICATA REGEX
    grep -o -E "[^(),:]+" "$tree" | grep "|" | sort -u > tmp_seqs.txt

    # Extract sequences in one go using awk - MODIFICATA ESTENSIONE .fa
    awk -v seqs=tmp_seqs.txt '
        BEGIN {
            while ((getline < seqs) > 0) wanted[$1] = 1
        }
        /^>/ {
            seq = substr($0, 2)
            keep = (seq in wanted)
        }
        keep { print }
    ' "$OG_folder/$OG.fa" >> "$output"
done

# Clean up temporary file
rm -f tmp_seqs.txt
