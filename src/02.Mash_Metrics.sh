#!/bin/bash

# Default k-mer values
kmers=(7 9 11 12 13)
genus=""

# Parse command-line options
while getopts "k:g:" opt; do
    case ${opt} in
        k )
            # Parse k-mers as a comma-separated string into an array
            IFS=',' read -r -a kmers <<< "${OPTARG}"
            ;;
        g )
            genus=${OPTARG}
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            exit 1
            ;;
    esac
done

# Define directories
parent_dir=$(dirname "$PWD")
source="$parent_dir/phallett/data/Taxa_Selected"
outdir="$parent_dir/phallett/test/Metrics_Results"
mkdir -p "$outdir"

# Create a list of subdirectories within the source directory
if [ -n "$genus" ]; then
  # If a genus name is provided, only process that genus
  subdirs=("$source/$genus")
else
  # Otherwise, process all subdirectories
  subdirs=($(find "$source" -mindepth 1 -type d))
fi

# Process each subdirectory
for subdir in "${subdirs[@]}"; do
  subdir_basename=$(basename "$subdir")
  genus_outdir="${outdir}/${subdir_basename}/signatures"
  mkdir -p "$genus_outdir" # Create directory for signatures
  echo "Processing directory: $subdir"

  # Ensure we are in the correct directory
  cd "$subdir" || { echo "Failed to cd into $subdir"; exit 1; }

  # Generate Sourmash signatures for each k-mer size
  for k in "${kmers[@]}"; do
    echo "Running sourmash sketch for k-mer size ${k}"

    # Process each FASTA file individually
    for fasta_file in *.fasta; do
      if [ -f "$fasta_file" ]; then
        # Extract base name of the FASTA file (e.g., genome1.fasta -> genome1)
        base_name=$(basename "$fasta_file" .fasta)
        sig_file="${base_name}_k${k}.sig"
        
        sourmash sketch dna -p k=${k} "$fasta_file" --singleton -o "$sig_file"
        
        if [ $? -eq 0 ]; then
          echo "Signature created: $sig_file"
        else
          echo "Error creating signature for $fasta_file with k-mer size $k"
          continue
        fi

        # Move signature to output directory
        mv "$sig_file" "$genus_outdir/"
      else
        echo "No FASTA files found in $subdir"
      fi
    done
  done

  # Run sourmash compare for each k-mer size
  for k in "${kmers[@]}"; do
    sig_files=("$genus_outdir"/*_k${k}.sig)
    if [ ${#sig_files[@]} -gt 0 ]; then
      echo "Comparing signatures with k-mer size $k for $subdir_basename"
      
      sourmash compare --ksize ${k} --csv "${outdir}/sourmash_distance_${subdir_basename}_k${k}.csv" "${sig_files[@]}"
      
      if [ $? -eq 0 ]; then
        echo "Sourmash compare matrix calculated for $subdir_basename with k-mer $k"
      else
        echo "Sourmash compare failed for $subdir_basename with k-mer $k"
      fi
    else
      echo "No signature files found in $genus_outdir for k-mer $k"
    fi
  done

  # Return to the parent directory
  cd - || { echo "Failed to return to the original directory"; exit 1; }
done
