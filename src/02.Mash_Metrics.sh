#!/bin/bash

# Default kmers values
kmers=(7 9 11 12 13)
genus=""

while getopts "k:g:" opt; do
    case ${opt} in
        k )
            # Parse kmers as a space-separated string into an array
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

# Create variables for paths
parent_dir=$(dirname "$PWD")
source="$parent_dir/phallett/data/Taxa_Selected"
mkdir -p "$parent_dir/phallett/test/Metrics_Results"
outdir="$parent_dir/phallett/test/Metrics_Results"

# Create a list of subdirectories within the working directory
if [ -n "$genus" ]; then
  # If a genus name is provided, only process that genus
  subdirs=("$source/$genus")
else
  # Otherwise, process all subdirectories
  subdirs=($(find "$source" -mindepth 1 -type d))
fi

# Create a single signatures directory
sig_subdir_basename="signatures"
mkdir -p "${outdir}/${sig_subdir_basename}"

# Loop through each subdirectory
for subdir in "${subdirs[@]}"; do
  subdir_basename=$(basename "$subdir")
  genus_outdir="${outdir}/${sig_subdir_basename}/${subdir_basename}"
  mkdir -p "$genus_outdir" # Create a directory for each genus
  echo "Processing directory: $subdir"

  # Ensure we are in the correct directory
  cd "$subdir" || { echo "Failed to cd into $subdir"; exit 1; }

  # MASH ALGORITHM
  for k in "${kmers[@]}"; do
    sketch_file="${subdir_basename}_sketch_mash_k${k}"
    if ls *.fasta 1> /dev/null 2>&1; then
      mash sketch -k "${k}" -o "${sketch_file}" *.fasta
      echo "Sketching complete for $subdir_basename with kmer $k"

      if [ -f "${sketch_file}.msh" ]; then
        mash dist "${sketch_file}.msh" "${sketch_file}.msh" > "${subdir_basename}_mash_distance_k${k}.tab"
        echo "Mash distance matrix calculated for $subdir_basename with kmer $k"
      else
        echo "Mash sketch file ${sketch_file}.msh does not exist"
      fi
    else
      echo "No fasta files found in ${subdir}"
    fi
  done

  # Move mash results
  mv sketch* "${outdir}"
  mv *.tab "${outdir}"

  # SOURMASH ALGORITHM
  for k in "${kmers[@]}"; do
    echo "Running sourmash sketch for k-mer size ${k}"
    sourmash sketch dna -p k=${k} *.fasta -o "${subdir_basename}_k${k}.sig"
    if [ $? -ne 0 ]; then
      echo "Sourmash sketching failed for ${subdir_basename} with kmer $k"
      continue
    fi

    # Verify that .sig files are generated
    sig_files=("${subdir_basename}_k${k}.sig")
    if [ -f "${sig_files}" ]; then
      echo "Sourmash sketching complete for ${subdir_basename} with kmer $k"
    else
      echo "No signature files generated for ${subdir_basename} with kmer $k"
    fi
  done

  # Move the signatures to the newly created subdirectory
  mv *.sig "$genus_outdir/"
  echo "All signatures moved to $genus_outdir"

  # Run sourmash compare for each k-mer value
  for k in "${kmers[@]}"; do
    sig_files=("$genus_outdir"/*.sig)
    if [ ${#sig_files[@]} -gt 0 ]; then
      echo "Comparing signatures with kmer $k for $subdir_basename"
      sourmash compare --ksize ${k} --csv "${outdir}/sourmash_distance_${subdir_basename}_k${k}.csv" "${sig_files[@]}"
      if [ $? -eq 0 ]; then
        echo "Sourmash compare matrix calculated for $subdir_basename with kmer $k"
      else
        echo "Sourmash compare failed for $subdir_basename with kmer $k"
      fi
    else
      echo "No signature files found in $genus_outdir for kmer $k"
    fi
  done

  # Move sourmash results
  mv sourmash* "${outdir}"

  # Return to the parent directory
  cd - || { echo "Failed to return to the original directory"; exit 1; }
done
