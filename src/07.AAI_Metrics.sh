

# Create variables for paths
parent_dir=$(dirname "$PWD")
source="$parent_dir/data/Taxa_Selected"
outdir="$parent_dir/test/Metrics_Results"

# Create a list of subdirectories within the working directory
if [ -n "$genus" ]; then
  # If a genus name is provided, only process that genus
  subdirs=("$source/$genus")
else
  # Otherwise, process all subdirectories
  subdirs=($(find "$source" -mindepth 1 -type d))
fi

#Loop through each subdirectory
for subdir in "${subdirs[@]}"; do
  subdir_basename=$(basename "$subdir")
  cd "$subdir"

  comparem_output="${outdir}/${subdir_basename}/${subdir_basename}_compare_results"
  mkdir -p "${comparem_output}"

  comparem aai_wf ./ "${comparem_output}" --file_ext fasta

  done 


