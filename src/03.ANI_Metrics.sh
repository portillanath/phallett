#!/bin/bash

# Default kmer values
kmers=(12 11 10 9 8)
frag_lengths=(500)  # Default fragment length
genus=""

while getopts "k:f:g:" opt; do
  case $opt in
    k)
      kmers=($OPTARG)
      ;;
    f)
      frag_lengths=($OPTARG)
      ;; 
    g)
      genus=($OPTARG)
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Create variables for paths
parent_dir=$(dirname "$PWD")
source="$parent_dir/phallett/data/Taxa_Selected"
outdir="$parent_dir/phallett/test/Metrics_Results"

cd "$source"

subdirs=($(find "$source" -mindepth 1 -type d))
for subdir in "${subdirs[@]}" ; do
cd $subdir
genusname=${subdir##*/}
mv "${subdir}/signatures_${genusname}" "${outdir}/"
echo "The moving signatures were sucessful"
done

# Create a list of subdirectories within the working directory
if [ -n "$genus" ]; then
  # If a genus name is provided, only process that genus
  subdirs=("$source/$genus")
else
  # Otherwise, process all subdirectories
  subdirs=($(find "$source" -mindepth 1 -type d))
fi

for subdir in "${subdirs[@]}"; do
  genusname=${subdir##*/}
  fasta_files=("${subdir}"/*.fasta)
  ls -d "${subdir}"/*.fasta > "${genusname}.list"
  mv *.list "$outdir"

  # Now we are going to calculate the ANI value. The file generated is interpreted as:
  # first column Ref_file, Query_file, ANI, Align_fraction_ref, Align_fraction_query, Ref_name, Query_Name

  output_file="${outdir}/skani_distance_${genusname}.txt"
  echo -e "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name" > "$output_file"
  
  for i in "${!fasta_files[@]}"; do
    fasta1="${fasta_files[$i]}"

    for j in "${!fasta_files[@]}"; do
      fasta2="${fasta_files[$j]}"

      echo "Comparing $fasta1 with $fasta2"
      skani dist "$fasta1" "$fasta2" | tail -n +2 >> "$output_file"
      skani dist "$fasta2" "$fasta1" | tail -n +2 >> "$output_file"
    done
  done

  # Remove exact duplicate rows from the output file using awk
  awk '!seen[$0]++' "$output_file" > "${output_file}.tmp"
  mv "${output_file}.tmp" "$output_file"

  mv *.txt "$outdir"
  
for frag_len in "${frag_lengths[@]}"; do
    for k in "${kmers[@]}"; do
       #average_nucleotide_identity.py -i "${subdir_basename}.list" -o "${subdir}/fastani_${subdir_basename}_frag_${frag_len}_${k}" --method ANIb
      fastANI --ql "$outdir/${genusname}.list" --rl "$outdir/${genusname}.list" -o "${outdir}/fastani_${genusname}_frag_${frag_len}_${k}" --fragLen "${frag_len}" --kmer "${k}"
      if [ -f "${outdir}/fastani_${genusname}_frag_${frag_len}_${k}" ]; then
        echo "Fastani was calculated for $subdir_basename with kmer $k and fragment length $frag_len"
      else
        echo "Error: fastANI calculation failed for $subdir_basename with kmer $k and fragment length $frag_len"
      fi
    done
  done
done

