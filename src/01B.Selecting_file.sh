#!/bin/bash

blastpor=0.75
evalue=1e-5
parent_dir=$(dirname "$PWD")
file="$parent_dir/phallett/GCF_000836945.fasta"
updatedb="false"

# Parse command-line options
while getopts "m:f:b:e:u:" opt; do  
  case $opt in
    m)
      module=$OPTARG
      ;;
    f)
      file=$OPTARG
      ;;
    b)
      blastpor=$OPTARG
      ;;
    e)
      evalue=$OPTARG
      ;;
    u)
      updatedb=$OPTARG
      ;;  
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

# Now running the Python script
python3 "$parent_dir/phallett/src/01B.Selecting_file.py" -file "$file" -updatedb "$updatedb" -blastpor "$blastpor" -evalue "$evalue"
