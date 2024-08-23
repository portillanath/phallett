#!/bin/bash
#For Memberships analysis 
parent_dir=$(dirname "$PWD")

# Parse command-line options
while getopts "m:f:" opt; do  
  case $opt in
    m)
      module=$OPTARG
      ;;
    f)
      file=$OPTARG
      ;;  
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Shift processed options
shift $((OPTIND -1))

# The first argument is the parent directory
source_dir="$1"
shift

# Remaining arguments are folder names
folder_names=("$@")

# Check if parent directory is provided
if [ -z "$source_dir" ]; then
    echo "Error: Parent directory is not provided."
    echo "Usage: $0 [-m module] [-f file] parent_directory folder1 folder2 ..."
    exit 1
fi

# Check if at least one folder name is provided
if [ ${#folder_names[@]} -eq 0 ]; then
    echo "Error: No folder names provided."
    echo "Usage: $0 [-m module] [-f file] parent_directory folder1 folder2 ..."
    exit 1
fi

# Run the Python script with the parent directory and folder names as arguments
python3 "$parent_dir/src/06.Bubbles.py" "$source_dir" "${folder_names[@]}"