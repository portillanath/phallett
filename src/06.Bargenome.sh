#!/bin/bash
# For Memberships analysis 

# Get the parent directory (one level up from the current working directory)
parent_dir=$(dirname "$PWD")

# Run the Python script with the fixed parent directory path
python3 "$parent_dir/src/06.Bargenome.py" "$parent_dir"
