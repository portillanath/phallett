#!/bin/bash

# Create a directory for taxa selected
parent_dir=$(dirname "$PWD")
mkdir -p "$parent_dir/phallett/data/Taxa_Selected"

# Extract --all-hosts flag if present
all_hosts_flag=""
genus_names=()

for arg in "$@"; do
    if [[ "$arg" == "--all-hosts" ]]; then
        all_hosts_flag="--all-hosts"
    else
        genus_names+=("$arg")
    fi
done

# Run Python script for each genus with the optional flag
for genus_name in "${genus_names[@]}"; do
    python3 "$parent_dir/phallett/src/01A.Taxa_Curation_Level.py" "$genus_name" $all_hosts_flag
done
