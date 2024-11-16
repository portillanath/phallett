import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  
from pathlib import Path
from Bio import SeqIO
import argparse

# Set the main directory
current_dir = Path.cwd() 
parent_dir = current_dir.parent
main_directory = (Path(parent_dir) / "data" / "Taxa_Selected")

# Initialize empty lists to store genus names, genome sizes, and file counts
genus_names = []
genome_sizes = []
genus_file_counts = []

# Function to calculate the genome size of a FASTA file
def get_genome_size(fasta_file_path):
    genome_size = 0
    for record in SeqIO.parse(fasta_file_path, "fasta"):
        genome_size += len(record.seq)
    return genome_size

# Iterate through subdirectories
for subdir in os.listdir(main_directory):
    # Create the full path to the subdirectory
    subdirectory = os.path.join(main_directory, subdir)

    # List FASTA files in the subdirectory
    fasta_files = [file for file in os.listdir(subdirectory) if file.endswith(".fasta")]
    
    # Get the genus name from the subdirectory name
    genus_name = os.path.basename(subdirectory)

    # Initialize list for storing genome sizes for this genus
    genus_genome_sizes = []

    # Iterate through the FASTA files and calculate genome sizes
    for fasta_file in fasta_files:
        fasta_file_path = os.path.join(subdirectory, fasta_file)
        genome_size = get_genome_size(fasta_file_path)
        genus_genome_sizes.append(genome_size)
    
    # Store the genus name and genome sizes
    genus_names.extend([genus_name] * len(genus_genome_sizes))  # Repeat genus name for each genome size
    genome_sizes.extend(genus_genome_sizes)

    # Store the count of FASTA files for this genus
    genus_file_counts.extend([len(fasta_files)] * len(genus_genome_sizes))  # Repeat count for each genome size

# Create a DataFrame with the collected data
data = pd.DataFrame({
    "Genus": genus_names,
    "Genome Size": genome_sizes,
    "Count": genus_file_counts
})

# Filter genera with more than 10 genomes
data_filtered = data[data["Count"] > 10]

# Save the filtered data as a CSV file
csv_file_path = (Path(parent_dir) / "Summary_feed_filtered.csv")
data_filtered.to_csv(csv_file_path, index=False)
print(f"Data saved as CSV: {csv_file_path}")

# Plot boxplot for genome sizes per genus
plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")
ax = sns.boxplot(x="Genus", y="Genome Size", data=data_filtered, palette="Set2")
plt.title("Genome Sizes by Genus")
plt.xlabel("Genus")
plt.ylabel("Genome Size (bp)")
plt.xticks(rotation=90, ha="right", fontsize=10)

# Add annotations for the number of genomes (FASTA files) per genus
for i, genus in enumerate(data_filtered["Genus"].unique()):
    # Find the number of genomes (FASTA files) for the current genus
    count = data_filtered[data_filtered["Genus"] == genus]["Count"].iloc[0]
    
    # Annotate the plot with the count value, positioned above each box
    ax.text(i, max(data_filtered[data_filtered["Genus"] == genus]["Genome Size"]) * 1.05, 
            f"n={count}", ha='center', fontsize=10, color='black')

plt.tight_layout()

# Save the plot as a PDF file
plt.savefig("genome_sizes_by_genus_boxplot_with_count.pdf", format="pdf")
print("Plot saved as PDF: genome_sizes_by_genus_boxplot_with_count.pdf")

