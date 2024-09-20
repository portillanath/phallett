import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  
from pathlib import Path 
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO

# Set the main directory
current_dir = Path.cwd() 
parent_dir = current_dir.parent
main_directory = (Path(parent_dir) / "data" / "Taxa_Selected")

# Initialize empty lists to store genus names and genome counts
genus_names = []
genome_counts = []

# Iterate through subdirectories
for subdir in os.listdir(main_directory):
    # Create the full path to the subdirectory
    subdirectory = os.path.join(main_directory, subdir)

    # List FASTA files in the subdirectory
    fasta_files = [file for file in os.listdir(subdirectory) if file.endswith(".fasta")]

    # Get the genus name from the subdirectory name
    genus_name = os.path.basename(subdirectory)

    # Append genus name and genome count to lists
    genus_names.append(genus_name)
    genome_counts.append(len(fasta_files))

# Create a DataFrame
data = pd.DataFrame({"Genus": genus_names, "Count": genome_counts})

# Filter genera with more than 10 genomes
data_filtered = data[data["Count"] > 10]
print(data_filtered)

# Set the PDF file path and name
csv_file_path = (Path(parent_dir)/ "Summary_feed_filtered.csv")

# Save the filtered data as a CSV file
data_filtered.to_csv(csv_file_path, index=False)

print("Data saved as CSV:", csv_file_path)

# Use a different color palette for the plot (e.g., "Set1" from seaborn)
colors = sns.color_palette("Set1", len(data_filtered))

# Create a bar plot using matplotlib with the custom color palette
plt.figure(figsize=(10, 6))
plt.bar(data_filtered["Genus"], data_filtered["Count"], color=colors)
plt.title("Genome Counts by Genus (Genus with >10 Genomes)")
plt.xlabel("Genus")
plt.ylabel("Count")
plt.xticks(rotation=90, ha="right", fontsize=7)  # Rotate labels by 90 degrees and set label font size
plt.tight_layout()

# Save the plot as a PDF file
plt.savefig("Summary_feed_filtered.pdf", format="pdf")

print("Plot saved as PDF:", "Summary_feed_filtered.pdf")

