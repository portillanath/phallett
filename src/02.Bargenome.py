import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

# Global variables for the plot function
df = None
output_directory = None
folder_name = None

def calculate_genome_size(fasta_file):
    genome_size = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_size += len(record.seq)
    return genome_size

def extract_data_from_folders(parent_dir):
    global df, output_directory, folder_name
    
    # Define paths
    data_taxa_selected = Path(parent_dir) / "phallett" / "data" / "Taxa_Selected"
    general_path = Path(parent_dir) / "phallett" / "test" / "Metrics_Results"
    vmr_file_path = Path(parent_dir) / "phallett" / "data" / "Virus_Metadata_Resource" / "VMR.csv"
    
    # Read the VMR CSV into a DataFrame
    try:
        vmr_df = pd.read_csv(vmr_file_path)
    except FileNotFoundError:
        print(f"Error: VMR.csv not found at {vmr_file_path}")
        return
    except pd.errors.EmptyDataError:
        print("Error: VMR.csv is empty.")
        return

    # Process each folder in the Taxa_Selected directory
    for folder_name in data_taxa_selected.iterdir():
        if folder_name.is_dir():
            # Initialize a list to store the data for the current folder
            result_data = []

            # Iterate through files in the folder
            for file in folder_name.iterdir():
                if file.is_file():
                    # Extract the Virus_GENBANK_accession from the file name
                    virus_accession = file.stem  # Assuming the file name matches the accession number
                    if file.suffix == '.fasta':
                        genome_size = calculate_genome_size(file)
                        genome_size_info = {'Virus GENBANK accession': virus_accession, 'Genome Size': genome_size}

                        # Append to the metadata 
                        if virus_accession in vmr_df['Virus GENBANK accession'].values:
                            metadata_row = vmr_df[vmr_df['Virus GENBANK accession'] == virus_accession].iloc[0].to_dict()
                            metadata_row.update(genome_size_info)
                            result_data.append(metadata_row)
                    else:
                        if virus_accession in vmr_df['Virus GENBANK accession'].values:
                            metadata_row = vmr_df[vmr_df['Virus GENBANK accession'] == virus_accession].iloc[0].to_dict()
                            result_data.append(metadata_row)

            # Concatenate all collected data into a single DataFrame for the current folder
            if result_data:
                df = pd.DataFrame(result_data)
                print(f"Data for {folder_name.name}:")
                print(df.head())  # Print the first few rows of the DataFrame for inspection

                # Check if the necessary columns are present
                if 'Genus' not in df.columns:
                    print(f"Error: 'Genus' column missing for folder {folder_name.name}.")
                    continue  # Skip this folder and move to the next one

                # Define the output file path for the current folder
                output_directory = general_path / folder_name.name
                output_directory.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist

                # Define the output file path for the current folder
                output_file = output_directory / f"metadata_{folder_name.name}.csv"
                df.to_csv(output_file, index=False)
                
                # Generate and save the bar plot
                plot_barplot(folder_name)
            else:
                print(f"No matching data found for folder '{folder_name.name}'.")

def plot_barplot(folder_name):
    global df, output_directory
    
    # Ensure 'Genus' column exists in DataFrame
    if 'Genus' not in df.columns:
        print(f"Error: 'Genus' column not found in DataFrame for {folder_name.name}. Skipping plot.")
        return
    
    print(f"Generating plot for {folder_name.name}...")
    
    try:
        # Create bar plot with hue for genera
        plt.figure(figsize=(14, 12))  # Increased height to provide more space for the legend
        barplot = sns.barplot(data=df, x='Virus GENBANK accession', y='Genome Size', hue='Genus', palette='Set2')
        plt.xticks(rotation=90)
        plt.title(f'Genome Sizes by Virus GENBANK Accession for {folder_name.name}')
        plt.xlabel('Virus GENBANK Accession')
        plt.ylabel('Genome Size')

        # Get the top 10 genera based on their count in the DataFrame
        top_10_genera = df['Genus'].value_counts().head(10).index.tolist()

        # Filter the legend to show only top 10 genera
        handles, labels = barplot.get_legend_handles_labels()
        filtered_handles = []
        filtered_labels = []

        for handle, label in zip(handles, labels):
            if label in top_10_genera:
                filtered_handles.append(handle)
                filtered_labels.append(label)

        # Set the legend to only display the top 10 genera and place it inside the upper-right corner
        plt.legend(handles=filtered_handles, labels=filtered_labels, title='Genus', loc='upper right', frameon=False, ncol=1)

        # Adjust layout to prevent clipping and ensure plot fits well
        plt.tight_layout()

        # Save plot as PDF with folder name
        plot_pdf_file = output_directory / f'{folder_name.name}_genome_sizes_plot.pdf'
        print(f"Saving plot to {plot_pdf_file}")  # Debugging: print output path
        plt.savefig(plot_pdf_file)
        plt.close()
    except Exception as e:
        print(f"Error generating plot for {folder_name.name}: {e}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process folders and extract data from files.')
    parser.add_argument('parent_directory', type=str, help='Path to the parent directory')

    # Parse the arguments
    args = parser.parse_args()
    
    # Call the extraction function
    extract_data_from_folders(args.parent_directory)

if __name__ == "__main__":
    main()
