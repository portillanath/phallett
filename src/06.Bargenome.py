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
    data_taxa_selected = Path(parent_dir) / "data" / "Taxa_Selected"
    general_path = Path(parent_dir) / "test" / "Metrics_Results"
    vmr_file_path = Path(parent_dir) / "data" / "Virus_Metadata_Resource" / "VMR.csv"
    
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
                
                # Define the output file path for the current folder
                output_directory = general_path / folder_name.name
                output_directory.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist

                # Define the output file path for the current folder
                output_file = output_directory / f"metadata_{folder_name.name}.csv"
                df.to_csv(output_file, index=False)
                
                # Generate and save the bar plot
                plot_barplot()
            else:
                print(f"No matching data found for folder '{folder_name.name}'.")

def plot_barplot():
    plt.figure(figsize=(14, 12))  # Increased height to provide more space for the legend
    
    # Ensure 'Genus' column exists in DataFrame
    if 'Genus' not in df.columns:
        print("Error: 'Genus' column not found in DataFrame.")
        return
    
    # Create bar plot with hue for genera
    barplot = sns.barplot(data=df, x='Virus GENBANK accession', y='Genome Size', hue='Genus', palette='Set2')
    plt.xticks(rotation=90)
    plt.title(f'Genome Sizes by Virus GENBANK Accession for {folder_name.name}')
    plt.xlabel('Virus GENBANK Accession')
    plt.ylabel('Genome Size')
    
    # Adjust legend position to be below the plot
    plt.legend(title='Genus', bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol=5, frameon=False)
    
    # Adjust layout to prevent clipping and ensure plot fits well
    plt.tight_layout()
    
    # Save plot as PDF with folder name
    plot_pdf_file = output_directory / f'{folder_name.name}_genome_sizes_plot.pdf'
    plt.savefig(plot_pdf_file)
    plt.close()

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
