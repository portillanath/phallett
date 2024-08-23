import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from itertools import product
import os as os 


def calcualte_genome_size(fasta_file):
    genome_size=0
    for record in SeqIO.parse(fasta_file,"fasta"):
        genome_size+=len(record.seq)
    return genome_size

def extract_data_from_folders(parent_dir, folder_names):
    # Define the paths
    data_taxa_selected = Path(parent_dir)
    data_path = data_taxa_selected.parent 
    general_path=data_path.parent
    vmr_file_path = Path(data_path) / "Virus_Metadata_Resource" / "VMR.csv"
    
    # Read the VMR CSV into a DataFrame
    try:
        vmr_df = pd.read_csv(vmr_file_path)
    except FileNotFoundError:
        print(f"Error: VMR.csv not found at {vmr_file_path}")
        return
    except pd.errors.EmptyDataError:
        print("Error: VMR.csv is empty.")
        return

    # Process each folder
    for folder_name in folder_names:
        folder_path = data_taxa_selected / folder_name
        if not folder_path.is_dir():
            print(f"Warning: The folder '{folder_path}' does not exist.")
            continue
        
        # Initialize a list to store the data for the current folder
        result_data = []

        # Iterate through files in the folder
        for file in folder_path.iterdir():
            if file.is_file():
                # Extract the Virus_GENBANK_accession from the file name
                virus_accession = file.stem  # Assuming the file name matches the accession number
                if file.suffix == '.fasta':
                   genome_size= calcualte_genome_size(file)
                   genome_size_info={'Virus GENBANK accession': virus_accession, 'Genome Size':genome_size}

                   #Append to the metada 
                   if virus_accession in vmr_df['Virus GENBANK accession'].values:
                      metadata_row= vmr_df[vmr_df['Virus GENBANK accession']==virus_accession].iloc[0].to_dict()
                      metadata_row.update(genome_size_info)
                      result_data.append(metadata_row)
                else:
                   if virus_accession in vmr_df['Virus GENBANK accession'].values:
                      metadata_row=vmr_df[vmr_df['Virus GENBANKS accession']==virus_accession].iloc[0].to_dict()
                      result_data.append(metadata_row)

        # Concatenate all collected data into a single DataFrame for the current folder
        if result_data:
            final_df = pd.DataFrame(result_data)
            # Define the output file path for the current folder

            # Define the output folder path
            output_directory = Path(general_path) /"test"/"Metrics_Results" / folder_name
            output_directory.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist

            # Define the output file path for the current folder
            output_file = output_directory / f"metadata_{folder_name}.csv"
            final_df.to_csv(output_file, index=False)
        else:
            print(f"No matching data found for folder '{folder_name}'.")

def create_bubble_data(parent_dir, folder_names):
   data_taxa_selected = Path(parent_dir)
   data_path = data_taxa_selected.parent 
   general_path=data_path.parent
   test_path=Path(general_path) / "test" / "Metrics_Results"

   for folder_name in folder_names:
        folder_path = test_path / folder_name
        filemetadata = f"metadata_{folder_name}.csv"
        csvmetadata = Path(folder_path)/filemetadata
        metadata_df=pd.read_csv(csvmetadata, low_memory=False)
        filesummary= f"summary_{folder_name}.csv"
        csvsummary=Path(folder_path)/filesummary
        summary_df=pd.read_csv(csvsummary, low_memory=False)
        merged_df = pd.merge(summary_df, metadata_df, left_on='GenomeA', right_on='Virus GENBANK accession', how='left')
        merged_df.drop(columns='Virus GENBANK accession', inplace=True)
        bubble_output_file = Path(folder_path) / f"bubble_{folder_name}.csv"
        merged_df.to_csv(bubble_output_file, index=False)

def plot_bubble(parent_dir, folder_names):
    data_taxa_selected = Path(parent_dir)
    data_path = data_taxa_selected.parent
    general_path = data_path.parent
    test_path = general_path / "test" / "Metrics_Results"
    
    for folder_name in folder_names:
        folder_path = test_path / folder_name
        folder_bubbles = f"bubble_{folder_name}.csv"
        csvbubble = folder_path / folder_bubbles
        bubble_df = pd.read_csv(csvbubble, low_memory=False)
        unique_ani=bubble_df['algorithm_ani'].unique()
        unique_mash= bubble_df['algorithm_mash'].unique()
        combinations=product(unique_ani,unique_mash)
        for ani,mash in combinations:
            subset=bubble_df[(bubble_df['algorithm_ani']==ani)&(bubble_df['algorithm_mash']==mash )]
            pdf_filename = f'bubble_{folder_name}_{ani}_{mash}.pdf'
            output_path = os.path.join(folder_path, pdf_filename)
            with PdfPages(output_path)as pdf:
                kmer_pairs=product(subset['kmer_ani'].unique(),subset['kmer_mash'].unique())
                for kmer_ani, kmer_mash in kmer_pairs:
                    plot_df=subset[(subset['kmer_ani']==kmer_ani)&(subset['kmer_mash']==kmer_mash)]
                    if not plot_df.empty:
                        plt.figure(figsize=(5,4))
                        sns.scatterplot(data=plot_df,
                                        x='ani_distance',
                                        y='mash_distance',
                                        hue='Species',
                                        size='Genome Size',
                                        sizes=(20,200),
                                        legend='full',
                                        palette='husl')
                        plt.title(f'Bubble kmer_ani={kmer_ani} kmer_mash={kmer_mash}')
                        plt.xlabel('ANI')
                        plt.ylabel('Mash disntance') 
                        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Species', title_fontsize='3', fontsize='2')                 
                        pdf.savefig()
                        plt.close()
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process folder names and extract data from files.')
    parser.add_argument('parent_directory', type=str, help='Path to the parent directory')
    parser.add_argument('folder_names', nargs='+', help='Names of folders to process')

    # Parse the arguments
    args = parser.parse_args()
    
    # Call the extraction function
    extract_data_from_folders(args.parent_directory, args.folder_names)
    create_bubble_data(args.parent_directory,args.folder_names)
    plot_bubble(args.parent_directory, args.folder_names)

    # Call the mergin function 

if __name__ == "__main__":
    main()