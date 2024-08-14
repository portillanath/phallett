import argparse
import pandas as pd
from pathlib import Path

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
                # Find the corresponding row in the VMR DataFrame
                matching_rows = vmr_df[vmr_df['Virus GENBANK accession'] == virus_accession]
                
                # If there are matching rows, add to result_data
                if not matching_rows.empty:
                    result_data.append(matching_rows)

        # Concatenate all collected data into a single DataFrame for the current folder
        if result_data:
            final_df = pd.concat(result_data, ignore_index=True)
            # Define the output file path for the current folder

            # Define the output folder path
            output_directory = Path(general_path) /"test"/"Metrics_Results" / folder_name
            output_directory.mkdir(parents=True, exist_ok=True)  # Create the directory if it doesn't exist

            # Define the output file path for the current folder
            output_file = output_directory / f"metadata_{folder_name}.csv"
            final_df.to_csv(output_file, index=False)
            print(f"Data extracted and saved to {output_file}")
        else:
            print(f"No matching data found for folder '{folder_name}'.")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process folder names and extract data from files.')
    parser.add_argument('parent_directory', type=str, help='Path to the parent directory')
    parser.add_argument('folder_names', nargs='+', help='Names of folders to process')

    # Parse the arguments
    args = parser.parse_args()
    
    # Call the extraction function
    extract_data_from_folders(args.parent_directory, args.folder_names)

if __name__ == "__main__":
    main()
