import argparse
import pandas as pd
from pathlib import Path

def extract_data_from_folders(parent_dir, folder_names):
    # Define the paths
    data_taxa_selected = Path(parent_dir) / "phallett" / "data" / "Taxa_Selected"
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

    # Initialize a list to store the data
    result_data = []

    # Process each folder
    for folder_name in folder_names:
        folder_path = data_taxa_selected / folder_name
        if not folder_path.is_dir():
            print(f"Warning: The folder '{folder_path}' does not exist.")
            continue
        
        # Iterate through files in the folder
        for file in folder_path.iterdir():
            if file.is_file():
                # Extract the Virus_GENBANK_accession from the file name
                virus_accession = file.stem  # Assuming the file name matches the accession number
                # Find the corresponding row in the VMR DataFrame
                matching_rows = vmr_df[vmr_df['Virus_GENBANK_accession'] == virus_accession]
                
                # If there are matching rows, add to result_data
                if not matching_rows.empty:
                    result_data.append(matching_rows)

    # Concatenate all collected data into a single DataFrame
    if result_data:
        final_df = pd.concat(result_data, ignore_index=True)
        # Save to a new CSV
        output_file = Path(parent_dir) / "phallett" / "data" / "Virus_Metadata_Resource" / "output.csv"
        final_df.to_csv(output_file, index=False)
        print(f"Data extracted and saved to {output_file}")
    else:
        print("No matching data found.")

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
