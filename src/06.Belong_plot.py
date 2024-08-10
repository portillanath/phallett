import os
import pandas as pd

def find_genbank_accessions(parent_dir, csv_path):
    # Read the CSV file
    ictv_data = pd.read_csv(csv_path)

    # Get the relevant columns from the DataFrame
    relevant_columns = ['Virus GENBANK accession', 'Order', 'Suborder', 
                        'Family', 'Subfamily', 'Genus', 'Subgenus', 'Species']
    
    # Create a new DataFrame to hold the results
    results = []

    # Loop through each file in the specified directory
    for root, dirs, files in os.walk(parent_dir):
        for file in files:
            # Check if the file name corresponds to any Virus GENBANK accession
            if file in ictv_data['Virus GENBANK accession'].values:
                # Get the row corresponding to the found accession
                row = ictv_data[ictv_data['Virus GENBANK accession'] == file]
                
                # Append the required data to the results list
                results.append(row[relevant_columns].values.flatten().tolist())

    # Create a DataFrame from the results and save to a CSV
    results_df = pd.DataFrame(results, columns=relevant_columns)
    
    # Specify the output file path
    output_file_path = os.path.join(parent_dir, 'genbank_results.csv')
    
    # Save the DataFrame to a CSV file
    results_df.to_csv(output_file_path, index=False)

    print(f"Results saved to {output_file_path}")

# Example usage
parent_dir = '/path/to/your/directory'
csv_path = os.path.join(parent_dir, 'phallett/data/Virus_Metadata_Resource/VMR.csv')
find_genbank_accessions(parent_dir, csv_path)
