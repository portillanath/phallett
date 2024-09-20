import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from pathlib import Path

def alignment_skani(file_path, output_dir, subdir_name):
    # Read the DataFrame from the file
    df = pd.read_csv(file_path, sep='\t')  # Adjust the separator if necessary

    # Display the first few rows of the DataFrame to check the column names
    print(f"Processing file: {file_path}")
    print(df.head())
    print("Available columns:", df.columns.tolist())  # Print available columns for debugging

    # Check if the necessary columns exist
    required_columns = ['Align_fraction_query', 'ANI']
    missing_columns = [col for col in required_columns if col not in df.columns]

    if not missing_columns:
        # Rename columns
        df.rename(columns={'Ref_file': 'GenomaA', 'Query_file': 'GenomeB'}, inplace=True)
        
        # Extract values before the last '/' and before the '.' for GenomaA and GenomeB
        df['GenomaA'] = df['GenomaA'].apply(lambda x: os.path.basename(x).split('.')[0])
        df['GenomeB'] = df['GenomeB'].apply(lambda x: os.path.basename(x).split('.')[0])

        # Create a figure and axis object
        plt.figure(figsize=(10, 6))
        
        # Plot the scatter plot
        sns.scatterplot(data=df, x='Align_fraction_query', y='ANI', alpha=0.7)

        # Fit the linear regression model
        X = df[['Align_fraction_query']].values.reshape(-1, 1)
        y = df['ANI'].values
        model = LinearRegression()
        model.fit(X, y)
        y_pred = model.predict(X)

        # Plot the regression line
        plt.plot(df['Align_fraction_query'], y_pred, color='red', linewidth=2, label='Linear Regression')

        # Update title with subdirectory name
        plt.title(f'Scatter Plot of Align_fraction_query vs ANI with Linear Regression Skani - {subdir_name}')
        plt.xlabel('Align_fraction_query')
        plt.ylabel('ANI')
        plt.grid(True)
        plt.legend()

        # Save the plot as a PDF file
        plot_pdf_path = os.path.join(output_dir, f'alignment_fraction_{subdir_name}.pdf')
        plt.savefig(plot_pdf_path, format='pdf')
        print(f'Saved the scatter plot as a PDF to {plot_pdf_path}')
        
        # Close the plot to avoid display and free up memory
        plt.close()

        # Create DataFrame with the linear regression values
        df['Linear_Regression'] = y_pred

        # Add the algorithm column
        df['algorithm'] = 'skani'
        df['kmer'] = 'static'

        # Save the new DataFrame as a CSV file
        alignment_fraction_df = df[['GenomaA', 'GenomeB', 'ANI', 'Align_fraction_query', 'Linear_Regression', 'algorithm', 'kmer']]
        output_csv_path = os.path.join(output_dir, 'alignment_fraction.csv')
        alignment_fraction_df.to_csv(output_csv_path, index=False)
        print(f'Saved the alignment fraction DataFrame to {output_csv_path}')
    else:
        print(f"Required columns are missing. Missing columns: {missing_columns}")

def alignment_fastani(directory):
    dataframes = {}
    
    # Loop through files in the directory
    for filename in os.listdir(directory):
        # Check if filename matches the new expected pattern
        if filename.startswith('fastani_'):
            # Adjust the regex pattern based on the filename structure
            match = re.match(r'^fastani_(.*?)_frag_500_(\d+)', filename)
            if match:
                genus_subdir = match.group(1)
                kmer_size = int(match.group(2))
                
                # Define the path to the file
                file_path = os.path.join(directory, filename)
                
                # Read the file into a DataFrame with manual column names
                try:
                    # Load the DataFrame without headers, assign headers manually
                    df = pd.read_csv(file_path, sep='\t', header=None, names=['GenomaA', 'GenomeB', 'ANI', 'mappings', 'total_fragments'])
                    print(f"Read {filename} into DataFrame with shape {df.shape}")
                except Exception as e:
                    print(f"Error reading {filename}: {e}")
                    continue

                # Calculate AF and add as a new column
                try:
                    df['AF'] = (df['mappings'] / df['total_fragments']) * 100
                    dataframes[kmer_size] = df
                except Exception as e:
                    print(f"Error calculating AF for {filename}: {e}")
            else:
                print(f"Filename does not match expected pattern: {filename}")

    if not dataframes:
        print("No valid fastani files found in the directory.")
        return

    # Plot scatter plots and perform linear regression for each k-mer size
    results = []
    num_plots = len(dataframes)
    fig, axes = plt.subplots(nrows=num_plots, ncols=1, figsize=(12, 6 * num_plots))
    if num_plots == 1:
        axes = [axes]  # Ensure axes is iterable if only one plot

    for idx, (kmer_size, df) in enumerate(dataframes.items()):
        ax = axes[idx]
        ax.scatter(df['ANI'], df['AF'], alpha=0.7)
        ax.set_title(f'Scatter Plot of ANI vs AF for k-mer size {kmer_size} - {genus_subdir}')
        ax.set_xlabel('ANI')
        ax.set_ylabel('AF')
        ax.grid(True)
        
        # Fit the linear regression model
        X = df[['ANI']].values.reshape(-1, 1)
        y = df['AF'].values
        model = LinearRegression()
        model.fit(X, y)
        y_pred = model.predict(X)

        # Add linear regression values to the DataFrame
        df['Linear_Regression'] = y_pred
        df['algorithm'] = 'fastani'
        df['kmer'] = kmer_size

        # Collect results
        results.append(df[['GenomaA', 'GenomeB', 'ANI', 'AF', 'Linear_Regression', 'algorithm', 'kmer']])
        
        # Plot the regression line
        ax.plot(df['ANI'], y_pred, color='red', linewidth=2, label='Linear Regression')
        ax.legend()

    # Save the grid of plots as a PDF
    plot_pdf_path = os.path.join(directory, 'ANI_vs_AF_grid.pdf')
    plt.tight_layout()
    plt.savefig(plot_pdf_path, format='pdf')
    plt.close()
    print(f'Saved the grid of scatter plots as a PDF to {plot_pdf_path}')

    # Concatenate all DataFrames and save to CSV
    if results:
        all_results = pd.concat(results)
        output_csv_path = os.path.join(directory, 'alignment_fraction_fastani.csv')
        all_results.to_csv(output_csv_path, index=False)
        print(f'Saved the alignment fraction DataFrame to {output_csv_path}')
    else:
        print("No results to save for fastani.")

def main():
    # Set paths
    current_dir = Path.cwd()
    parent_dir = current_dir.parent
    workdir = parent_dir / "test" / "Metrics_Results"

    # Check if the directory exists
    print(f"Work directory: {workdir}")
    print("Directory exists:", workdir.exists())

    # List subdirectories in the source directory
    subdirectories = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d))]

    # Process files in each subdirectory
    for genus in subdirectories:
        genus_dir = os.path.join(workdir, genus)

        # Ensure the output directory exists
        output_dir = genus_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Process skani files
        skani_file_path = os.path.join(genus_dir, f'skani_distance_{genus}.txt')
        if os.path.isfile(skani_file_path):
            alignment_skani(skani_file_path, output_dir, genus)
        else:
            print(f"Skani file {skani_file_path} does not exist.")

        # Process fastani files
        alignment_fastani(genus_dir)

if __name__ == "__main__":
    main()
