from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram


# Set working directory
current_dir = Path.cwd()
parent_dir = current_dir.parent
workdir = parent_dir / "test" / "Metrics_Results"

# List subdirectories in the work directory
subdirectories = [d for d in workdir.iterdir() if d.is_dir()]

# Process each subdirectory
for genus in subdirectories:
    genus_name = genus.name
    aai_data_dir = genus / f"{genus_name}_compare_results" / "aai"
    aai_tsv = aai_data_dir / "aai_summary.tsv"

    try:
        # Read the TSV file into a DataFrame
        aai_df = pd.read_csv(aai_tsv, sep='\t')
  
        # Create a heatmap
        heatmap_data = aai_df.pivot(index='#Genome A', columns='Genome B', values='Mean AAI')
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, cmap='flare', cbar_kws={'label': 'Mean AAI'})
        plt.title('Heatmap of Average Amino Acid Identity (AAI)')
        plt.xlabel('Genome B')
        plt.ylabel('Genome A')
        heatmap_file = genus / f"{genus_name}_AAI_heatmap.pdf"
        plt.savefig(heatmap_file, format='pdf')
        plt.close()  # Close the plot to free memory

        print(f"Heatmap saved as {heatmap_file}")

    except FileNotFoundError:
        print(f"File not found: {aai_tsv}. Skipping to the next directory.")
        continue
    except Exception as e:
        print(f"An error occurred while processing {aai_tsv}: {e}")
        continue

# Process each subdirectory
for genus in subdirectories:
    genus_name = genus.name
    aai_data_dir = genus / f"{genus_name}_compare_results" / "aai"
    aai_tsv = aai_data_dir / "aai_summary.tsv"

    try:
        # Read the TSV file into a DataFrame
        aai_df = pd.read_csv(aai_tsv, sep='\t')
  
        # Create a heatmap
        heatmap_data = aai_df.pivot(index='#Genome A', columns='Genome B', values='Orthologous fraction (OF)')
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, cmap='viridis', cbar_kws={'label': 'Orthologous fraction (OF)'})
        plt.title('Heatmap of Orthologous fraction (OF)')
        plt.xlabel('Genome B')
        plt.ylabel('Genome A')
        heatmap_file = genus / f"{genus_name}_Orthologous fraction (OF)_heatmap.pdf"
        plt.savefig(heatmap_file, format='pdf')
        plt.close()  # Close the plot to free memory

        print(f"Heatmap saved as {heatmap_file}")

    except FileNotFoundError:
        print(f"File not found: {aai_tsv}. Skipping to the next directory.")
        continue
    except Exception as e:
        print(f"An error occurred while processing {aai_tsv}: {e}")
        continue



