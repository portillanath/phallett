import os
import glob
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

# Set working directory
current_dir = Path.cwd()
parent_dir = current_dir.parent
workdir = os.path.expanduser((Path(parent_dir) / "phallett" / "test" / "Metrics_Results"))
subdirectories = [os.path.join(workdir, name) for name in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, name))]

# Create the argument parser
parser = argparse.ArgumentParser(description='Process some arguments.')

# Add arguments
parser.add_argument('-mx', type=str, help='The mx argument')
parser.add_argument('-kmersx', type=str, help='The kmersx argument')
parser.add_argument('-my', type=str, help='The my argument')
parser.add_argument('-kmersy', type=str, help='The kmersy argument')

# Parse the arguments
args = parser.parse_args()

# Access the arguments
print(f'mx: {args.mx}, kmersx: {args.kmersx}, my: {args.my}, kmersy: {args.kmersy}')

mx = args.mx
my = args.my
kmersx = [int(k) for k in args.kmersx.split(",")]
kmersy = [int(k) for k in args.kmersy.split(",")]

# Define tools based on arguments
tool_my = []
tool_mx = []

if my == "mash":
    tool_my = ["mash", "sourmash"]
elif my == "ani":
    tool_my = ["fastani", "skani"]
elif my == "aai":
    tool_my = ["comparem"]
elif my == "viridic":
    tool_my = ["viridic"]
elif my == "vcontact2":
    tool_my = ["vcontact2"]

if mx == "mash":
    tool_mx = ["mash", "sourmash"]
elif mx == "ani":
    tool_mx = ["fastani", "skani"]
elif mx == "aai":
    tool_mx = ["comparem"]
elif mx == "viridic":
    tool_mx = ["viridic"]
elif mx == "vcontact2":
    tool_mx = ["vcontact2"]

# Initialize DataFrames
mx_data = pd.DataFrame()
my_data = pd.DataFrame()

for genus in subdirectories:
    os.chdir(genus)
    genus_name = os.path.basename(genus)
    genus_dir = os.path.join(workdir, genus_name)
    summaries = glob.glob(os.path.join(genus_dir, "*metrics*.csv"))

    # Load mx data
    if mx in [mx, my]:
        mx_files = [s for s in summaries if f"{mx}_metrics_" in s]
        if mx_files:
            mx_data = pd.read_csv(mx_files[0])
            print(f"mx_data for {genus_name} loaded:\n{mx_data.head()}")
        else:
            print(f"No {mx} data for {genus_name}")

    # Load my data
    if my in [mx, my]:
        my_files = [s for s in summaries if f"{my}_metrics_" in s]
        if my_files:
            my_data = pd.read_csv(my_files[0])
            print(f"{genus_name} {my} data loaded:\n{my_data.head()}")
        else:
            print(f"No {my} data for {genus_name}")

    for algorithm_mx in tool_mx:
        for algorithm_my in tool_my:
            df_subset_mx = mx_data[mx_data[f"algorithm_{mx}"] == algorithm_mx]
            df_subset_my = my_data[my_data[f"algorithm_{my}"] == algorithm_my]

            # Check for empty DataFrames
            if df_subset_mx.empty or df_subset_my.empty:
                print(f"Empty DataFrames for mx: {algorithm_mx} or my: {algorithm_my} in genus: {genus_name}.")
                continue

            kmer_values_mx = sorted(pd.unique(df_subset_mx[f"kmer_{mx}"]))
            kmer_values_my = sorted(pd.unique(df_subset_my[f"kmer_{my}"]))

            if not kmer_values_my or not kmer_values_mx:
                print(f"No kmer values for {my} with algorithm {algorithm_my} in genus {genus_name}")
                continue

            fig, axes = plt.subplots(len(kmer_values_mx), len(kmer_values_my), figsize=(13, 9), sharex=True)
            axes = axes.reshape(len(kmer_values_mx), len(kmer_values_my))
            fig.text(0.04, 0.04, algorithm_my, va='center', rotation='vertical', fontsize=9)

            for j, my_kmer in enumerate(kmer_values_my):
                for i, mx_kmer in enumerate(kmer_values_mx):
                    df_subset_mx_kmer = df_subset_mx[df_subset_mx[f"kmer_{mx}"] == mx_kmer]
                    df_subset_my_kmer = df_subset_my[df_subset_my[f"kmer_{my}"] == my_kmer]

                    merged_df = pd.merge(df_subset_mx_kmer, df_subset_my_kmer, on=["GenomeA", "GenomeB"], how="outer")
                    merged_df.sort_values(by=[f"{mx}_distance", f"{my}_distance"], inplace=True)
                    merged_df = merged_df.dropna()

                    if not merged_df.empty:
                        colors = []
                        for x_value in merged_df[f"{mx}_distance"]:
                            if x_value <= 88.00:
                                colors.append('red')
                            elif 88.00 < x_value <= 94.00:
                                colors.append('orange')
                            else:
                                colors.append('pink')

                        ax = axes[i, j]
                        ax.set_title(f"{algorithm_mx} kmer: {kmer_values_mx[i]}", fontsize=9)
                        ax.tick_params(axis='both', which='major', labelsize=9)

                        scatter = ax.scatter(merged_df[f"{mx}_distance"], merged_df[f"{my}_distance"],
                                             c=colors, alpha=0.5, s=1)

                        ax.set_ylabel(f"{algorithm_my} kmer: {kmer_values_my[j]}", rotation=90, ha='center', fontsize=10)
                        ax.yaxis.set_label_coords(-0.25, 0.5)

            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.4, hspace=0.4)
            handles = [
                plt.Line2D([0], [0], marker='o', color='w', label='Family >=88.00', markerfacecolor='red', markersize=7),
                plt.Line2D([0], [0], marker='o', color='w', label='Genus 88.00-94.00', markerfacecolor='orange', markersize=7),
                plt.Line2D([0], [0], marker='o', color='w', label='Species >=94.00', markerfacecolor='pink', markersize=7)
            ]
            fig.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.0, 1.0), fontsize=9, frameon=False)

            pdf_filename = f"{genus_name}_{algorithm_mx}_{algorithm_my}.pdf"
            print(f"Saving PDF: {pdf_filename}")
            with PdfPages(pdf_filename) as pdf:
                pdf.savefig(fig)
                plt.close(fig)
