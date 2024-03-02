# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:48:12 2023

@author: rkerrid
"""

import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from dia_sis.pipeline.utils import manage_directories


def create_report(df, contams, filtered_out, path, params):
    # Count rows for each dataframe
    counts_df = df['Run'].value_counts()
    counts_contams = contams['Run'].value_counts()
    counts_filtered_out = filtered_out['Run'].value_counts()

    # Construct the description string
    strict_params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_strict_filters'].items()])
    loose_params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_loose_filters'].items()])

    description = f"Strict parameters used:\n{strict_params_str} \n\n Loose parameters used:\n{loose_params_str}"
    
    #Group by run for analysis of each
    df_grouped = df.groupby('Run')
    contams_grouped = contams.groupby('Run')
    filtered_out_grouped = filtered_out.groupby('Run')

    # Set up the PDF
    manage_directories.create_directory(path,'reports')
    output_dir = path + '/reports'
    pdf_path = os.path.join(output_dir, 'filtering_report.pdf')

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(8, 11))
        plt.axis('off')
        plt.text(0.5, 0.98, "Filtering QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.75, description, ha='center', va='center', wrap=True)
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        
        # Display table of counts
        table_data = [["Run", "DF Count", "Contaminants Count", "Filtered Out Count"]]
        for run in counts_df.keys():
            table_data.append([run, counts_df.get(run, 0), counts_contams.get(run, 0), counts_filtered_out.get(run, 0)])
        
        plt.table(cellText=table_data, cellLoc = 'center', loc='center', colWidths=[0.2,0.2,0.3,0.3])
        plt.axis('off')
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

        # Get the list of unique runs
        runs = df_grouped.groups.keys()

        # For each run, plot the histograms and save to PDF
        # Suppress runtime warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning) 
        for run in runs:
            plot_histograms_for_run(run, [df_grouped, contams_grouped, filtered_out_grouped], ['df', 'contams', 'filtered_out'])
            # )
            pdf.savefig()  # Saves the current figure into the PDF
            plt.close()


def plot_histograms_for_run(run, dfs, labels, column='Precursor.Quantity'):
    """Plot histograms for a given run from multiple dataframes."""
    colors = sns.color_palette("husl", len(dfs))  # Get a color palette

    for df, label, color in zip(dfs, labels, colors):
        # Calculate log2 values
        data = np.log2(df.get_group(run)[column].dropna())
        # Filter out -inf and inf values
        data = data[np.isfinite(data)]
        plt.hist(data, alpha=0.5, label=label, bins=300, color=color)  # Removed edgecolor and added color parameter

    plt.title(f'Histogram for {run}')
    plt.xlabel(f"log2({column})")
    plt.ylabel('Frequency')
    plt.legend()


