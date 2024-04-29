# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 17:37:24 2023

@author: rkerrid
"""
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import json

from dia_sis.workflow.utils import manage_directories



def create_report(df, path, params):
    # Construct the description string
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_loose_filters'].items()])
    description = f"Parameters used:\n{params_str}"
    
    counts_df = df['Run'].value_counts()
    
    # Set up the PDF
    # create_reports_directory(path)
    manage_directories.create_directory(path,'reports')
    output_dir = path + '/reports'
    pdf_path = os.path.join(output_dir, 'silac_precursors_report.pdf')
    
    # # Construct the description string
    # # Load filtering parameters from JSON
    # CONFIG_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'configs')
    # json_path = os.path.join(CONFIG_DIR, 'filtering_parameters.json')
    # with open(json_path, 'r') as f:
    #     params = json.load(f)
    # params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_strict_filters'].items()])
    # description = f"Parameters used:\n{params_str}"
    
    df_grouped = df.groupby('Run')
    
    
    
    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(8, 11))
        plt.axis('off')
        plt.text(0.5, 0.98, "Silac Precursors QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)

        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        

                

        # Display table of counts
        table_data = [["Run", "DF Count"]]
        for run in counts_df.keys():
            table_data.append([run, counts_df.get(run, 0)])
        
        plt.table(cellText=table_data, cellLoc = 'center', loc='center', colWidths=[0.2,0.2])
        plt.axis('off')
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

        # Get the list of unique runs
        runs = df_grouped.groups.keys()

        # For each run, plot the histograms and save to PDF
        for run in runs:
            plot_histograms_for_run(run, df_grouped, ['Precursor.Translated H','Precursor.Translated L'])
            
            # )
            pdf.savefig()  # Saves the current figure into the PDF
            plt.close()
             
def plot_histograms_for_run(run, df_grouped, labels):
    """Plot histograms for a given run from multiple dataframes."""
    colors = sns.color_palette("husl", 3)  # Get a color palette

    # Extract the DataFrame for the specific run
    df_run = df_grouped.get_group(run)

    for label, color in zip(labels, colors):
        # Calculate log2 values
        data = np.log2(df_run[label].dropna())  # Add 1 to avoid log2(0)
        
        # Filter out -inf and inf values
        data = data[np.isfinite(data)]
        plt.hist(data, alpha=0.5, label=label, bins=300, color=color)

    plt.title(f'Histogram for {run}')
    plt.xlabel("log2 Intensity")
    plt.ylabel('Frequency')
    plt.legend()

