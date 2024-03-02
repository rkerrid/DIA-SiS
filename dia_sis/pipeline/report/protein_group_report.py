# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 13:33:49 2023

@author: rkerrid
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import json 


def create_protein_groups_report(df, path):
    # Count rows for each dataframe
    counts_df = df['Run'].value_counts().reset_index()
    counts_df.columns = ['Run', 'Counts']  # Naming the columns

# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import seaborn as sns
# import os
# from silac_dia_tools.pipeline.utils import manage_directories


# def create_report(df, path, params):
#     # Construct the description string
#     params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
#     description = f"Parameters used:\n{params_str}"
    
#     # Count rows for each dataframe
#     counts_df = df['Run'].value_counts().reset_index()
#     counts_df.columns = ['Run', 'Counts']  # Naming the columns
    
#     # Set up the PDF
#     # create_reports_directory(path)
#     manage_directories.create_directory(path,'reports')
#     output_dir = path + '/reports'
#     pdf_path = os.path.join(output_dir, 'protein_groups_report.pdf')

#     with PdfPages(pdf_path) as pdf:
#         # Title and introduction
#         plt.figure(figsize=(11, 8))
#         plt.axis('off')
#         plt.text(0.5, 0.98, "Protein Ratios QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
#         plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)
        
#         pdf.savefig()  # Saves the current figure into the PDF
#         plt.close()
        
#         # Creating the barplot
#         plt.figure(figsize=(10, 7))
#         barplot = sns.barplot(x='Run', y='Counts', data=counts_df, orient='v')
#         plt.title('Counts per Run')
#         plt.xlabel('Run')
#         plt.ylabel('Counts')
        
#         # Rotating x-axis labels
#         barplot.set_xticklabels(barplot.get_xticklabels(), rotation=45, horizontalalignment='right')
        
#         # Adding counts on the bars
#         for index, row in counts_df.iterrows():
#             barplot.text(row.name, row.Counts, row.Counts, color='black', ha="center")
        
#         pdf.savefig()  # Saves the current figure into the PDF
#         plt.close()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import os
from dia_sis.pipeline.utils import manage_directories



def create_report(df, path, params):
    # Construct the description string
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_strict_filters'].items()])
    description = f"Parameters used:\n{params_str}"

    
    # Construct the description string
    # Load filtering parameters from JSON
    CONFIG_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'configs')
    json_path = os.path.join(CONFIG_DIR, 'filtering_parameters.json')
    with open(json_path, 'r') as f:
        params = json.load(f)
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
    description = f"Parameters used:\n{params_str}"
    
    # Set up the PDF
    manage_directories.create_directory(path, 'reports')
    output_dir = os.path.join(path, 'reports')
    pdf_path = os.path.join(output_dir, 'protein_groups_report.pdf')

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(11, 8))
        plt.axis('off')

        plt.text(0.5, 0.98, "Protein Groups QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        
        # # Creating the barplot
        # plt.figure(figsize=(10, 7))
        # barplot = sns.barplot(x='Run', y='Counts', data=counts_df, orient='v')
        # plt.title('Counts per Run')
        # plt.xlabel('Run')
        # plt.ylabel('Counts')
        
        # Rotating x-axis labels
        # barplot.set_xticklabels(barplot.get_xticklabels(), rotation=45, horizontalalignment='right')
        
        # Adding counts on the bars
        # for index, row in counts_df.iterrows():
        #     barplot.text(row.name, row.Counts, row.Counts, color='black', ha="center")
        
        # pdf.savefig()  # Saves the current figure into the PDF
        # plt.text(0.5, 0.98, "Protein Ratios QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        # plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)
        # pdf.savefig()

        # plt.close()

        # Process the data for each 'Run'
        for run in df['Run'].unique():
            run_df = df[df['Run'] == run]
            counts = run_df.groupby('Run')[['H', 'L']].apply(lambda x: (x > 0).sum()).reset_index()
            
            # Plotting
            plt.figure(figsize=(10, 7))
            melted_counts = counts.melt(id_vars=['Run'], value_vars=['H', 'L'], var_name='Type', value_name='Counts')
            sns.barplot(x='Run', y='Counts', hue='Type', data=melted_counts)
            plt.title(f'Non-Zero Counts for Run {run}')
            plt.xlabel('Sample')
            plt.ylabel('Counts')
            plt.xticks(rotation=45)

            pdf.savefig()
            plt.close()

