# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 13:33:49 2023

@author: rkerrid
"""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os


from icecream import ic
import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from dia_sis.pipeline.utils import manage_directories



def create_correlation_heatmap(file_path, pdf):

    # Load the dataset
    df = pd.read_csv(file_path, index_col=0)  # Assuming the first column is 'Protein Group'

    # Calculate the correlation matrix
    corr_matrix = df.corr(method='pearson')

    # Plot the heatmap
    f, ax = plt.subplots(figsize=(11, 9))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm')
  
    # Add 1 to last_slash_index to start after the "/", and csv_index is the end point
    title = file_path[file_path.rfind('/') + 1 : file_path.find('.csv')]
     
    plt.title(f"Heatmap of Sample Correlations in {title}")
    # Save the heatmap to the PDF
    pdf.savefig(f)
    plt.close(f)


def create_report(path, params, href=True):
    # Construct the description string
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
    description = f"Parameters used:\n{params_str}"

    # Set up the PDF
    manage_directories.create_directory(path, 'reports')
    output_dir = os.path.join(path, 'reports')

    pdf_path = os.path.join(output_dir, 'protein_intensities_report.pdf')

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(11, 8))
        plt.axis('off')

        plt.text(0.5, 0.98, "Protein Intensities QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        
        path = f'{path}protein intensities/'
        if href:
            file_list = [f'{path}light_href.csv', f'{path}light_unnorm.csv', f'{path}nsp_href.csv', f'{path}nsp_unnorm.csv'] 
        else:
            file_list = [f'{path}light_dlfq.csv', f'{path}light_unnorm.csv', f'{path}nsp_dlfq.csv', f'{path}nsp_unnorm.csv']        

        for file in file_list:
            
            # Add correlation heatmap page
            create_correlation_heatmap(file, pdf)




