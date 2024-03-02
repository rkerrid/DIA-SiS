# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 08:20:35 2024

@author: rkerrid
"""



import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import os
from dia_sis.pipeline.utils import manage_directories



def create_report(path, params):
    # Construct the description string
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_strict_filters'].items()])
    description = f"Parameters used:\n{params_str}"
    
    # Set up the PDF
    manage_directories.create_directory(path, 'reports')
    output_dir = os.path.join(path, 'reports')
    pdf_path = os.path.join(output_dir, 'protein_groups_report.pdf')
    
    # Read in each dataframe  
    # Load the data from all three CSV files
    file_paths = {
        'reference': f'{path}/protein_groups/href_href.csv',
        'light_proteome': f'{path}/protein_groups/light_href.csv',
        'nsp_proteome': f'{path}/protein_groups/nsp_href.csv'
    }
    
    # Initialize a dictionary to hold the dataframes
    dataframes = {}
    
    # Load each CSV into a dataframe and store it in the dictionary
    for key, path in file_paths.items():
        dataframes[key] = pd.read_csv(path)
    
   

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(11, 8))
        plt.axis('off')

        plt.text(0.5, 0.98, "Protein Groups QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        
        # Function to count non-zero values in dataframe columns, excluding the first column
        def count_non_zeros(df):
            return df[df.columns[1:]].apply(lambda x: x[pd.notnull(x) | (x != 0)].count())
            
        # Count non-zero values for each dataset
        non_zero_counts = {key: count_non_zeros(df) for key, df in dataframes.items()}
        
        # Since we assume the column names are consistent across files, we'll use the column names from one of the datasets for plotting
        sample_columns = non_zero_counts['reference'].index
        
        # Prepare data for plotting
        plot_data = pd.DataFrame({key: non_zero_counts[key] for key in non_zero_counts.keys()}, index=sample_columns)
        
        # Plotting
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Colors for each dataset
        colors = ['blue', 'orange', 'green']
        labels = list(non_zero_counts.keys())
        
        # Create bar plots
        plot_data.plot(kind='bar', ax=ax, color=colors, width=0.8)
        
        # Customize the plot
        plt.title('Comparison of Non-Zero Value Counts Across Datasets')
        plt.xlabel('Sample Columns')
        plt.ylabel('Non-Zero Values Count')
        plt.xticks(rotation=45, ha="right")
        plt.legend(labels)
        
        # Tight layout for better spacing
        plt.tight_layout()
        
        pdf.savefig()
        plt.close()
        
        # Now plot cummulative intensities 
        # Function to sum values in dataframe columns, excluding the first column
        def sum_values(df):
            return df[df.columns[1:]].apply(lambda x: x.sum())
        
        # Sum values for each dataset
        sums = {key: sum_values(df) for key, df in dataframes.items()}
        
        # Prepare data for plotting
        sum_data = pd.DataFrame({key: sums[key] for key in sums.keys()}, index=sample_columns)
        
        # Plotting
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Create bar plots
        sum_data.plot(kind='bar', ax=ax, color=colors, width=0.8)
        
        # Customize the plot
        plt.title('Comparison of Summed Values Across Datasets')
        plt.xlabel('Sample Columns')
        plt.ylabel('Summed Values')
        plt.xticks(rotation=45, ha="right")
        plt.legend(labels)
        
        # Tight layout for better spacing
        plt.tight_layout()
        
        pdf.savefig()
        plt.close()

      


