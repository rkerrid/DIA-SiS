# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 08:41:02 2023

@author: rkerrid
"""
import pandas as pd
pd.set_option('display.max_columns', None)
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

from icecream import ic
from matplotlib.backends.backend_pdf import PdfPages
from fpdf import FPDF
import os

import sys
print(sys.path)

from silac_dia_tools.pipeline.pipeline import Pipeline

test_data_bm = 'C:/data/silac_dia_tools_files/data/BM data/' 
test_data_bm = 'C:/phd projects/silac_dia_tools/test data/'                  

def filter_and_rename(df, string, rep_string):
    cols = [col for col in df.columns if string in col and rep_string in col]
    df_filtered = df[['Protein.Group'] + cols]
    return df_filtered.rename(columns={df_filtered.columns[1]: 'Measurement'})

def calculate_ratios(df_S4, df_S5):
    df_S5['ratio'] = df_S4['Measurement'] / df_S5['Measurement']
    df_S5['ratio'] = df_S5['ratio'].replace([np.inf, -np.inf], np.nan)
    df_S5['log_Measurement'] = np.log10(df_S5['Measurement'])
    df_S5['log_ratio'] = np.log10(df_S5['ratio'])
    return df_S5

def shift_ratios(df, mask):
    median_log_human = df.loc[mask, 'log_ratio'].median()
    df['log_shifted_ratio'] = df['log_ratio'] - median_log_human
    return df, median_log_human

def plot_data(ax, df, mask, color, label, x_col, y_col, plot_median=False):
    ax.scatter(df.loc[mask, x_col], df.loc[mask, y_col], color=color, alpha=0.5)
    ax.set_xlabel('log10(Intensity)', fontsize=16)
    ax.set_ylabel('log10(S4:S5)', fontsize=16)
    ax.grid(True)
    
    if plot_median:
        median_value = df.loc[mask, y_col].median()
        ax.axhline(y=median_value, color='grey', linestyle='-', linewidth=4, alpha=0.5, label=f'Median {label}')
    ax.legend()

def add_href_intensity(df):
    df = df.copy()
    href = pd.read_csv(f'{test_data_bm}protein intensities/reference_href.csv')
    href = href[['Protein.Group','S5']] # get the x asiy for all plots which will be the href from S5 rep3, df is Protein.Group, 1a
    href.rename(columns = {'S5':'href'}, inplace=True)
    
    df = df.merge(href[['Protein.Group', 'href']],
                  on='Protein.Group',
                  how='left')
    df['href'] = np.log10(df['href'])
     
    return df
    
def plot_ratios(df, expected_human, expected_ecoli, title):
 
    df_S4 = df[['Protein.Group','S4']]
    df_S4 = df_S4.rename(columns={'S4':'Measurement'})
    df_S5 = df[['Protein.Group','S5']]
    df_S5 = df_S5.rename(columns={'S5':'Measurement'})
    df_S5 = calculate_ratios(df_S4, df_S5)
    
    mask_ecoli = df_S5['Protein.Group'].str.contains('ECOLI_')
    mask_human = ~mask_ecoli
    df_S5, median_log_human = shift_ratios(df_S5, mask_human)
   
    df_S5 = add_href_intensity(df_S5)
    
    fig, ax = plt.subplots(1, 2, figsize=(20, 8), sharey=True)
    fig.suptitle(title, fontsize=22, fontweight='bold')
    
    ax[0].set_title('Ecoli Subset', fontsize=16)
    plot_data(ax[0], df_S5, mask_ecoli, 'red', 'Ecoli', 'href', 'log_shifted_ratio', plot_median=True)
    
    ax[1].set_title('Human Subset', fontsize=16)
    plot_data(ax[1], df_S5, mask_human, 'blue', 'Human', 'href', 'log_shifted_ratio')
    
    ax[0].set_ylim(-1.5, 1)
    ax[0].set_xlim(2.5, 7)
    ax[1].set_xlim(2.5, 7)
    for a in ax:
        a.axhline(y=np.log10(expected_human) , color='blue', linestyle='--', label='Expected Human')
        a.axhline(y=np.log10(expected_ecoli) , color='red', linestyle='--', label='Expected Ecoli')
        a.legend()

    plt.show()
    return fig  # return the figure object

def save_plot_to_pdf(df,  expected_human, expected_ecoli, title, filename):
    fig = plot_ratios(df, expected_human, expected_ecoli, title)  # Get the figure object
    os.makedirs(os.path.dirname(filename), exist_ok=True)  
    fig.savefig(filename)  # Save the figure object directly
    plt.close(fig)  # Close the specific figure





def analyze_results():
    # # import dlfq normalized proteomes and compare s5:s4 ecoli and human ratios rep3
    df_light_lfq = pd.read_csv(f'{test_data_bm}protein intensities/light_dlfq.csv')
    df_heavy_lfq = pd.read_csv(f'{test_data_bm}protein intensities/nsp_dlfq.csv')
    df_total_lfq = pd.read_csv(f'{test_data_bm}protein intensities/total_dlfq.csv')
    
    # # import href normalized proteomes and compare s5:s4 ecoli and human ratios rep3
    df_light_href = pd.read_csv(f'{test_data_bm}protein intensities/light_href.csv')
    df_heavy_href = pd.read_csv(f'{test_data_bm}protein intensities/reference_href.csv')
    df_total_href = df_heavy_href.copy()
    df_total_href.iloc[:,1:] = df_heavy_href.iloc[:,1:].add(df_light_href.iloc[:,1:])
    
    
    # # import unnormalized proteomes and compare s5:s4 ecoli and human ratios rep3
    df_light = pd.read_csv(f'{test_data_bm}protein intensities/light_unnorm.csv')
    df_heavy = pd.read_csv(f'{test_data_bm}protein intensities/reference_unnorm.csv')
    df_total = df_heavy.copy()
    df_total.iloc[:,1:] = df_heavy.iloc[:,1:].add(df_light.iloc[:,1:])
    
    # workings for human abundance in each sample
    human_light_s5 = 66.8 
    human_light_s4 = 66.8 
    human_heavy = 60
    human_total_s5 = human_light_s5 + human_heavy
    human_total_s4 = human_light_s4 + human_heavy
    # workings for ecoli abundance in each sample
    ecoli_light_s5 = 66.8 
    ecoli_light_s4 = 6.8
    ecoli_heavy = 6
    ecoli_total_s5 = ecoli_light_s5 + ecoli_heavy
    ecoli_total_s4 = ecoli_light_s4 + ecoli_heavy
    # workings for expected ratios in each channel (human)
    expected_human_light = human_light_s4/human_light_s5
    expected_human_heavy = human_heavy/human_heavy
    expected_human_total = human_total_s4/human_total_s5
    # workings for expected ratios in each channel (ecoli)
    expected_ecoli_light = ecoli_light_s4/ecoli_light_s5
    expected_ecoli_heavy = ecoli_heavy/ecoli_heavy
    expected_ecoli_total = ecoli_total_s4/ecoli_total_s5
    
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size = 12)
    pdf.cell(200, 10, txt = "Bench Mark Ratio Report Comparing LFQ, Heavy Reference, and Unnormalized Data", ln = True, align = 'C')
    
    # Add the variables block
    vars_block = """
    # workings for human abundance in each sample
    human_light_s5 = 66.8 
    human_light_s4 = 66.8 
    human_heavy = 60
    human_total_s5 = human_light_s5 + human_heavy
    human_total_s4 = human_light_s4 + human_heavy
    # workings for ecoli abundance in each sample
    ecoli_light_s5 = 66.8 
    ecoli_light_s4 = 6.8
    ecoli_heavy = 6
    ecoli_total_s5 = ecoli_light_s5 + ecoli_heavy
    ecoli_total_s4 = ecoli_light_s4 + ecoli_heavy
    # workings for expected ratios in each channel (human)
    expected_human_light = human_light_s4/human_light_s5
    expected_human_heavy = human_heavy/human_heavy
    expected_human_total = human_total_s4/human_total_s5
    # workings for expected ratios in each channel (ecoli)
    expected_ecoli_light = ecoli_light_s4/ecoli_light_s5
    expected_ecoli_heavy = ecoli_heavy/ecoli_heavy
    expected_ecoli_total = ecoli_total_s4/ecoli_total_s5
    """
    pdf.multi_cell(0, 10, txt = vars_block)
    
    # Save the plots to individual files and add them to the PDF
    plots = [
        ("df_light", "expected_human_light", "expected_ecoli_light", "Light (no norm)"),
        ("df_heavy", "expected_human_heavy", "expected_ecoli_heavy", "Heavy (no norm)"),
        ("df_total", "expected_human_total", "expected_ecoli_total", "Total (no norm)"),
        
        ("df_light_lfq", "expected_human_light", "expected_ecoli_light", "Light (lfq)"),
        ("df_heavy_lfq", "expected_human_heavy", "expected_ecoli_heavy", "Heavy (lfq)"),
        ("df_total_lfq", "expected_human_total", "expected_ecoli_total", "Total (lfq)"),
        
        ("df_light_href", "expected_human_light", "expected_ecoli_light", "Light (href)"),
        ("df_heavy_href", "expected_human_heavy", "expected_ecoli_heavy", "Heavy (href)"),
        ("df_total_href", "expected_human_light", "expected_ecoli_total", "Total (href)")
        
    ]
    
    
    pdf_output_path = os.path.join(test_data_bm, "benchmark_report.pdf")
    
    for df_name, human_exp, ecoli_exp, title in plots:
        plot_filename = os.path.join(test_data_bm, f"{title.replace(' ', '_')}.png")
        
        try:
            # Save plot
            save_plot_to_pdf(eval(df_name), eval(human_exp), eval(ecoli_exp), title, plot_filename)
            
            # Debugging: Check if file exists and print its size
            if os.path.exists(plot_filename):
                print(f"Plot {title} saved successfully at {plot_filename}, Size: {os.path.getsize(plot_filename)} bytes")
            else:
                print(f"Failed to save the plot: {title}")
                continue  # Skip to next iteration if the plot was not saved
            
            # Add to PDF
            pdf.add_page()
            pdf.image(plot_filename, x = 10, y = 10, w = 180)  # Specify size and position
    
            
            # Remove plot file
            os.remove(plot_filename)  
        except Exception as e:
            print(f"An error occurred while processing {title}: {str(e)}")
    
    
    
    pdf.output(pdf_output_path)
    
    
 
# Process diann output files: filtering, formatting silac precursors, ratios, intensities (directLFQ) with 'H pulse'
pipeline = Pipeline(test_data_bm, 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="H")

pipeline.run_dlfq_pipeline()
pipeline.run_href_pipeline()
pipeline.generate_reports()

analyze_results()