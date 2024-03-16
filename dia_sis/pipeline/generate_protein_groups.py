# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:49:28 2024

@author: robbi
"""
import pandas as pd
import numpy as np
import time
from icecream import ic
from .utils import manage_directories


class HrefRollUp:
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        
        self.formatted_precursors = None
        self.protein_groups = None
        self.href_df = None
        
    def generate_protein_groups(self):
        start_time = time.time()
        # formatting and ratios
        self.formatted_precursors = self.format_silac_channels(self.filtered_report) ###
        # ic(self.formatted_precursors)
        # remove rows where precursors dont have matching valid values
        self.formatted_precursors = self.keep_only_shared_intensities(self.formatted_precursors)
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors)
        # ic(self.formatted_precursors)
        self.href_df = self.calculate_precursor_href_intensities(self.formatted_precursors)
        self.protein_groups = self.compute_protein_level(self.formatted_precursors)

        self.protein_groups = self.href_normalization(self.protein_groups, self.href_df) #uses precursor median

        self.output_protein_groups(self.protein_groups, self.path)
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        return self.formatted_precursors, self.protein_groups
    
    def format_silac_channels(self, df):
        print('Formatting SILAC channels')
        # Pivot for each label
        pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
        pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' H')
        
        # Merge the pivoted DataFrames
        merged_df = pd.concat([pivot_L, pivot_H], axis=1)
        
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        merged_df.reset_index(inplace=True)
        return merged_df
    
    
    def keep_only_shared_intensities(self, formatted_precursors):
        # Replace 0, inf, -inf with NaN for the specified columns
        formatted_precursors['Precursor.Translated H'] = formatted_precursors['Precursor.Translated H'].replace([0, np.inf, -np.inf], np.nan).dropna()
        formatted_precursors['Precursor.Translated L'] = formatted_precursors['Precursor.Translated L'].replace([0, np.inf, -np.inf], np.nan).dropna()
        
        formatted_precursors['Ms1.Translated H'] = formatted_precursors['Ms1.Translated H'].replace([0, np.inf, -np.inf], np.nan).dropna()
        formatted_precursors['Ms1.Translated L'] = formatted_precursors['Ms1.Translated L'].replace([0, np.inf, -np.inf], np.nan).dropna()
    
        # Keep only rows where both 'Precursor Translated H' and 'Precursor Translated L' are present (non-NaN)
        filtered_df = formatted_precursors.dropna(subset=['Precursor.Translated H','Ms1.Translated H'])
        filtered_df = formatted_precursors.dropna(subset=['Precursor.Translated L','Ms1.Translated L'])
           
        # can also filter so that light and heavy are in the same row
        return filtered_df

    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        df['Precursor.Translated L/H'] = df['Precursor.Translated L'] / df['Precursor.Translated H'] 
        df['Ms1.Translated L/H'] = df['Ms1.Translated L'] / df['Ms1.Translated H'] 
      
        return df

    def calculate_precursor_href_intensities(self, df):
        def combined_median(ms1_series, precursor_series):
            # if len(ms1_series) >= 1 and len(precursor_series) >= 1:
            combined_series = np.concatenate([ms1_series, precursor_series])
            combined_series = np.log10(combined_series)  # Log-transform the combined series
            return np.median(combined_series)  # Return the median of the log-transformed values
       
        # Group by protein group and apply the custom aggregation
        grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
            'href': combined_median(x['Ms1.Translated H'], x['Precursor.Translated H']) 
        })).reset_index()
       
        return grouped[['Protein.Group', 'href']]
    
   
    def compute_protein_level(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        for run in runs:
            run_df = df[df['Run'] == run]
    
            def combined_median_ratios(ms1_series, precursor_series):
                
                valid_ms1 = ms1_series.replace([0, np.inf, -np.inf], np.nan).dropna()
                valid_precursor = precursor_series.replace([0, np.inf, -np.inf], np.nan).dropna()
                
                # Ensure at least 1 valid value in either series before combining
                # if len(valid_ms1) >= 1 and len(valid_precursor) >= 1:
                combined_series = np.concatenate([valid_ms1, valid_precursor])
                combined_series = np.log10(combined_series)  # Log-transform the combined series
                return np.median(combined_series)  # Return the median of the log-transformed values
                # else:
                #     return np.nan
    
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby('Protein.Group').apply(lambda x: pd.Series({
                'L/H ratio': combined_median_ratios(x['Ms1.Translated L/H'], x['Precursor.Translated L/H'])
            })).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)
        cols = ['Run', 'Protein.Group', 'L/H ratio']
     
        # result[cols].to_csv('G:/My Drive/Data/main experiments/protein_groups_unnormalized_log10.csv', sep=',')
        return result[cols]

    
    def href_normalization(self, protein_groups, href):
        print('Calculating adjusted intensities using reference')
        # Merge the href_df onto protein groups containing optimized ratios
        # ic(protein_groups)
        # ic(href)
        merged_df = protein_groups.merge(href, on='Protein.Group', how='inner')
        
        # Obtain normalized light intensities by adding the L/H ratio to the heavy refference in log space
        merged_df['L_norm'] = merged_df['L/H ratio'] + merged_df['href']
        
        
        # # testing section
        
        # merged_df = protein_groups
        # merged_df['L_norm'] == merged_df['L/H ratio'] + 3
        return merged_df
    
    
    def output_protein_groups(self, df, path):
        manage_directories.create_directory(self.path, 'protein_groups')
        print(f'Outputing normalized protein intensities to {path}/protein_groups')
        
        # Subset and rename columns
        df = df[['Run', 'Protein.Group', 'href', 'L_norm']]
        df = df.rename(columns={'href': 'H', 'L_norm': 'L'})

        # Pivoting for 'H' to produce href output in wide format
        h_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='H')
        
        # Pivoting for 'L' to produce normalized light intensites in wide format
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')

        # then output each table to csv 
        h_pivot_df.to_csv(f'{path}/protein_groups/href.csv', sep=',')
        l_pivot_df.to_csv(f'{path}/protein_groups/light.csv', sep=',')

        return h_pivot_df, l_pivot_df










