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
        
    def generate_protein_groups(self):
        start_time = time.time()
        # Will use Precursor.Translated for quantification
        quantification = 'Precursor.Translated'
        # formatting and ratios
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors, quantification)
        self.protein_groups = self.compute_protein_level_test(self.formatted_precursors)
        # self.protein_groups = self.compute_protein_level_test(self.formatted_precursors)

        # Adjusting intensities and outputing data
        self.protein_groups = self.calculate_href_intensities(self.protein_groups)
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
    
        # remove all rows that contain a NaN under the Label H column (i.e., no H precursor is present for that row)
        # Apply dropna on merged_df instead of df
        merged_df = merged_df.dropna(subset=['Precursor.Translated H'])
        # replace precursor quantity with summed silac channels as input for direct lefq and as 'total intensity' for href quantification
        merged_df['Precursor.Quantity'] = merged_df['Precursor.Translated H'] + merged_df['Precursor.Translated L'] 
 
        return merged_df
    
    def calculate_precursor_ratios(self, df, quantification):
        print(f'Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        # Calculate ratios for all chanels (Precursor.Quantity is the total intensity of all 3 chanels, the default diann value has been overwritten at this point)
        df['Precursor.Translated H/T'] = df['Precursor.Translated H'] / (df['Precursor.Translated L'] + df['Precursor.Translated H'])
        df['Ms1.Translated H/T'] = df['Ms1.Translated H'] / (df['Ms1.Translated L'] + df['Ms1.Translated H'])
        return df
    
    def compute_protein_level(self, df): # this function looks for at least 3 valid values for each ratio and sums Precursor.Quantity (which is the sum of precursor translated values) for total intensity
        print('Rolling up to protein level')
        # Function to filter values and compute median
        def valid_median(series):
            valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
            valid_series = np.log2(valid_series)
            return 2 **valid_series.median() if len(valid_series) >= 2 else np.nan
        
        def valid_sum(series):
            valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
            return valid_series.sum() 
        
        result = df.groupby(['Protein.Group', 'Run']).agg({
            'Precursor.Translated H/T': valid_median,
            'Precursor.Translated L/T': valid_median,
            'Precursor.Quantity': valid_sum        
        })
        result['H'] = result['Precursor.Translated H/T']*result['Precursor.Quantity']
        result['L'] = result['Precursor.Translated L/T']*result['Precursor.Quantity']
        result = result.reset_index()
        
        cols = ['Run', 'Protein.Group', 'Precursor.Quantity', 'H', 'L'] 
        return result[cols]
    

    def compute_protein_level_test(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        for run in runs:
            run_df = df[df['Run'] == run]
    
            def combined_median(ms1_series, precursor_series):
                # Replace invalid values with NaN and drop them
                valid_ms1 = ms1_series.replace([0, np.inf, -np.inf], np.nan).dropna()
                valid_precursor = precursor_series.replace([0, np.inf, -np.inf], np.nan).dropna()
    
                # Ensure at least 3 valid values in each series before combining
                if len(valid_ms1) >= 2 and len(valid_precursor) >= 2:
                    combined_series = np.concatenate([valid_ms1, valid_precursor])
                    combined_series = np.log2(combined_series)  # Log-transform the combined series
                    return 2 ** np.median(combined_series)  # Return the median of the log-transformed values
                else:
                    return np.nan
    
            def valid_sum(series):
                valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
                return valid_series.sum()
    
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
                'H/T ratio': combined_median(x['Ms1.Translated H/T'], x['Precursor.Translated H/T']),
                'Precursor.Quantity': valid_sum(x['Precursor.Quantity'])
            })).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)
        result['H'] = result['H/T ratio']*result['Precursor.Quantity']
        result['L'] = result['Precursor.Quantity'] - result['H']
        cols = ['Run','Protein.Group', 'H', 'L', 'H/T ratio', 'Precursor.Quantity']
    
        # Returning the dataframe with specified columns
        return result[cols]





    # Adjust unnormalized intensities
    def calculate_href_intensities(self, df):
        print('Calculating adjusted intensities using reference')
        df_copy = df.copy(deep=True)
        
        # Calculate median H value and reset index to make it a DataFrame
        h_ref = df_copy.groupby('Protein.Group')['H'].median().reset_index()
        
        # Rename the median column to 'h_ref'
        h_ref = h_ref.rename(columns={'H': 'h_ref'})
        
        # Merge the original DataFrame with the h_ref DataFrame
        merged_df = df.merge(h_ref, on='Protein.Group', how='inner')
        
        # calculate factor to multiply other chanels by dividing href by original H intensity for each PG
        merged_df['factor'] = merged_df['h_ref']/merged_df['H']
        
        # Normalize other chanels with this factor
        merged_df['H_norm'] = merged_df['H']*merged_df['factor']
        merged_df['L_norm'] = merged_df['L']*merged_df['factor']
        
        return merged_df
    
    def output_protein_groups(self, df, path):
        manage_directories.create_directory(self.path, 'protein_groups')
        print(f'Outputing normalized protein intensities to {path}/protein_groups')
        cols = ['Run', 'Protein.Group', 'H_norm', 'L_norm']
        df = df[cols]
        df = df.rename(columns={'H_norm': 'H', 'L_norm': 'L'})

        # Pivoting for 'H'
        h_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='H')
        
        
        # Pivoting for 'L'
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')
        
        # then output each table to csv for h.href, l.href, m.href
        h_pivot_df.to_csv(f'{path}/protein_groups/href_href.csv', sep=',')
        l_pivot_df.to_csv(f'{path}/protein_groups/light_href.csv', sep=',')

        return h_pivot_df, l_pivot_df
