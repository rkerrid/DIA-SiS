# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:37:19 2024

@author: robbi
"""


from icecream import ic
# from dia_sis.pipeline.pipeline import Pipeline as pileline
import pandas as pd 
import numpy as np

df = pd.read_csv('G:/My Drive/data/main experiments/20240219 baby benchmark for pydia_sis/preprocessing/formatted_precursors.tsv', sep = '\t')
protein = 'ECOLI_P0A790;P0A790'
df = df[df['Protein.Group'].str.contains('ECOLI_P0A790;P0A790')]


def calculate_precursor_href_intensities( df):
    
    def combined_median(ms1_series, precursor_series):
        # Replace invalid values with NaN and drop them
        valid_ms1 = ms1_series.replace([0, np.inf, -np.inf], np.nan).dropna()
        valid_precursor = precursor_series.replace([0, np.inf, -np.inf], np.nan).dropna()
   
        # Ensure at least 3 valid values in each series before combining
        if len(valid_ms1) >= 1 or len(valid_precursor) >= 1:
            combined_series = np.concatenate([valid_ms1, valid_precursor])
            combined_series = np.log10(combined_series)  # Log-transform the combined series
            return np.median(combined_series)  # Return the median of the log-transformed values
        else:
            return np.nan
   
    # Group by protein group and apply the custom aggregation
    grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
        'href': combined_median(x['Ms1.Translated H'], x['Precursor.Translated H']) 
    })).reset_index()
   
    return grouped[['Protein.Group', 'href']]

result = calculate_precursor_href_intensities(df)