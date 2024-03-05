# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 18:27:52 2024

@author: robbi
"""

import pandas as pd
import numpy as np
import operator
import time 

from tqdm import tqdm
import os
from icecream import ic

class Preprocessor:
    def __init__(self, path, params, filter_cols, meta_data=None):
        self.path = path
        self.meta_data = meta_data
        self.params = params
        self.chunk_size = 180000
        self.update = True
        self.filter_cols = filter_cols 
        
    def import_report(self):
        print('Beginning import report.tsv')
        start_time = time.time()
        # use dask to read in file quickly to count rows for tqdm loading bar
        
        chunks = []
        filtered_out = []
        contaminants = []
        file_path = f"{self.path}report.tsv"
        count = 1
        # Estimate rows in file size
        file_size_bytes = os.path.getsize(file_path)
        average_row_size_bytes = 1000  # This is an example; you'll need to adjust this based on your data
        # Estimate the number of rows
        estimated_rows = file_size_bytes / average_row_size_bytes
        total_chunks = estimated_rows/self.chunk_size
        
        for chunk in tqdm(pd.read_table(file_path, sep="\t", chunksize=self.chunk_size), 
                      total=total_chunks, desc='Estimated loading of report.tsv based on file size'):
                # reduce data size by subsetting report.tsv based on metadata, and removing columns not needed for further analysis
                # in the following loop we also annotate the silac chanels and append genes to Protein.Groups for downstream useage
                pd.options.mode.chained_assignment = None  # Turn off SettingWithCopyWarning since adaptions are being made to original df during import
                if self.meta_data is not None:
                    
                    chunk = self.subset_based_on_metadata(chunk)
                    chunk = self.relabel_run(chunk)
                    
                chunk['Genes'] = chunk['Genes'].fillna('')
                chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
                chunk['Label'] = ""
                chunk = self.add_label_col(chunk)
                # chunk = self.drop_non_valid_h_rows(chunk) probably dont need this func because filtering will remove these rows
                chunk = self.remove_cols(chunk)
                
                # annotate df with SILAC chanel then apply strict filters to H by droping the precursor, or adding NaN for L and M channels if they dont pass loose filters
                chunk, chunk_filtered_out = self.filter_channel(chunk, "H")
                chunk = self.apply_nan_by_loose_filtering(chunk,"L")
                
                contam_chunk = self.identify_contaminants(chunk)
                
                #remove filter cols before concatinating all dfs and returning
                chunk.drop(self.filter_cols, axis=1, inplace=True)
                chunks.append(chunk)
                filtered_out.append(chunk_filtered_out)
                contaminants.append(contam_chunk)
                
                # if self.update:
                #     print(f'Chunk {count} processed')
                # if count == 1:
                #     break
            
        # append chunks to respective dfs and return  
        df = pd.concat(chunks, ignore_index=True)
        filtered_out_df = pd.concat(filtered_out, ignore_index=True)
        contaminants_df = pd.concat(contaminants, ignore_index=True)
        print('Finished import')
        end_time = time.time()
        print(f"Time taken for import: {end_time - start_time} seconds")
        return df, filtered_out_df, contaminants_df


    def subset_based_on_metadata(self, chunk):
        filtered_chunk = chunk[chunk['Run'].isin(self.meta_data['Run'])]
        return filtered_chunk
    
    def relabel_run(self, chunk):
        run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))

        # Apply the mapping to df2['Run'] and raise an error if a 'Run' value doesn't exist in df1
        chunk['Run'] = chunk['Run'].map(run_to_sample)
        if chunk['Run'].isna().any():
            raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
            
        return chunk

    def add_label_col(self, chunk):
        # Extract the label and add it as a new column
        chunk['Label'] = chunk['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]
    
        # Remove the '(SILAC-K|R-([HML]))' part from the 'Precursor.Id' string
        chunk['Precursor.Id'] = chunk['Precursor.Id'].str.replace(r'\(SILAC-(K|R)-[HML]\)', '', regex=True)
    
        return chunk

    def remove_cols(self, chunk):
        cols = ['Run', 'Protein.Group', 'Precursor.Id', 'Label',  'Precursor.Quantity','Ms1.Translated','Precursor.Translated'] + self.filter_cols 
        chunk = chunk[cols]
        return chunk

    def filter_channel(self, chunk, label):
        ops = {
            "==": operator.eq, "<": operator.lt, "<=": operator.le,
            ">": operator.gt, ">=": operator.ge
        }
    
        # Check if 'H' labeled rows are present
        if label in chunk['Label'].values:
            # Start with a mask that selects all chanel rows
            h_rows_mask = chunk['Label'] == label
    
            for column, condition in self.params['apply_loose_filters'].items():
                op = ops[condition['op']]
                # Update the mask to keep chanel rows that meet the condition
                h_rows_mask &= op(chunk[column], condition['value'])
    
            # Filter out chanel rows that do not meet all conditions
            filtered_chunk = chunk[h_rows_mask | (chunk['Label'] != label)]
            chunk_filtered_out = chunk[~h_rows_mask & (chunk['Label'] == label)]
        else:
            # If the label is not present, return the whole chunk and an empty 'filtered' chunk
            filtered_chunk = chunk
            chunk_filtered_out = pd.DataFrame(columns=chunk.columns)
    
        return filtered_chunk, chunk_filtered_out

    def apply_nan_by_loose_filtering(self, chunk, label):
        # Create boolean mask for length of chunk
        filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        ops = {
            "==": operator.eq, "<": operator.lt, "<=": operator.le,
            ">": operator.gt, ">=": operator.ge
        }
        filtered_chunk = chunk.copy()
    
        if label in chunk['Label'].values:
            for column, condition in self.params['apply_loose_filters'].items():
                op = ops[condition['op']]
                # Update the condition for each column if any of the filtering criterea don't pass
                filtering_condition &= op(chunk[column], condition['value'])
    
            nan_cols = ['Precursor.Translated', 'Precursor.Quantity', 'Ms1.Translated']
            for col in nan_cols:
                if col in chunk.columns:
                    # Set NaN where the condition is False
                    filtered_chunk.loc[~filtering_condition, col] = np.nan
    
        return filtered_chunk

    def identify_contaminants(self, chunk):
         chunk_copy = chunk.copy(deep=True)
         contams_mask = chunk_copy['Protein.Group'].str.contains('Cont_', case=False, na=False)
         self._validate_boolean_mask(contams_mask)
              
         contaminants = chunk_copy[contams_mask]
         return contaminants
    
    def _validate_boolean_mask(self, mask):
        if not all(isinstance(x, bool) for x in mask):
            invalid_values = mask[~mask.isin([True, False])]
            print(f"Non-boolean values in mask: {invalid_values}")
            
            
   