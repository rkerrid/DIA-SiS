import pandas as pd
# Turn off the SettingWithCopyWarning
pd.options.mode.chained_assignment = None
import operator
import time 
from tqdm import tqdm
import os


class Preprocessor:
    def __init__(self, path, params, filter_cols, requantify, meta_data=None):
        self.path = path
        self.meta_data = meta_data
        self.params = params
        self.chunk_size = 180000
        self.update = True
        self.filter_cols = filter_cols 
        self.requantify = requantify
        
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
        average_row_size_bytes = 800  
        # Estimate the number of rows
        estimated_rows = file_size_bytes / average_row_size_bytes
        total_chunks = estimated_rows/self.chunk_size
        
        for chunk in tqdm(pd.read_table(file_path, sep="\t", chunksize=self.chunk_size), 
                      total=total_chunks, desc='Estimated loading of report.tsv based on file size'):
                # reduce data size by subsetting report.tsv based on metadata, and removing columns not needed for further analysis
                # in the following loop we also annotate the silac chanels and append genes to Protein.Groups for downstream useage
                if self.meta_data is not None:
                    chunk = self._subset_based_on_metadata(chunk) # double check
                    chunk = self._relabel_run(chunk)
                    
                chunk['Genes'] = chunk['Genes'].fillna('')
                chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
                chunk['Label'] = ""
                chunk = self._add_label_col(chunk)
                chunk = self._remove_cols(chunk)
                
                chunk, chunk_filtered_out = self._filter_channel(chunk, "H")
                
                # requantify settings set to False will filter light precursors and concat the filtered out precursors with chunck_filtered_out 
                if not self.requantify:
                    chunk, filtered_out_light = self._filter_channel(chunk, "L")
                    chunk_filtered_out = pd.concat([chunk_filtered_out, filtered_out_light], ignore_index=True)
                    print('not requantify')

                # ID contaminants for report
                # contam_chunk = self._identify_contaminants(chunk)
                chunk, contam_chunk = self._remove_contaminants(chunk)
                
                #remove filter cols before concatinating all dfs and returning filtered df for reports and protein roll up
                chunk.drop(self.filter_cols, axis=1, inplace=True)
                chunks.append(chunk)
                filtered_out.append(chunk_filtered_out)
                contaminants.append(contam_chunk)
                
                # if count == 1: # unncomment for troubleshooting
                #     break
            
        # append chunks to respective dfs and return  
        df = pd.concat(chunks, ignore_index=True)
        filtered_out_df = pd.concat(filtered_out, ignore_index=True)
        contaminants_df = pd.concat(contaminants, ignore_index=True)
        print('Finished import')
        end_time = time.time()
        print(f"Time taken for import: {end_time - start_time} seconds")
        return df, filtered_out_df, contaminants_df


    def _subset_based_on_metadata(self, chunk):
        filtered_chunk = chunk[chunk['Run'].isin(self.meta_data['Run'])]
        return filtered_chunk
    
    def _relabel_run(self, chunk):
        run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))

        # Apply the mapping to the current chunk and raise an error if a 'Run' value from the metadata doesn't exist in this chunk
        chunk['Run'] = chunk['Run'].map(run_to_sample)
        if chunk['Run'].isna().any():
            raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
            
        return chunk

    def _add_label_col(self, chunk):
        # Extract the label and add it as a new column
        chunk['Label'] = chunk['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]
    
        # Remove the '(SILAC-K|R-([HML]))' part from the 'Precursor.Id' string
        chunk['Precursor.Id'] = chunk['Precursor.Id'].str.replace(r'\(SILAC-(K|R)-[HML]\)', '', regex=True)
    
        return chunk

    def _remove_cols(self, chunk):
        cols = ['Run', 'Protein.Group', 'Precursor.Id', 'Label',  'Precursor.Quantity','Ms1.Translated','Precursor.Translated'] + self.filter_cols 
        chunk = chunk[cols]
        return chunk

    def _filter_channel(self, chunk, label):
        ops = {
            "==": operator.eq, "<": operator.lt, "<=": operator.le,
            ">": operator.gt, ">=": operator.ge
        }
    
        # Check if label is present in the chunk
        if label in chunk['Label'].values:
            # Start with a mask that selects all chanel rows
            channel_rows_mask = chunk['Label'] == label
    
            for column, condition in self.params['apply_loose_filters'].items(): # can set this to 'strict filtering', see params.json in configs folder
                op = ops[condition['op']]
                # Update the mask to keep chanel rows that meet the condition
                channel_rows_mask &= op(chunk[column], condition['value'])
    
            # Filter out chanel rows that do not meet all conditions
            filtered_chunk = chunk[channel_rows_mask | (chunk['Label'] != label)]
            chunk_filtered_out = chunk[~channel_rows_mask & (chunk['Label'] == label)]
        else:
            # If the label is not present, return the whole chunk and an empty 'filtered' chunk
            filtered_chunk = chunk
            chunk_filtered_out = pd.DataFrame(columns=chunk.columns)
    
        return filtered_chunk, chunk_filtered_out


    def _identify_contaminants(self, chunk):
         chunk_copy = chunk.copy(deep=True)
         contams_mask = chunk_copy['Protein.Group'].str.contains('Cont_', case=False, na=False)
         # self._validate_boolean_mask(contams_mask)
              
         contaminants = chunk_copy[contams_mask]
         return contaminants
     
    def _remove_contaminants(self, chunk):
         # chunk_copy = chunk.copy(deep=True)
         contams_mask = chunk['Protein.Group'].str.contains('Cont_', case=False, na=False)
         # self._validate_boolean_mask(contams_mask)
              
         contaminants = chunk[contams_mask]
         chunk = chunk[~contams_mask]
         return chunk, contaminants
    
    # def _validate_boolean_mask(self, mask):
    #     if not all(isinstance(x, bool) for x in mask):
    #         invalid_values = mask[~mask.isin([True, False])]
    #         print(f"Non-boolean values in mask: {invalid_values}")
            
            
   