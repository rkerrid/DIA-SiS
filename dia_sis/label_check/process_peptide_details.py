# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:15:43 2023

@author: robbi
"""
import pandas as pd
from icecream import ic
ic.disable()

#combine data from MQ and thermo raw files with the option to save to a temp csv when adjusting plotting
def process_details(peptides, ms1_dfs, path):
    for peptide_details in peptides:
        raw_file = peptide_details['Raw_file'].item() # get key for raw files dictionary
        ic(peptide_details)
        ic(raw_file)
        ms1_df = ms1_dfs[raw_file]
        mz_values,intensity_values = extract_values(ms1_df, peptide_details)
        sequence = peptide_details['Sequence'].values[0]
        print(sequence)
        data = {"mz": mz_values,
                "intensity": intensity_values
                }
        df = pd.DataFrame(data)
        df.to_csv(f"{path}{sequence}.csv", index=False)
    
    
## Extracting raw file data 
def extract_values(ms1_df, peptide_details): #will break here
    scan_number = peptide_details["scan_number"].values[0]
    scan = ms1_df[ms1_df['scanms1']==scan_number][["mz_values_ms1","intensity_values_ms1", "rtms1"]]

    # Extract lists from DataFrame
    x_values = scan['mz_values_ms1'].tolist()
    y_values = scan['intensity_values_ms1'].tolist()
    x_values = [item for sublist in x_values for item in sublist]
    y_values = [item for sublist in y_values for item in sublist]
   
    return x_values, y_values

