# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 13:09:52 2023

@author: rkerrid
"""

## Extracting raw file data 
def extract_values(ms1_df, peptide_details):
    scan_number = peptide_details["scan_number"].values[0]
    scan = ms1_df[ms1_df['scanms1']==scan_number][["mz_values_ms1","intensity_values_ms1", "rtms1"]]

    # Extract lists from DataFrame
    x_values = scan['mz_values_ms1'].tolist()
    y_values = scan['intensity_values_ms1'].tolist()
    x_values = [item for sublist in x_values for item in sublist]
    y_values = [item for sublist in y_values for item in sublist]
   
    return x_values, y_values


##Getting peptide details from MQ files
# Locate the most abundant peptides containing R, K, and P
def get_most_abundant_krp_peptide(msms, evidence):
    msms = msms[msms['Proteins'].notna()]
    df_sorted = msms.sort_values(by='Precursor_Intensity', ascending=False)

    conditions = [
        df_sorted['Sequence'].str.contains('P'),
        (df_sorted['Sequence'].str.contains('K') & ~df_sorted['Sequence'].str.contains('P')),
        (df_sorted['Sequence'].str.contains('R') & ~df_sorted['Sequence'].str.contains('P'))
    ]

    peptides = []
    for condition in conditions:
        sequence = df_sorted[condition].iloc[0]
        df_sorted = df_sorted[df_sorted['Sequence'] != sequence['Sequence']]
        peptides.append(get_peptide_details(sequence, evidence))

    return peptides #list of peptide details

# Extract peptide details
def get_peptide_details(peptide, evidence):
    full_scan_number = peptide["Precursor_full_scan_number"]
    evidence_id = peptide['Evidence_ID']
    labeling_state = peptide['Labeling_state']
    
    peptide_details = evidence[evidence["id"] == evidence_id][['Raw_file','Sequence', 'm/z', 'Charge', 'Retention_time']]    
    peptide_details['scan_number'] = full_scan_number
    peptide_details["labeling_state"] = labeling_state
    
    sequence = peptide_details["Sequence"].values[0] 
    peptide_details["AA_mass"] = get_AA_mass(sequence)
    
    return peptide_details  

def get_AA_mass(sequence):
    if "K" in sequence:
        return 8
    elif "R" in sequence:
        return 10
    return 1

# Locate the most abundant heavy peptide
def get_most_abundant_heavy(msms, evidence):
    msms = msms[msms['Proteins'].notna()]
    msms = msms[msms["Labeling_state"]==1]
    df_sorted = msms.sort_values(by='Precursor_Intensity', ascending=False)
    largest_heavy = df_sorted.iloc[0]
    peptide = get_peptide_details(largest_heavy, evidence)
    return peptide
    
















