# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 13:00:02 2023

@author: rkerrid
"""

import pandas as pd
import numpy as np
from alphapept.pyrawfilereader import RawFileReader
import logging
from tqdm import tqdm
import os


def list_raw_files(directory_path):
    """
    Function to list all files ending with '.raw' in a directory
    
    Parameters:
    directory_path (str): Path to the directory
    
    Returns:
    list: List of filenames ending with '.raw'
    """
    raw_files = []  # Create an empty list to store the filenames
    
    # Iterate over the filenames in the directory
    for filename in os.listdir(directory_path):
        if filename.endswith(".raw"):  # Check if the filename ends with ".raw"
            raw_files.append(filename)  # Add the filename to the list
    print(f"found the following raw files for label check {raw_files}")
    return raw_files  # Return the list of filenames

##Import MQ
def import_file(file):
    df = pd.read_csv(file, sep='\t')
    all_cols = list(df.columns)
    df.columns = [x.replace(" ", "_") for x in all_cols ]
    return df

##Import thermo
# this code is either directly coppied or heavily based on on code from 
# the proteomics visualization package from the Mann lab https://github.com/MannLabs/ProteomicsVisualization and 
# Alphapept from the same lab https://github.com/MannLabs/alphapept

def get_raw_data(path, raw_files, no_of_peaks):
    ms1_dfs = {}

    for raw_file in raw_files:
        print(f"Begin processing {raw_file}")
        # Extracting filename without extension to use as dictionary key
        file_name = raw_file.split('.')[0]  
        ms1_dfs[file_name] = load_thermo_raw(path + raw_file, no_of_peaks)
    return ms1_dfs
    
def load_thermo_raw(
    raw_file,
    most_abundant=1000
    ):
    """
    Load a Thermo raw file and extract all spectra
    """
    
    rawfile = RawFileReader(raw_file)

    spec_indices = np.array(
        range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1)
        
    )
    print(spec_indices)
    scan_list = []
    rt_list = []
    mass_list = []
    int_list = []
    ms_list = []
    prec_mzs_list = []
    mono_mzs_list = []
    charge_list = []

    for i in tqdm(range(len(spec_indices))): 
        try:
            ms_order = rawfile.GetMSOrderForScanNum(i)
        
            rt = rawfile.RTFromScanNum(i)

            if ms_order == 2:
                prec_mz = rawfile.GetPrecursorMassForScanNum(i, 0)

                # mono_mz, charge = rawfile.GetMS2MonoMzAndChargeFromScanNum(i)
            else:
                prec_mz, mono_mz, charge = 0,0,0

            masses, intensity = rawfile.GetCentroidMassListFromScanNum(i)
          
            # if ms_order == 2:
            #     masses, intensity = get_most_abundant(masses, intensity, most_abundant)
        
            scan_list.append(i)
            rt_list.append(rt)
            mass_list.append(np.array(masses))
            int_list.append(np.array(intensity, dtype=np.int64))
            ms_list.append(ms_order)
            prec_mzs_list.append(prec_mz)
            mono_mzs_list.append(mono_mz)
            charge_list.append(charge)
        except KeyboardInterrupt as e:
            raise e
        except SystemExit as e:
            raise e
        except Exception as e:
            logging.info(f"Bad scan={i} in raw file '{raw_file}'")

    scan_list_ms1 = [scan_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    rt_list_ms1 = [rt_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    mass_list_ms1 = [mass_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    int_list_ms1 = [int_list[i] for i, _ in enumerate(ms_list) if _ == 1]
    ms_list_ms1 = [ms_list[i] for i, _ in enumerate(ms_list) if _ == 1]

    check_sanity(mass_list)

    data = {}
    
    data["scan_list_ms1"] = np.array(scan_list_ms1)
    data["rt_list_ms1"] = np.array(rt_list_ms1)
    data["mass_list_ms1"] = np.array(mass_list_ms1, dtype=object)
    data["int_list_ms1"] = np.array(int_list_ms1, dtype=object)
    data["ms_list_ms1"] = np.array(ms_list_ms1)

    rawfile.Close()
    
    ms1_df = pd.DataFrame(
        {'scanms1':data["scan_list_ms1"],
         'rtms1':data["rt_list_ms1"],
         'mz_values_ms1':data["mass_list_ms1"],
         'intensity_values_ms1':data['int_list_ms1']
            }
        ) 
    
    return ms1_df

def check_sanity(mass_list):
    """
    Sanity check for mass list to make sure the masses are sorted
    """
    if not all(
        mass_list[0][i] <= mass_list[0][i + 1] for i in range(len(mass_list[0]) - 1)
    ):
        raise ValueError("Masses are not sorted.")
        
        
        
