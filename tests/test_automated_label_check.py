# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 16:51:01 2023

@author: rkerrid
"""
from dia_sis.label_check import file_io as io
from dia_sis.label_check import get_peptide_details 
from dia_sis.label_check import plot_spectra
from dia_sis.label_check import process_peptide_details

# work
path = 'C:/data/silac_dia_tools_files/data/spikein data/'
# home
path = 'G:/My Drive/Data/data/label check proof of concept/'
path = 'G:/My Drive/Data/data/label check test/'

# # path of raw files
raw_files = io.list_raw_files(path)

#import MQ and thermo raw file data
msms = io.import_file(f"{path}msms.txt")
evidence = io.import_file(f"{path}evidence.txt")
ms1_dfs = io.get_raw_data(path, raw_files, 1000)
# ic(ms1_dfs)

# get peptide details
peptides = get_peptide_details.get_most_abundant_krp_peptide(msms,evidence) 

# process details and save temp csv
process_peptide_details.process_details(peptides, ms1_dfs, path)
# import and plot data
plot_spectra.import_ms_data_for_plotting(path, peptides)