# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 17:19:53 2023

@author: robbi
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from icecream import ic
# # Load the datasets
light_data_path = 'G:/My Drive/Data/data/testing pipeline dev/bm whole set/new/protein intensities/light_href.csv'
nsp_data_path = 'G:/My Drive/Data/data/testing pipeline dev/bm whole set/new/protein intensities/reference_href.csv'

light_data_path = 'G:/My Drive/Data/data/testing pipeline dev/bm whole set/old/protein intensities/light_href.csv'
nsp_data_path = 'G:/My Drive/Data/data/testing pipeline dev/bm whole set/old/protein intensities/reference_href.csv'

light_data = pd.read_csv(light_data_path)
nsp_data = pd.read_csv(nsp_data_path)
# ic(nsp_data)
# # Merging the datasets on the 'Protein.Group' column
merged_data = pd.merge(light_data, nsp_data, on='Protein.Group', suffixes=('_light', '_nsp'))
# Filtering to keep only rows where 'Protein.Groups' contains 'ECOLI_'
merged_data = merged_data[merged_data['Protein.Group'].str.contains('ECOLI_')]

# Now, 'filtered_data' contains only the rows where 'Protein.Groups' has 'ECOLI_'


# Filtering out columns for sample S3_1
s3_1_light = merged_data['S3_1_light']
s3_1_nsp = merged_data['S3_1_nsp']

# # Avoiding division by zero and excluding missing values
# valid_indices = (s3_1_nsp != 0) & ~s3_1_nsp.isna() & ~s3_1_light.isna()
# s3_1_light = s3_1_light[valid_indices]
# s3_1_nsp = s3_1_nsp[valid_indices]

# Calculating the ratio
ratio_s3_1 = s3_1_light / s3_1_nsp
ic(ratio_s3_1)

# Creating the scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(np.log2(ratio_s3_1), np.log2(s3_1_light), alpha=0.5)
plt.title('Intensity vs. Ratio (Light/Nsp)')
plt.xlabel('Ratio (Light/Nsp)')
plt.ylabel('Intensity (Light)')
plt.grid(True)
plt.show()

# Extracting the first replicate (Rep 1) from each sample
samples = [col.split('_')[0] for col in light_data.columns if '_' in col]  # Extract sample names
samples = list(set(samples))  # Remove duplicates

# Initialize lists to store ratios and non-zero counts
ratios = {}
non_zero_counts = {}

# Loop through each sample to calculate the ratios and non-zero counts
for sample in samples:
    sample_light = merged_data[f'{sample}_1_light']
    sample_nsp = merged_data[f'{sample}_1_nsp']

    # Ratio calculation
    valid_indices = (sample_nsp != 0) & ~sample_nsp.isna() & ~sample_light.isna()
    ratio = sample_light[valid_indices] / sample_nsp[valid_indices]
    ratios[sample] = ratio

    # Counting non-zero values
    non_zero_counts[sample] = sample_light[sample_light != 0].count()

expected_5 = np.log2(100/10)
expected_4 = np.log2(10/10)
expected_3 = np.log2(2/10)
expected_2 = np.log2(1/10)
hline_values = [expected_5, expected_4, expected_3, expected_2]
# Convert the ratio dictionary to a DataFrame for the box plot
ratios_df = pd.DataFrame.from_dict(ratios, orient='index').transpose()
ic(ratios_df)
ratios_df = np.log2(ratios_df)
ordered_columns = ['S5', 'S4', 'S3', 'S2', 'S1']
ratios_df = ratios_df[ordered_columns]
ratios_df.boxplot(figsize=(12, 6))
plt.title('Box Plot of Ratios Between Light and Nsp Samples (Rep 1)')
plt.ylabel('Ratio log2(L/H)')
plt.xticks(rotation=45)
for value in hline_values:
    plt.hlines(value, xmin=-0.5, xmax=len(ordered_columns), colors='red', linestyles='dashed')
plt.grid(True)
plt.show()

# Bar plot of non-zero counts

# Ordering the samples for the bar plot
ordered_samples = ['S5', 'S4', 'S3', 'S2', 'S1']
ordered_non_zero_counts = [non_zero_counts[sample] for sample in ordered_samples]

# Bar plot of non-zero counts in the specified order
plt.figure(figsize=(12, 6))
plt.bar(ordered_samples, ordered_non_zero_counts, color='skyblue')
plt.title('Non-Zero Values Count in Light Dataset for Each Sample (Rep 1)')
plt.ylabel('Count of Non-Zero Values')
plt.xlabel('Sample')
plt.xticks(rotation=45)
plt.grid(axis='y')
plt.show()



