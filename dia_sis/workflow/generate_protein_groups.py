import pandas as pd
import numpy as np
import time
from .utils import manage_directories
from tqdm import tqdm


class DiaSis:
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        
        self.formatted_precursors = None
        self.protein_groups = None
        self.href_df = None
        
    def generate_protein_groups(self):
        start_time = time.time()
        # formatting channels
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        
        # calculate L/H ratios for each precursor       
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors)
        
        # Calculate global href df
        self.href_df = self.calculate_precursor_href_intensities_sum(self.formatted_precursors)
        
        # Compute protein level ratios
        self.protein_ratios = self.compute_protein_level(self.formatted_precursors)
        self.protein_ratios.to_csv(f'{self.path}/preprocessing/ratios.csv', sep=',')
        
        # Merge href df with protein level ratios to normalize light intensities
        self.protein_groups = self.href_normalization(self.protein_ratios, self.href_df) 
        
        # Output protein groups table of Light intensities and corrosponding href df
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
        formatted_precursors = pd.concat([pivot_L, pivot_H], axis=1)
        
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        formatted_precursors.reset_index(inplace=True)
        
        # Replace 0, inf, -inf with NaN for the specified columns
        formatted_precursors['Precursor.Translated H'] = formatted_precursors['Precursor.Translated H'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Precursor.Translated L'] = formatted_precursors['Precursor.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        
        formatted_precursors['Ms1.Translated H'] = formatted_precursors['Ms1.Translated H'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Ms1.Translated L'] = formatted_precursors['Ms1.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        
        return formatted_precursors

    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        df['Precursor.Translated L/H'] = df['Precursor.Translated L'] / df['Precursor.Translated H'] 
        df['Ms1.Translated L/H'] = df['Ms1.Translated L'] / df['Ms1.Translated H'] 
      
        return df

    # def calculate_precursor_href_intensities(self, df): # Older method where median H channel precursors were used to generate H reference
    #     print('Calculate href df')
    #     df = df.copy(deep = True)
    #     df = df.dropna(subset=['Precursor.Translated H','Ms1.Translated H'])
        
    #     def combined_median(ms1_series, precursor_series):
    #         combined_series = np.concatenate([ms1_series, precursor_series])
    #         combined_series = np.log10(combined_series)  # Log-transform the combined series
    #         return np.median(combined_series)  # Return the median of the log-transformed values
       
    #     # Group by protein group and apply the custom aggregation
    #     grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
    #         'href': combined_median(x['Ms1.Translated H'], x['Precursor.Translated H']) 
    #     })).reset_index()
       
    #     return grouped[['Protein.Group', 'href']]
    
    def calculate_precursor_href_intensities_sum(self, df): 
        print('Calculate href df')
        single_cell = True
        df = df.copy(deep = True)
            
        # df = df.dropna(subset=['Precursor.Translated H','Ms1.Translated H'])
        df = df.dropna(subset=['Ms1.Translated L/H'])
        runs = df['Run'].unique()
        runs_list = []
        for run in tqdm(runs, desc='Computing heavy intensities for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_sum(ms1_series, precursor_series):
                total_intensity = ms1_series + precursor_series
                
                total_intensity = total_intensity.dropna()
                
                total_intensity = np.sum(total_intensity)
                return np.log10(total_intensity)  # Return the median of the log-transformed values
           
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
                'href': combined_sum(x['Ms1.Translated H'], x['Precursor.Translated H']) 
            })).reset_index()
           
            grouped['Run'] = run
            runs_list.append(grouped)
        
        # combine runs into dataframe
        combined_df = pd.concat(runs_list, axis=0, ignore_index=True)

        # Group by 'Protein.Group' and calculate the median of 'href'
        result = combined_df.groupby('Protein.Group')['href'].median().reset_index()     
        return result
   
    def compute_protein_level(self, df):
        single_cell = True
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
        
       
        # df = df.dropna(subset=['Precursor.Translated L/H','Ms1.Translated L/H'])
        df = df.dropna(subset=['Ms1.Translated L/H'])
        
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_median_ratios(ms1_series, precursor_series):
                combined_series = np.concatenate([ms1_series, precursor_series])
                combined_series = combined_series[~np.isnan(combined_series)]
                combined_series = np.log10(combined_series)  # Log-transform the combined series
                return np.median(combined_series)  # Return the median of the log-transformed values
    
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby('Protein.Group').apply(lambda x: pd.Series({
                'L/H ratio': combined_median_ratios(x['Ms1.Translated L/H'], x['Precursor.Translated L/H'])
            })).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)
        cols = ['Run', 'Protein.Group', 'L/H ratio']
        return result[cols]

    
    def href_normalization(self, protein_groups, href):
        print('Calculating adjusted intensities using reference')
        # Merge the href_df onto protein groups containing protein level ratios
        merged_df = protein_groups.merge(href, on='Protein.Group', how='inner')
        
        # Obtain normalized light intensities by adding the L/H ratio to the heavy refference in log space
        merged_df['L_norm'] = merged_df['L/H ratio'] + merged_df['href']
        
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










