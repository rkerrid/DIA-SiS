"""
Created on Sun Aug 13 13:13:37 2023

@author: rkerrid
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# import files to plot
def import_ms_data_for_plotting(path, peptides):    
    for peptide_details in peptides:
        sequence = peptide_details["Sequence"].values[0]
        df = pd.read_csv(f"{path}{sequence}.csv")
        mz = df.mz
        intensity = df.intensity
        #plot each peptide
        scan_plot(mz, intensity, peptide_details)

def plot_points_and_lines(mz_values, intensity_values, expected_light_peak, expected_heavy_peak):
    for x, y in zip(mz_values, intensity_values):
        color = 'blue' if (np.isclose(x, expected_light_peak, atol=1e-2) or np.isclose(x, expected_heavy_peak, atol=1e-2)) else 'red'
        plt.scatter([x], [y], color=color, marker='o', s=10)
        plt.vlines(x=x, ymin=0, ymax=y, colors=color, linestyles='dashed', linewidths=0.8)


def set_axes_limits(mz_values, intensity_values, expected_heavy_peak):
    plt.xlim(expected_heavy_peak - 10, expected_heavy_peak + 10)
    max_intensity_in_range = np.max(
        intensity_values[(mz_values >= expected_heavy_peak - 10) & (mz_values <= expected_heavy_peak + 10)])
    y_axis_margin = max_intensity_in_range * 0.1
    plt.ylim(0, max_intensity_in_range + y_axis_margin)


def draw_error_bar(expected_light_peak, expected_heavy_peak):
    y_limit_quarter = (plt.ylim()[1] - plt.ylim()[0]) / 1.4  # Adjusted the position
    error_diff = expected_heavy_peak - expected_light_peak
    plt.errorbar([expected_light_peak, expected_heavy_peak],
                 [y_limit_quarter, y_limit_quarter],
                 xerr=[0, 0], capsize=5, color='gray', linewidth=0.5)  # Faint color and line
    plt.annotate(f'{error_diff}',
                 ((expected_light_peak + expected_heavy_peak) / 2, y_limit_quarter),
                 textcoords="offset points", xytext=(0,5), ha='center', color='gray')  # Labeling the error


def annotate_peaks(mz_values, intensity_values, expected_light_peak, expected_heavy_peak):
    for x, y in zip(mz_values, intensity_values):
        if np.isclose(x, expected_light_peak, atol=1e-2):
            plt.annotate(f'{expected_light_peak:.2f} L',
                         (expected_light_peak, y),
                         textcoords="offset points", xytext=(0, 5), ha='center')
        elif np.isclose(x, expected_heavy_peak, atol=1e-2):
            plt.annotate(f'{expected_heavy_peak:.2f} H', (expected_heavy_peak, y),
                         textcoords="offset points", xytext=(0, 5), ha='center')


def scan_plot(mz_values, intensity_values, peptide_details):
    if peptide_details.empty:
        print("No peptide details found. Cannot create the plot.")
        return

    # Extract details...
    sequence = peptide_details["Sequence"].values[0]
    retention_time = peptide_details["Retention_time"].values[0]
    charge = peptide_details["Charge"].values[0]
    expected_heavy_peak = peptide_details["m/z"].values[0]
    labeling_state = peptide_details["labeling_state"].values[0]
    AA_mass = peptide_details["AA_mass"].values[0]
    scan_number = peptide_details["scan_number"].values[0]
    raw_file = peptide_details['Raw_file'].item()
    if labeling_state == 1:
        expected_heavy_peak = expected_heavy_peak + (AA_mass / charge)

    expected_light_peak = expected_heavy_peak - (AA_mass / charge)

    plot_points_and_lines(mz_values, intensity_values, expected_light_peak, expected_heavy_peak)
    set_axes_limits(mz_values, intensity_values, expected_heavy_peak)
    # y_limit_half = (plt.ylim()[1] - plt.ylim()[0]) / 1.1
    draw_error_bar(expected_light_peak, expected_heavy_peak)#y_limit_half
    annotate_peaks(mz_values, intensity_values, expected_light_peak, expected_heavy_peak)

    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.title(f'MS1 of {sequence} Charge: {charge} \n Raw file: {raw_file} \n RT: {retention_time}, Scan no: {scan_number}')
    plt.grid(True)
    plt.legend()
    plt.show()
