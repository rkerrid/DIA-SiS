# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 05:04:13 2023

@author: robbi
"""
import os

def create_directory(path, directory):
    # Combine the paths
    new_folder_path = os.path.join(f'{path}{directory}')
    
    # Create the new folder
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        print(f"Folder created successfully at {new_folder_path}")
    else:
        print(f"Folder already exists at {new_folder_path}")