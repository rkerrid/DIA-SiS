# # -*- coding: utf-8 -*-
# """
# Created on Wed Nov 29 11:52:35 2023

# @author: rkerrid
# """

from icecream import ic
from dia_sis.pipeline.pipeline import Pipeline as pileline
import pandas as pd 

'''attempting to format the silac channels first then filter afterwards. Filter columns to keep in this step are:
    Parameters used:
        'Lib.PG.Q.Value'
Precursor.Charge > 1
Mass.Evidence > 0.5
Global.PG.Q.Value < 0.01
Channel.Q.Value < 0.03
Translated.Q.Value < 0.03
Translated.Quality >= 0.05

additional columns may be required


'''
if __name__ == "__main__":
    
 
   
    
    path = 'G:/My Drive/Data/main experiments/20240219 baby benchmark for pydia_sis/'
    # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    pipeline = pileline( f'{path}', 'params.json', meta='meta.csv')
    # pipeline.make_metadata()
    pipeline.execute_pipeline() # in href mode
  
    


