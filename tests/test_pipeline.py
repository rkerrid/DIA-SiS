from dia_sis.workflow.pipeline import Pipeline as pileline


if __name__ == "__main__":
    
 
    path = 'G:/My Drive/Data/main experiments/20240219 baby benchmark for pydia_sis/'
    pipeline = pileline( f'{path}', 'params.json', meta='meta.csv')
    # pipeline.make_metadata()
    pipeline.execute_pipeline()
  
    


