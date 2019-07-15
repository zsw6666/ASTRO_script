import numpy as np
import pandas as pd

def Read_dat(pathname,sep='\s+',engine='python',header=None):

    datafile=pd.read_csv(pathname,sep,engine=engine,header=header)

    return datafile.values()