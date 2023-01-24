# Author of this program -- Narendra N. Das (JPL)

#This python program is to perform ingestion 
# of Fertilizer for rice or maize
#edited for use on command line (no GUI)
#Import libraries

from rheas.dssat import utils as dutils
from datetime import datetime
import os

#define fertilizer dates and amounts
#format: {days_after_planting_1:amount1,days_after_planting_2:amount2}
#example: {0.0: 30.0, 15.0: 45.0, 30.0: 20.0} would be 30kg/ha at planting, 45 kg/ha 15 days after planting, and 20 kg/ha 30 days after planting
params = {0.0: 10.0, 30.0: 15.0}
#database name
db_name = 'rheasupdate'
#number of ensembles, shouldn't change from 1 as only the first is read
ensem_no = int(1)
#shapefile 
shp_filename = '/data/RHEAS/kenya/shp/gadm36_KEN_0.shp' #warning- if your shapefile has an attribute called 'name' it will cause an error
#name (only dynamic)
name='dynamic'
#crop - can be 'rice' or 'maize'
crop='maize'

def PopulateFertilizer():    

    now=datetime.now()
    dt=str(now)[:16]

    
    dutils.addfertilizer(db_name, shp_filename, params, ensem_no, dt,name,crop)
    print('Ingestion finished successfully!')
    

PopulateFertilizer()

