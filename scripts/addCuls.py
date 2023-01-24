# script for importing cultivars for multiple shapefiles at once

from rheas.dssat import utils


params = []
params.append({'p1':320, 'p2':0.52, 'p5':940, 'g2':620, 'g3':6.0, 'phint':38.9}) #example cultivar 1
params.append({'p1':320, 'p2':0.52, 'p5':940, 'g2':620, 'g3':6.0, 'phint':38.9}) #example cultivar 2


shps = []

dbname = 'rheas'


shps.append('kenya/shp/TransNzoia.shp') #example shapefile 1
shps.append('kenya/shp/UasinGishu.shp') #example shapefile 2

for i in range(len(shps)):
	print(shps[i])
	#print params[i]	
	utils.addCultivar(dbname, shps[i], params[i], nens=50, crop="maize")
