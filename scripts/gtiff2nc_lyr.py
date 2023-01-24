'''
this script converts a time series of tiffs into a netcdf
this works with layered variables such as soil moisture and soil temperature
some paths may need changed
'''

#from netCDF4 import Dataset, date2num, num2date
import os, sys, glob
import numpy as np
import datetime as dt
#import rasterio as rio


countries = ['ethiopia','kenya','malawi','rwanda','tanzania','uganda','zambia']




var = ['soil_moist','soil_temp']
iyear = 2017
eyear = 2022
#for yy in range(iyear,eyear+1):

for country in countries:
        # PATH TO WORKING DIR
	wdir = '/data/RHEAS/{0}/nc/25km/'.format(country)
	for v in var:
		flist = glob.glob(wdir+v+'/tiffs/*.tiff')
		flist = list(np.sort(flist))
		cwdir = wdir+v
		if not os.path.exists(cwdir):
			os.makedirs(cwdir)
		nctest = glob.glob(cwdir+'/*.nc')
		if len(nctest)>0:
			print('Deleting old files...')
			os.system('rm '+cwdir+'/*.nc')
	
		if len(flist)>0:
			print('Current variable: ',v)
						
			# create temp nc file from gtiff
			for f in flist:
				tmp1 = f.split('/')[-1].split('.')[0]
				tmp2 = tmp1.split('_')
				lyr = tmp2[6]
				dt = tmp2[5]
				outname = cwdir+'/'+v+'_'+dt+'_'+lyr+'.nc'
				print(outname)
				cmd1 = 'gdal_translate -of NetCDF '+f+' '+outname
				os.system(cmd1)
				#sys.exit()
		
			#cwdir = wdir+'nc/'+v
			flist2 = glob.glob(cwdir+'/*.nc')
			#sys.exit()
			for f in flist2:
				tmp1 = f.split('/')[-1].split('.')[0]
				tmp2 = tmp1.split('_')
				lyr = tmp2[3]
				yyyy = tmp2[2][:4]
				mm = tmp2[2][4:6]
				dd = tmp2[2][6:]
				ymd = tmp2[2]
				out = cwdir+'/tmp_'+ymd+'_'+lyr+'.nc'
				#print(tmp1,tmp2)
				print(out)
				#sys.exit()
				cmd2 = 'cdo -f nc4 -setcalendar,standard -settaxis,'+ymd+',00:00,1day '+f+' '+out
				os.system(cmd2)
		
			for year in range(iyear, eyear+1): 
				#merge different netcdfs into 1 using as time stack 
				cmd13 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_'+str(year)+'*_1.nc '+cwdir+'/tmp_'+str(year)+'_merged_1.nc'
				cmd23 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_'+str(year)+'*_2.nc '+cwdir+'/tmp_'+str(year)+'_merged_2.nc'
				cmd33 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_'+str(year)+'*_3.nc '+cwdir+'/tmp_'+str(year)+'_merged_3.nc'
				os.system(cmd13)
				os.system(cmd23)
				os.system(cmd33)
			cm1 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_20*_merged_1.nc '+cwdir+'/tmp_merged_1.nc'
			cm2 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_20*_merged_2.nc '+cwdir+'/tmp_merged_2.nc'
			cm3 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_20*_merged_3.nc '+cwdir+'/tmp_merged_3.nc'
			os.system(cm1)
			os.system(cm2)
			os.system(cm3) 

			#change attributes of the nc file
			os.system('ncrename -v Band1,'+v+' '+cwdir+'/tmp_merged_1.nc')
			os.system('ncrename -v Band1,'+v+' '+cwdir+'/tmp_merged_2.nc')
			os.system('ncrename -v Band1,'+v+' '+cwdir+'/tmp_merged_3.nc')
			os.system('ncatted -a long_name,'+v+',m,c,"'+v+'" '+cwdir+'/tmp_merged_1.nc')
			os.system('ncatted -a long_name,'+v+',m,c,"'+v+'" '+cwdir+'/tmp_merged_2.nc')
			os.system('ncatted -a long_name,'+v+',m,c,"'+v+'" '+cwdir+'/tmp_merged_3.nc')
			#os.system('ncatted -a units,SM,c,c,"m^3/m^3" '+cwdir+'_tmp_merged.nc')

			#strip crs and grid_mapping info
			os.system('ncatted -a grid_mapping,'+v+',d,, '+cwdir+'/tmp_merged_1.nc')
			os.system('ncatted -a grid_mapping,'+v+',d,, '+cwdir+'/tmp_merged_2.nc')
			os.system('ncatted -a grid_mapping,'+v+',d,, '+cwdir+'/tmp_merged_3.nc')

			os.system('ncks -O -v '+v+' '+cwdir+'/tmp_merged_1.nc '+cwdir+'/'+v+'_final.nc')
			os.system('ncks -O -v '+v+' '+cwdir+'/tmp_merged_2.nc '+cwdir+'/'+v+'_final_2.nc')
			os.system('ncks -O -v '+v+' '+cwdir+'/tmp_merged_3.nc '+cwdir+'/'+v+'_final_3.nc')
			#sys.exit()
	
			#delete temp file
			#os.system('rm '+cwdir+'*tmp*.nc')
			#os.system('rm '+cwdir+'tmp_*.nc')
			#'''
        

