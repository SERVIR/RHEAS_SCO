# CODE TO CONVERT GTIFFS TO NC FILES 
# USINF CDO AND NCKS LIBRARIRES
#####################################

import os, sys, glob
import numpy as np
import datetime as dt

countries = ['ethiopia','kenya','malawi','rwanda','tanzania','uganda','zambia']



# LIST OF VARIABLES TO BE CONVERTED
#var = ['albedo','baseflow','cdi','dryspells','evap','latent','par',
#	'prec','rainf','rootmoist','runoff','sensible','severity',
#	'smdi','spi1','spi3','spi6','spi12','sri1','sri3','sri6','sri12','swe',
#	'surfstor','tmax','tmin']
var = ['albedo','baseflow','evap','latent','par',
	'prec','rainf','rootmoist','runoff','sensible',
	'swe','surfstor','tmax','tmin']
#var = ['evap']


# SET START AND END YEARS
iyear = 2017
eyear = 2022


for country in countries:
	# PATH TO WORKING DIR
	wdir = '/data/RHEAS/{0}/nc/25km/'.format(country)

	# LOOP THROUGH ALL VARIABLES
	for v in var:

		# FIRST CREATE SINGLE DAY NC FILE FROM GTIFF
		flist = glob.glob(wdir+v+'/tiffs/*.tiff')
		flist = list(np.sort(flist))

		# SET CURRENT WORKING DIR, IF DOESNOT EXISTS THEN CREATE ONE
		cwdir = wdir+v
		if not os.path.exists(cwdir):
			os.makedirs(cwdir)
		print(cwdir)
		# DELETE ANY OLD EXISTING FILE
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

				outname = cwdir+'/'+v+'_'+tmp2[4]+'.nc'
				print(outname)
				#sys.exit()
				cmd1 = 'gdal_translate -of NetCDF '+f+' '+outname
				os.system(cmd1)
			# LOOP 2 
			# GO THROUGH ALL THE NC FILES CREATED IN LOOP 1
			flist2 = glob.glob(cwdir+'/*.nc')
			for f in flist2:
				tmp1 = f.split('/')[-1].split('.')[0]
				tmp2 = tmp1.split('_')
				tmp3 = tmp2[1].split('-')
				ymd = tmp3[0]+tmp3[1]+tmp3[2]
				out = cwdir+'/tmp_'+ymd+'.nc'
				print(out)
				cmd2 = 'cdo -f nc4 -setcalendar,standard -settaxis,'+ymd+',00:00,1day '+f+' '+out
				os.system(cmd2)

			#LOOP 3 - GOING THROUGH ALL THE YEARS
			for year in range(iyear, eyear): 
				#merge different netcdfs into 1 using as time stack 
				cmd13 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_'+str(year)+'*.nc '+cwdir+'/tmp_'+str(year)+'_merged_1.nc'
				os.system(cmd13)

			cm1 = 'cdo -f nc4 mergetime '+cwdir+'/tmp_*_merged_1.nc '+cwdir+'/tmp_merged_1.nc'

			os.system(cm1) 

			#change attributes of the nc file
			os.system('ncrename -v Band1,'+v+' '+cwdir+'/tmp_merged_1.nc')
			os.system('ncatted -a long_name,'+v+',m,c,"'+v+'" '+cwdir+'/tmp_merged_1.nc')
			#os.system('ncatted -a units,SM,c,c,"m^3/m^3" '+cwdir+'_tmp_merged.nc')

			#strip crs and grid_mapping info
			os.system('ncatted -a grid_mapping,'+v+',d,, '+cwdir+'/tmp_merged_1.nc')

			os.system('ncks -O -v '+v+' '+cwdir+'/tmp_merged_1.nc '+cwdir+'/'+v+'_final.nc')
			#sys.exit()

			#delete temp file
			#os.system('rm '+cwdir+'*tmp*.nc')
			#os.system('rm '+cwdir+'tmp_*.nc')
			#'''
        

