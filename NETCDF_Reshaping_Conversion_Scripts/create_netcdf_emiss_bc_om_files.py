# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 12:55:11 2016
Modified on Sat Feb 11 2018

This script is for creating NetCDF files, with emissions from the 2011 version 
of the TNO_MACC inventory, which can be read by the UoM TNO emissions preprocessor.
 
The script will produce the netcdf files for emissions of:
bc1
ec_1_25
ec_25_10
oc_25
oc_25_10

These will be calculated using the fractions of PM2.5 and PM10 mass 
for oc_fine, oc_coarse, bc_fine and bc_coarse that are provided by TNO. 
And will group the emissions by sector
(pow, res, inc, pei, exf, sol, tra1-2-3-4-5, nrt, was, and agr).

We will assume that the oc_25_10 and ec_25_10 mass fractions will be the
fraction of the (PM10 - PM2.5) mass (*not* simply of the PM10 mass, which
includes the PM2.5, and so would include bc1, ec_1_25, and ec_25_10).

Also, we will assume that oc_25 and oc_25_10 mass is the *OM* (not OC only)
mass. It will have to be reduced by a scaling factor to get an OC only
mass if that usage is required.

The splits between bc1 and ec_1_25 for each emission sector have been (roughly!)
estimated from the TNO 2005 data that we had previously. The mass fractions
of bc1 for each sector will be:
pow  = 0.6
res  = 0.75
inc  = 0.4
pei  = 0.4
exf  = 0.4
sol  =     0.0 (not used)
tra1 = 0.8
tra2 = 0.8
tra3 = 0.8
tra4 = 0.8
tra5 = 0.8
nrt  = 0.8
was  = 0.95
agr  = 0.95
(the mass fraction values of ec_1_25 will be (1 - bc1))


@author: fqu12abu & douglowe
"""

##Libraries used
import netCDF4 as nc
from netCDF4 import Dataset as ds
import numpy as np
import xray
#from monthdelta import monthdelta
import calendar
from datetime import datetime
from netCDF4 import num2date, date2num

import pandas as pd


## alternative to monthdelta, from https://stackoverflow.com/questions/4130922/how-to-increment-datetime-by-custom-months-in-python-without-using-library
def add_months(sourcedate,months):
	month = sourcedate.month - 1 + months
	year = sourcedate.year + month // 12
	month = month % 12 + 1
	day = min(sourcedate.day,calendar.monthrange(year,month)[1])
	return datetime(year,month,day)

## initialise the dates array
dates = []

## initialise the sector list, BC1 fractionation information, and SNAP value
sector_list = ('pow','res','inc','exf','sol','tra1','tra2','tra3','tra4','tra5','nrt','was','agr')

sect_bc1_frac = {'pow': 0.6, 'res': 0.75, 'inc': 0.4, 'pei': 0.4,
				 'exf': 0.4, 'sol': 0.0, 'tra1': 0.8, 'tra2': 0.8,
				 'tra3':0.8, 'tra4':0.8, 'tra5': 0.8, 'nrt': 0.8,
				 'was':0.95, 'agr':0.95 }
				 
sect_snap = {'pow': 1, 'res':  2, 'inc': 34,  'pei': -999,
			'exf':  5, 'sol':  6, 'tra1': 71, 'tra2': 72,
			'tra3':73, 'tra4':74, 'tra5': 75, 'nrt': 8,
			'was':  9, 'agr': 10 }



############(A) Open original TNO emissions inventory################
tno_2011 = './TNO_MACC_III_emissions_2011.nc4'  #input file
tno = nc.Dataset(tno_2011,'r')


############(A) Open the PM split excel spreadsheet ################
xl = pd.ExcelFile('PM_split_for_TNO_MACC_III.xlsx')
#### pull out the fine and coarse datasheets
fine_data = xl.parse('fine')
coarse_data = xl.parse('coarse')

#### set the index order - this makes searching through the database a *lot* faster later on...
fdi = fine_data.set_index(['Year','ISO3','SNAP'])
cdi = coarse_data.set_index(['Year','ISO3','SNAP'])

#### old method for attempting to speed things up - without much success
#fine_data_2011 = fine_data.loc[fine_data['Year'] == 2011]
#coarse_data_2011 = coarse_data.loc[coarse_data['Year'] == 2011]



############(B)Function for creating new netcdf file ################

def create_new_netcdf_file(file):
	dataset = ds(file,'w','r', format='NETCDF4_CLASSIC')
	
	#2)Create the dimensions
	lat = dataset.createDimension('lat', 672) #number of latitudes
	lon = dataset.createDimension('lon', 720) #number of longitudes
	time = dataset.createDimension('time', None)
	
	#3)Create variables
	######1D variables first######
	
	lat = dataset.createVariable('lat',np.float32, ('lat'), fill_value=False)
	lon = dataset.createVariable('lon',np.float32, ('lon'), fill_value=False)
	time = dataset.createVariable('time',np.float32, ('time'), fill_value=False)
	
	######3D variables (SNAP sectors) ########
	
	pow = dataset.createVariable('pow',np.float64, ('time','lat','lon'), fill_value=False)
	res = dataset.createVariable('res',np.float32, ('time','lat','lon'), fill_value=False)
	inc = dataset.createVariable('inc',np.float32, ('time','lat','lon'), fill_value=False)
	pei = dataset.createVariable('pei',np.float32, ('time','lat','lon'), fill_value=False)
	exf = dataset.createVariable('exf',np.float32, ('time','lat','lon'), fill_value=False)
	sol = dataset.createVariable('sol',np.float32, ('time','lat','lon'), fill_value=False)
	tra1 = dataset.createVariable('tra1',np.float32, ('time','lat','lon'), fill_value=False)
	tra2 = dataset.createVariable('tra2',np.float32, ('time','lat','lon'), fill_value=False)
	tra3 = dataset.createVariable('tra3',np.float32, ('time','lat','lon'), fill_value=False)
	tra4 = dataset.createVariable('tra4',np.float32, ('time','lat','lon'), fill_value=False)
	tra5 = dataset.createVariable('tra5',np.float32, ('time','lat','lon'), fill_value=False)
	nrt = dataset.createVariable('nrt',np.float32, ('time','lat','lon'), fill_value=False)
	was = dataset.createVariable('was',np.float32, ('time','lat','lon'), fill_value=False)
	agr = dataset.createVariable('agr',np.float32, ('time','lat','lon'), fill_value=False)
	
	
	#4) Add attributes units to 1D and 3D variables
	lat.units = 'degrees_north'
	lat.long_name = 'latitude'
	lon.units = 'degrees_east'
	lon.long_name = 'longitude'
	time.units = 'days since 1900-01-01 00:00'
	time.calendar ='gregorian'
	time.long_name = 'Time'
	
	pow.units = 'Kg yr-1'
	pow.long_name = 'Power generation'
	res.units = 'Kg yr-1'
	res.long_name= 'Residential, comercial and other combustion'
	inc.units = 'Kg yr-1'
	inc.long_name = 'Industrial combustion'
	pei.units = 'Kg yr-1'
	pei.long_name= 'Processed emission industrial'
	exf.units = 'Kg yr-1'
	exf.long_name = 'Extraction and distribution of fossil fuels'
	sol.units = 'Kg yr-1'
	sol.long_name = 'Solvent use'
	tra1.units = 'Kg yr-1'
	tra1.long_name = 'Road transport, gasoline'
	tra2.units = 'Kg yr-1'
	tra2.long_name = 'Road transport, diesel' 
	tra3.units = 'Kg yr-1'
	tra3.long_name = 'Road trasnport, LPG'
	tra4.units = 'Kg yr-1'
	tra4.long_name = 'Road trasnport, non-exhaust, volatilisation'
	tra5.units = 'Kg yr-1'
	tra5.long_name = 'Road transport, non-exhaust, wear'
	nrt.units = 'Kg yr-1'
	nrt.long_name = 'Non-road transport'
	was.units = 'Kg yr-1'
	was.long_name = 'Waste tratment and disposal'
	agr.units = 'Kg yr-1'
	agr.long_name = 'Agriculture'
	
	return dataset



############(B_pt2) Create new netcdf files with function above ######################
ds_bc1   = create_new_netcdf_file('netcdf4/tno_bc1.nc')
ds_ec25  = create_new_netcdf_file('netcdf4/tno_ec_1_25.nc')
ds_ec10  = create_new_netcdf_file('netcdf4/tno_ec_25_10.nc')
ds_oc25  = create_new_netcdf_file('netcdf4/tno_oc_25.nc')
ds_oc10  = create_new_netcdf_file('netcdf4/tno_oc_25_10.nc')
ds_pm25  = create_new_netcdf_file('netcdf4/tno_pm_25.nc')
ds_pm10  = create_new_netcdf_file('netcdf4/tno_pm10.nc')




 
                          
############# (C) extracting PM data from the original file, and organising it for the new files #######

#1)Extract and add latitude  and longitude values to our files

data = xray.open_dataset(tno_2011)  #Open NetCDF file (TNO-MACCII 2011) with x-ray library
print data

latitude_tno = tno['latitude'][:]  # Extract latitude from the original file
longitude_tno = tno['longitude'][:] 


def append_lat_lon_data(lat_tno,long_tno,file):
	file['lat'][:] = lat_tno
	file['lon'][:] = long_tno   

# append to the data files
append_lat_lon_data(latitude_tno,longitude_tno,ds_bc1)
append_lat_lon_data(latitude_tno,longitude_tno,ds_ec25)
append_lat_lon_data(latitude_tno,longitude_tno,ds_ec10)
append_lat_lon_data(latitude_tno,longitude_tno,ds_oc25)
append_lat_lon_data(latitude_tno,longitude_tno,ds_oc10)
append_lat_lon_data(latitude_tno,longitude_tno,ds_pm25)
append_lat_lon_data(latitude_tno,longitude_tno,ds_pm10)



#2)Extract PM2_5 and PM10 using the NETCDF4 library########
tno_net = ds(tno_2011) 
emis_cat_index = tno_net.variables['emission_category_index'][:]
lat_index =tno_net.variables['latitude_index'][:]
lon_index =tno_net.variables['longitude_index'][:]
pm25_data = tno_net.variables['pm2_5'][:]  
pm10_data = tno_net.variables['pm10'][:] 
### get the country ID info from NETCDF4 library
country_index = tno_net.variables['country_index'][:]
country_id    = tno_net.variables['country_id'][:]


#3)Loop through every emission category (emiss_cat) and pick the emissions values by sector, latitude, and longitude 
##### create a list of the sectors
## NOTE: this must be the same order as the sectors in the source file!!!!!!!

##### i) Create a 3D array with zeros 
pm25_arrays = np.zeros(shape=(13,12,672,720))
pm10_arrays = np.zeros(shape=(13,12,672,720))
bc_fine_arrays = np.zeros(shape=(13,12,672,720))
bc_coarse_arrays = np.zeros(shape=(13,12,672,720))
oc_fine_arrays = np.zeros(shape=(13,12,672,720))
oc_coarse_arrays = np.zeros(shape=(13,12,672,720))



####ii)Pick values from each sector and fill the 3D arrays with the new emission values
for i in range(len(lat_index)):
	
	if(i%10000 == 0): print 100.0*i/len(lat_index)
	
	sector = emis_cat_index[i]
	
	# save the sector data for PM25 and PM10
	pm25_arrays[sector-1,:,lat_index[i]-1, lon_index[i]-1] += pm25_data[i]
	pm10_arrays[sector-1,:,lat_index[i]-1, lon_index[i]-1] += pm10_data[i]
	

####iii) loop through the domain and sectors, building up arrays of scaling factors for EC and OC
for i in range(len(lat_index)):
	
	if(i%10000 == 0): print 100.0*i/len(lat_index)
	
	sector = emis_cat_index[i]
	c_id = ''.join(country_id[country_index[i]-1])
	
	# identify the fractional data for fine and coarse modes - using the indexes we set above
	temp_fine   = fdi.loc(axis=0)[2011,c_id,sect_snap[sector_list[sector-1]]]
	temp_coarse = cdi.loc(axis=0)[2011,c_id,sect_snap[sector_list[sector-1]]]
	
	### old method for searching for fractional data - very slow!
	#temp_fine = fine_data_2011.loc[(fine_data_2011['ISO3'] == c_id) &
	#							(fine_data_2011['SNAP'] == sect_snap[sector_list[sector-1]])]
	#temp_coarse = coarse_data_2011.loc[(coarse_data_2011['ISO3'] == c_id) &
	#								(coarse_data_2011['SNAP'] == sect_snap[sector_list[sector-1]])]
	
	bc_fine_arrays[sector-1,:,lat_index[i]-1, lon_index[i]-1] = temp_fine['EC_fine']
	oc_fine_arrays[sector-1,:,lat_index[i]-1, lon_index[i]-1] = temp_fine['OC_fine']
	bc_coarse_arrays[sector-1,:,lat_index[i]-1, lon_index[i]-1] = temp_coarse['EC_coarse']
	oc_coarse_arrays[sector-1,:,lat_index[i]-1, lon_index[i]-1] = temp_coarse['OC_coarse']



#####3) mapping data to BC and OC arrays #########


#### loop through the sectors, doing the math and saving the data to the files
for i in range(len(sector_list)):
	sector = sector_list[i]
	ds_bc1[sector][:] = pm25_arrays[i,:,:,:] * bc_fine_arrays[i,:,:,:] * sect_bc1_frac[sector]
	ds_ec25[sector][:] = pm25_arrays[i,:,:,:] * bc_fine_arrays[i,:,:,:] * (1.0 - sect_bc1_frac[sector])
	ds_oc25[sector][:] = pm25_arrays[i,:,:,:] * oc_fine_arrays[i,:,:,:]
	ds_ec10[sector][:] = (pm10_arrays[i,:,:,:] - pm25_arrays[i,:,:,:]) * bc_coarse_arrays[i,:,:,:]
	ds_oc10[sector][:] = (pm10_arrays[i,:,:,:] - pm25_arrays[i,:,:,:]) * oc_coarse_arrays[i,:,:,:]

for i in range(len(sector_list)):
	sector = sector_list[i]
	ds_pm25[sector][:] = pm25_arrays[i,:,:,:]
	ds_pm10[sector][:] = pm10_arrays[i,:,:,:]

### set all data in the "pei" sector to 0.0 (as we do not use that sector)
pei_zero_array = np.zeros(shape=(12,672,720))

ds_bc1['pei'][:] = pei_zero_array
ds_ec25['pei'][:] = pei_zero_array
ds_oc25['pei'][:] = pei_zero_array
ds_ec10['pei'][:] = pei_zero_array
ds_oc10['pei'][:] = pei_zero_array
ds_pm25['pei'][:] = pei_zero_array
ds_pm10['pei'][:] = pei_zero_array


####iii)Append values to the variable time (12 months)
#dates = [datetime(2000,01,01)+n*monthdelta(1) for n in range(nox_pow.shape[0])]     
dates = [add_months(datetime(2000,01,01),n) for n in range(12)]

def add_dates_to_files(dates,file):
	file['time'][:] = date2num(dates,units=file['time'].units,calendar=file['time'].calendar)

add_dates_to_files(dates,ds_bc1)
add_dates_to_files(dates,ds_ec25)
add_dates_to_files(dates,ds_ec10)
add_dates_to_files(dates,ds_oc25)
add_dates_to_files(dates,ds_oc10)
add_dates_to_files(dates,ds_pm25)
add_dates_to_files(dates,ds_pm10)



ds_bc1.close()
ds_ec25.close()
ds_ec10.close()
ds_oc25.close()
ds_oc10.close()

ds_pm25.close()
ds_pm10.close()











