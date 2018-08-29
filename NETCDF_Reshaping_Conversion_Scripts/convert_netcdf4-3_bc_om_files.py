# -*- coding: utf-8 -*-
"""

Simple script for converting NetCDF4 data files to NetCDF3 format
(for use with our fortran code, which uses only NetCDF3 format).
It also renames the files, to match the expected filenames.


@author: fqu12abu & douglowe
"""

import netCDF4 as nc
from netCDF4 import Dataset as ds



def convert_netcdf(file_in,file_out):
	#############(D)Convert from netCDF4 to netCDF3_Classic  #######
	#input file
	dsin = ds(file_in) 

	#output file
	dsout = ds(file_out, "w", format="NETCDF3_CLASSIC")

	#Copy dimensions
	for dname, the_dim in dsin.dimensions.iteritems():
		print dname, len(the_dim)
		dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


	# Copy variables
	for v_name, varin in dsin.variables.iteritems():
		outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
	  #  print varin.datatype
	
	# Copy variable attributes
		outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
	
		outVar[:] = varin[:]
	# close the output file
	dsout.close()

	#Check netCDF version
	print dsout.file_format


convert_netcdf('netcdf4/tno_bc1.nc',     'netcdf3/tno_ecoc_bc1.nc')
convert_netcdf('netcdf4/tno_ec_1_25.nc', 'netcdf3/tno_ecoc_ec_1_25.nc')
convert_netcdf('netcdf4/tno_ec_25_10.nc','netcdf3/tno_ecoc_ec_25_10.nc')
convert_netcdf('netcdf4/tno_oc_25.nc',   'netcdf3/tno_ecoc_oc_25.nc')
convert_netcdf('netcdf4/tno_oc_25_10.nc','netcdf3/tno_ecoc_oc_25_10.nc')
convert_netcdf('netcdf4/tno_pm_25.nc',   'netcdf3/tno_pm_25.nc')
convert_netcdf('netcdf4/tno_pm10.nc',    'netcdf3/tno_pm10.nc')
