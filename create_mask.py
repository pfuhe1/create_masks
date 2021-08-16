# Script to take polygon(s) and turn it(them) into a mask on a particular lat-lon grid
# Peter Uhe
# 16 / 06 /2017
import numpy as np
from matplotlib.path import Path
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import unicodedata

# Use osgeo (gdal library) to load shapefile
# Installed by '$conda install gdal'
from osgeo import ogr

####################################################################################

# Load a text file containing lists of vertices (space separated, one vertex per line)
# Creates matplotlib.path.Path objects representing the polygons. 
# Lines starting with '#' are ignored
# A blank line indicates the end of a polygon, so multiple polygons can be defined with new lines in between
def load_polygons(fname):
	polygons=[]
	tmp=[]
	# Import Polygons of mask in fname (separate polygons are separated by a blank line)
	for line in open(fname,'r'):
		if line[0]=='#': 
			# skip commented lines
			continue
		elif line.strip()!='':
			# Add coordinates to list
			tuple=line.strip().split()
			tmp.append(tuple) # lon,lat
		else: # blank line
			# Finished reading data for this polygon
			# create polygon path out of list of vertices
			if tmp !=[]:
				polygons.append(Path(np.array(tmp)))
				tmp=[]
	# If the file didn't end in a blank line, add the final polygon
	if tmp !=[]:
		polygons.append(Path(np.array(tmp)))
	return polygons

###################################################################################

def load_grid(fname,latname='latitude0',lonname='longitude0'):
	# load region grid and returns list of points (lon,lat)
	with Dataset(fname,'r') as f:
		# Load 1d arrays of lon and lat
		lat=f.variables[latname][:]
		lon=f.variables[lonname][:]
		
		if len(lat.shape)==2:
			# 2D lat and lon:
			lonxx=lon
			latyy=lat
		else:
			# Create 2D arrays of lon and lat
			lonxx,latyy=np.meshgrid(lon,lat)
	return lonxx,latyy

##################################################################################

# Function to create mask given polygons object and points array
def create_mask(polygons,points,nlat,nlon):
	# Convert polygons to mask (true if inside the polygon region)
	# add the masks for multiple polygons together
	for i,polygon in enumerate(polygons):
		# Determine if  points inside polygon
		tmp_mask = polygon.contains_points(points)
		# Reshape mask to dimensions of the grid
		mask=np.reshape(tmp_mask,[nlat,nlon])

	return ~mask # Invert the mask so true is outside the region

#################################################################################

# Wrapper function to return the mask givent the polygons file and grid file
def load_and_create_mask(f_polygons,f_grid,latname='latitude0',lonname='longitude0'):
	# Load inputs and create mask
	polygons=load_polygons(f_polygons)
	# Load 2D lon and lat arrays for grid
	lonxx,latyy=load_grid(f_grid,latname,lonname)
	nlat,nlon=lonxx.shape
	# Stack points into a N x 2 array (where N = nlat x nlon)
	points = np.vstack((lonxx.flatten(),latyy.flatten())).T
	# Call create_mask function for polygons and grid points
	return create_mask(polygons,points,nlat,nlon)

#################################################################################

def add_to_text(fileh,polygon):
	for coord in polygon.vertices:
		fileh.write(str(coord[0])+' '+str(coord[1])+'\n')
	fileh.write('\n')

#################################################################################
# 
# Function to read in a shapefile and output the polygons of each feature
# input 'fieldname' needs to be a unique identifier for each polygon e.g. ID or NAME
def load_shapefile(shapefile,fieldname,field_list=None):
	
	print('Loading Shapefile')
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(shapefile, 0)
	layer = dataSource.GetLayer()
	counties={}
	boundaries=[]
	types=[]
	for feature in layer:
		try:
			region=feature.GetField(fieldname)
		except ValueError:
			print(feature.items())
			raise Exception('Error, field "'+fieldname+'" does not exist in the shapefile')
		print(region)
		if field_list is not None and region not in field_list:
			# Skip this region
			continue


		geometry=feature.GetGeometryRef()
		boundary=eval(feature.geometry().Boundary().ExportToJson())
		#geometry = json['geometry']

		if boundary['type']=='LineString':
			polygons=[Path(np.array(boundary['coordinates']))]
		elif boundary['type']=='MultiLineString':
			polygons=[]
			for p in boundary['coordinates']:
				polygons.append(Path(np.array(p)))
		else:
			print('Error: unknown geometry')
			continue

		if region is not None and boundary is not None:
			counties[region]=polygons
	return counties
	
#################################################################################
# 
# Function to read in a shapefile and output the polygons of each feature
# input 'fieldname' needs to be a unique identifier for each polygon e.g. ID or NAME
# Also returns attributes of each feature
def load_shapefile_attrs(shapefile,fieldname,field_list=None):
	
	print('Loading Shapefile')
	driver = ogr.GetDriverByName("ESRI Shapefile")
	dataSource = driver.Open(shapefile, 0)
	layer = dataSource.GetLayer()
	counties={}
	attributes = {}
	boundaries=[]
	types=[]
	for feature in layer:
		try:
			region=feature.GetField(fieldname)
		except ValueError:
			print(feature.items())
			raise Exception('Error, field "'+fieldname+'" does not exist in the shapefile')
		print(region)
		attributes[region] = feature.items()
		if field_list is not None and region not in field_list:
			# Skip this region
			continue


		geometry=feature.GetGeometryRef()
		boundary=eval(feature.geometry().Boundary().ExportToJson())
		#geometry = json['geometry']

		if boundary['type']=='LineString':
			polygons=[Path(np.array(boundary['coordinates']))]
		elif boundary['type']=='MultiLineString':
			polygons=[]
			for p in boundary['coordinates']:
				polygons.append(Path(np.array(p),closed=True))
		else:
			print('Error: unknown geometry')
			continue

		if region is not None and boundary is not None:
			counties[region]=polygons
	return counties,attributes

################################################################################
		
def create_netcdf(template,data,outname,template_var='pr'):
	# create outfile object
	outfile=Dataset(outname,'w')

	# Create dimensions copied from template file
	temp=template.variables[template_var]
	for dim in temp.dimensions:
		if dim[:3]=='lat' or dim[:3] =='lon':
			leng=len(template.dimensions[dim])
		
			outfile.createDimension(dim[:3],leng)
			outfile.createVariable(dim[:3],'f',(dim[:3],))
			outfile.variables[dim[:3]][:]=template.variables[dim][:]
			#print template.variables[dim].__dict__
			for att in template.variables[dim].ncattrs():
				outfile.variables[dim[:3]].__setattr__(att,template.variables[dim].__getattribute__(att))

	# Create data variable (named region_mask)
	outfile.createVariable('region_mask','f',['lat','lon'])
	outfile.variables['region_mask'][:]=(data-1)*-1
	
	#outfile.flush()
	outfile.close()

############################################################################

# Create a number of masks, from the shapefile, output text file
#
# Input Arguments:
# shapefile: path of shapefile containing polygons
# fieldname: attribute name in shapefile used to identify each field. 
# field_list: (optional)- specify a list (subset) of fields to create masks for, otherwise masks will be created covering all fields
#
def create_polygon_textfiles(shapefile,fieldname,field_list=None):

	# first create folders (if needed)
	if not os.path.exists('masks_text'):
		os.mkdir('masks_text')

	# Load Shape file
	regions=load_shapefile(shapefile,fieldname,field_list=field_list)

	print('Looping over regions and creating text files of masks')
	# Either loop over all regions, or list of regions specified by 'field_list'
	if field_list == None:
		field_list = regions.keys()
	for region in field_list:
		#region_ascii = unicodedata.normalize('NFKD',str(region).decode('utf-8')).encode('ascii','ignore')
		#region_ascii = str(unicodedata.normalize('NFKD',region).encode('ascii','ignore'))
		region_ascii = unicodedata.normalize('NFKD',region)

		print(fieldname,'=',region)
		polygons = regions[region]
		# Write polygons to text file
		with open('masks_text/mask_'+region_ascii+'.txt','w') as text_polygons:
			for p in polygons:
				add_to_text(text_polygons,p)

#############################################################################

# Create a textfile of polygons, combining multiple fields from the shapefile
#
# Input Arguments:
# shapefile: path of shapefile containing polygons
# fieldname: attribute name in shapefile used to identify each field. 
# field_list: (optional)- specify a list (subset) of fields to include in the mask, otherwise the mask combines all fields
# region_name: (optional)- name of region used in output filename
#
def create_combined_textfiles(shapefile, fieldname, field_list=None, region_name='region'):

	# first create folder (if needed)
	if not os.path.exists('masks_text'):
		os.mkdir('masks_text')

	# Load Shape file
	regions=load_shapefile(shapefile,fieldname,field_list=field_list)

	print('Looping over regions and combining masks')
	# Either loop over all regions, or list of regions specified by 'field_list'
	if field_list == None:
		field_list = regions.keys()
	
	with open('masks_text/mask_'+region_name+'.txt','w') as text_polygons:
		for region in field_list:
			print(fieldname,'=',region)
			polygons = regions[region]

			# Add polygon to text file
			for p in polygons:
				add_to_text(text_polygons,p)
			
###############################################################################

# Create a number of masks, from the shapefile, for a specific grid
# Area outside the polygons is True/1, area inside the polygons is False/0
#
# Input Arguments:
# f_grid: filename of netcdf file contatining grid information
# latname, lonname, template_var: variable names for latitude, longitude and a template variable in f_grid netcdf file
# shapefile: path of shapefile containing polygons
# fieldname: attribute name in shapefile used to identify each field. 
# field_list: (optional)- specify a list (subset) of fields to create masks for, otherwise masks will be created covering all fields
# plot, netcdf_out: (optional) booleans- whether or not to create output plot and/or output netcdf file
#
# Returns: dictionary of region_name:mask_array pairs
#
def create_masks(f_grid, shapefile, fieldname, field_list=None, latname='lat',lonname='lon',template_var='pr', plot=False, netcdf_out = False):

	# first create folders (if needed)
	if plot and not os.path.exists('plots'):
		os.mkdir('plots')
	if netcdf_out and not os.path.exists('masks_netcdf'):
		os.mkdir('masks_netcdf')

	# Load Shape file
	regions=load_shapefile(shapefile,fieldname,field_list=field_list)

	# Load lat lon grid (for mask)
	lonxx,latyy=load_grid(f_grid,latname=latname,lonname=lonname)
	nlat,nlon=lonxx.shape
	# Update lon to be from -180 to 180 
	# NOTE: (this is only if the shapefile uses lat coordinates from -180-180 )
	# Comment out otherwise
	lonxx[lonxx>180]=lonxx[lonxx>180]-360
	# Turn lat and lon into a list of coordinates
	points = np.vstack((lonxx.flatten(),latyy.flatten())).T

	if plot:
		from mpl_toolkits.basemap import Basemap,cm		
		# Set up Basemap projection (may need fine tuning)
		m = Basemap(projection = 'robin',lon_0=180)
		xx,yy=m(lonxx,latyy) # basemap coordinates

	# Either loop over all regions, or list of regions specified by 'field_list'
	if field_list == None:
		field_list = regions.keys()

	# Dictionary of masks
	masks={}

	# Do the loop
	print('Looping over regions and creating gridded masks')
	for region in field_list:
		region_ascii = unicodedata.normalize('NFKD',str(region).decode('utf-8')).encode('ascii','ignore')

		print(fieldname,'=',region)
		polygons = regions[region]
		# Create mask out of polygon, matching points from grid
		mask = create_mask(polygons,points,nlat,nlon)
		
		# Add to dictionary
		masks[region_ascii] = mask

		if netcdf_out:
			create_netcdf(Dataset(f_grid,'r'),mask,'masks_netcdf/mask_'+region_ascii+'.nc', template_var=template_var)

		if plot:
			plt.clf()
			m.contourf(xx,yy,mask)
			plt.colorbar()
			m.drawcoastlines(linewidth=0.2)
			m.drawcountries(linewidth=0.2)
			plt.title('Mask: '+region_ascii)
			plt.savefig('plots/mask_'+region_ascii+'.png')

	return masks

#############################################################################

# Create a mask, combining multiple fields from a shapefile, for a specific grid
# Area outside the polygons is True/1, area inside the polygons is False/0
#
# Input Arguments:
# f_grid: filename of netcdf file contatining grid information
# latname, lonname, template_var: variable names for latitude, longitude and a template variable in f_grid netcdf file
# shapefile: path of shapefile containing polygons
# fieldname: attribute name in shapefile used to identify each field. 
# field_list: (optional)- specify a list (subset) of fields to include in the mask, otherwise the mask combines all fields
# plot, netcdf_out: (optional) booleans- whether or not to create output plot and/or output netcdf file
# region_name: (optional)- name of region in output files
# 
# Returns array of the combined mask
#
def create_mask_combined(f_grid,shapefile,fieldname,field_list=None,region_name='region',latname='lat',lonname='lon',template_var='pr', plot=False,netcdf_out=False):

	# first create folders (if needed)
	if plot and not os.path.exists('plots'):
		os.mkdir('plots')
	if netcdf_out and not os.path.exists('masks_netcdf'):
		os.mkdir('masks_netcdf')

	# Load Shape file
	regions=load_shapefile(shapefile,fieldname,field_list=field_list)

	# Load lat lon grid (for mask)
	lonxx,latyy=load_grid(f_grid,latname=latname,lonname=lonname)
	nlat,nlon=lonxx.shape
	# Update lon to be from -180 to 180 
	# NOTE: (this is only if the shapefile uses lat coordinates from -180-180 )
	# Comment out otherwise
	lonxx[lonxx>180]=lonxx[lonxx>180]-360
	# Turn lat and lon into a list of coordinates
	points = np.vstack((lonxx.flatten(),latyy.flatten())).T


	combined_mask = np.zeros([nlat,nlon])
	if plot:
		from mpl_toolkits.basemap import Basemap,cm
		# Set up Basemap projection (may need fine tuning)
		m = Basemap(projection = 'robin',lon_0=180)
		xx,yy=m(lonxx,latyy) # basemap coordinates
		plot_mask = np.zeros([nlat,nlon])

	print('Looping over regions and combining masks')
	# Either loop over all regions, or list of regions specified by 'field_list'
	if field_list == None:
		field_list = regions.keys()

	i=1
	for region in field_list:
		print(fieldname,'=',region)
		polygons = regions[region]
		
		# Create mask out of polygon, matching points from grid
		mask = create_mask(polygons,points,nlat,nlon)
		print(mask.shape)
		combined_mask = combined_mask + (mask-1)*-1
		if plot:
			plot_mask = plot_mask + (mask-1)*-i
		i+=1
	
	# Create netcdf for combined mask
	if netcdf_out:
		create_netcdf(Dataset(f_grid,'r'),combined_mask,'masks_netcdf/mask_'+region_name+'.nc',template_var=template_var)

	if plot:
		plt.clf()
		m.contourf(xx,yy,plot_mask)
		plt.colorbar()
		m.drawcoastlines(linewidth=0.2)
		m.drawcountries(linewidth=0.2)
		plt.title('Mask: '+region_name)
		plt.savefig('plots/mask_'+region_name+'.png')

	return mask

#############################################################################

# Create a mask, from textfile for a specific grid
# Area outside the polygons is True/1, area inside the polygons is False/0
#
# Input Arguments:
# f_grid: (filename of netcdf file contatining grid information)
# textfile: path of text file containing coordinates of polygons
# latname, lonname, template_var: variable names for latitude, longitude and a template variable in f_grid netcdf file
# plot, netcdf_out: (optional) booleans- whether or not to create output plot and/or output netcdf file
# 
# Returns: mask array 
#
def create_mask_fromtext(f_grid, textfile, region_name='region', latname='lat', lonname='lon', template_var='pr', plot=False, netcdf_out=False):

	# first create folders (if needed)
	if plot and not os.path.exists('plots'):
		os.mkdir('plots')
	if netcdf_out and not os.path.exists('masks_netcdf'):
		os.mkdir('masks_netcdf')

	# Load Shape file
	polygons=load_polygons(textfile)

	# Load lat lon grid (for mask)
	lonxx,latyy=load_grid(f_grid,latname=latname,lonname=lonname)
	nlat,nlon=lonxx.shape
	# Update lon to be from -180 to 180 
	# NOTE: (this is only if the shapefile uses lat coordinates from -180-180 )
	# Comment out otherwise
	lonxx[lonxx>180]=lonxx[lonxx>180]-360
	# Turn lat and lon into a list of coordinates
	points = np.vstack((lonxx.flatten(),latyy.flatten())).T

	if plot:
		from mpl_toolkits.basemap import Basemap,cm
		# Set up Basemap projection (may need fine tuning)
		m = Basemap(projection = 'robin',lon_0=180)
		xx,yy=m(lonxx,latyy) # basemap coordinates
		plot_mask = np.zeros([nlat,nlon])

	# Create mask out of polygon, matching points from grid
	mask = create_mask(polygons,points,nlat,nlon)

	# Create netcdf for combined mask
	if netcdf_out:
		create_netcdf(Dataset(f_grid,'r'),mask,'masks_netcdf/mask_'+region_name+'.nc',template_var=template_var)

	if plot:
		plt.clf()
		m.contourf(xx,yy,mask)
		plt.colorbar()
		m.drawcoastlines(linewidth=0.2)
		m.drawcountries(linewidth=0.2)
		plt.title('Mask: '+region_name)
		plt.savefig('plots/mask_'+region_name+'.png')

	return mask

#################################################################################
#
# Examples are below
#
if __name__=='__main__':
	# Set up input files
	f_grid='/export/silurian/array-01/pu17449/processed_data_clim_2deg/NorESM1-HAPPI.pr.All-Hist_monclim_ensmean.nc'
	#f_grid = '/export/silurian/array-01/pu17449/processed_data_clim/CAM5-1-2-025degree.pr.All-Hist_monclim_ensmean.nc'
	
	f_text = '../../mask_sjoukje.txt'
	create_mask_fromtext(f_grid, f_text, region_name='test', latname='lat', lonname='lon', template_var='pr', plot=True, netcdf_out=False)

	# Shapefile with polygon vertices
	#continent = 'as'
	#lev = 3
	#shapefile='/export/silurian/array-01/pu17449/shapefiles/river_basins/hybas/'+continent+str(lev)+'/hybas_'+continent+'_lev'+str(lev).zfill(2)+'_v1c.shp'
	# Name of key to use to identify each region in the shapefile
	#fieldname = 'SORT'
	#area = [16,17,18]

#	create_mask_combined(f_grid,shapefile,fieldname,field_list=area,latname='lat',lonname='lon',template_var='pr', region_name='rivertest' plot=True, netcdf_out=True)

# Example: countries
#
#	shapefile = '/export/silurian/array-01/pu17449/shapefiles/countries/ne_10m_admin_0_countries.shp'
#	fieldname = 'NAME'
#	# Create mask for each country in shapefile
#	create_masks(f_grid,shapefile,fieldname,latname='lat',lonname='lon',template_var='pr',plot=True)

# Example for Rob:
#
# Create mask for bangladesh
#	mask_array = create_mask_combined(f_grid,shapefile,fieldname,field_list=['Bangladesh'],region_name='Bangladesh',latname='lat',lonname='lon',plot=True,netcdf_out=True, template_var='pr')
	

# Example for Fredi:
#	
	# Template grid for output mask
#	f_grid='/export/silurian/array-01/pu17449/ie2sga.pdl4dec.nc'

	# Text file with region boundaries
#	f_text = '/home/bridge/pu17449/src/happi_analysis/river_basins/masks_text/mask_ganges.txt'

#	create_mask_fromtext(f_grid,f_text, region_name='ganges', latname='global_latitude0',lonname='global_longitude0', template_var='field8', plot=True, netcdf_out=False)

# Example for Sjoukje (creating text files)
#
#	shapefile = '/export/anthropocene/array-01/pu17449/basins_sjoukje/Aqueduct_river_basins_GANGES - BRAHMAPUTRA.shp'
#	fieldname = 'BASIN_NAME'
#	create_polygon_textfiles(shapefile,fieldname)

