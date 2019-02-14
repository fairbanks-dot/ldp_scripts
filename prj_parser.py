#
# PRJ Parser - prj_parser.py
#
# author - Troy Hicks,  Nov 2018  troy.hicks@alaska.gov   
#        - Jason Grimes, Feb 2019 jason@gina.alaska.edu
#

import os.path
import osr

################################################
# Read a .prj projection file and return a dictionary of the contents.
#
# {filename} is the path and filename of the prj file to load
################################################

def prj_reader(filename):
  # check if filename exists
  if not os.path.isfile(filename):
  	print ("I can't find .prj file " + filename + " to read.")
  	return 0

  # setup dictionary
  prj_dict = {}

  # read and parse .prj file
  handle = open(filename, "r")
  prjdata = handle.read().replace('\n', '')
  handle.close()

  srs = osr.SpatialReference()
  srs.ImportFromWkt(prjdata)

  # build projection dictionary
  #
  # check for allowable projections
  prj_dict['PROJECTION'] = srs.GetAttrValue('projection')
  if prj_dict['PROJECTION'] != 'Lambert_Conformal_Conic' and prj_dict['PROJECTION'] != 'Hotine_Oblique_Mercator_Azimuth_Natural_Origin' and prj_dict['PROJECTION'] != 'Transverse_Mercator':

     print('The ' + filename + ' file is not in an allowable projection:')
     print('  ' + prj_dict['PROJECTION'])
     print('The allowable projections are as follows...')
     print('  Transverse_Mercator')
     print('  Hotine_Oblique_Mercator_Azimuth_Natural_Origin')
     print('  Lambert_Conformal_Conic')
     return 0

  # the projection is good, load the rest of the values
  prj_dict['PROJCS'] = srs.GetAttrValue('projcs')
  prj_dict['GEOGCS'] = srs.GetAttrValue('geogcs')
  prj_dict['DATUM'] = srs.GetAttrValue('datum')
  prj_dict['SPHEROID'] = [srs.GetAttrValue('spheroid', 0), srs.GetAttrValue('spheroid', 1), srs.GetAttrValue('spheroid', 2)]
  prj_dict['PRIMEM'] = [srs.GetAttrValue('primem', 0), srs.GetAttrValue('primem', 1)]
  prj_dict['PROJECTION'] = srs.GetAttrValue('projection')
  prj_dict['PARAMETER_FE'] = srs.GetProjParm(osr.SRS_PP_FALSE_EASTING)
  prj_dict['PARAMETER_FN'] = srs.GetProjParm(osr.SRS_PP_FALSE_NORTHING)
  prj_dict['PARAMETER_CM'] = srs.GetProjParm(osr.SRS_PP_CENTRAL_MERIDIAN)
  prj_dict['PARAMETER_SP1'] = srs.GetProjParm(osr.SRS_PP_STANDARD_PARALLEL_1)
  prj_dict['PARAMETER_SP2'] = srs.GetProjParm(osr.SRS_PP_STANDARD_PARALLEL_2)
  prj_dict['PARAMETER_SF'] = srs.GetProjParm(osr.SRS_PP_SCALE_FACTOR)
  prj_dict['PARAMETER_LO'] = srs.GetProjParm(osr.SRS_PP_LATITUDE_OF_ORIGIN)
  prj_dict['UNIT'] = [srs.GetAttrValue('unit', 0), srs.GetAttrValue('unit', 1)]

  print(filename + ' has been loaded')
  print(prj_dict)

  # return the populated dictionary
  return prj_dict

################################################
# Write a .prj projection file from a dictionary
#
# {filename} is the path and filename of the prj file to write
# {prj_dict} is the dictionary of the projection attributes and parameters
################################################

def prj_writer(filename, prj_dict):
  srs = osr.SpatialReference()

  # populate the srs with the values from the dictionary
  srs.SetGeogCS(prj_dict['GEOGCS'], prj_dict['DATUM'], prj_dict['SPHEROID'][0], float(prj_dict['SPHEROID'][1]), float(prj_dict['SPHEROID'][2]), prj_dict['PRIMEM'][0], float(prj_dict['PRIMEM'][1]), prj_dict['UNIT'][0], float(prj_dict['UNIT'][1]))
  srs.SetProjCS(prj_dict['PROJCS'])
  srs.SetProjection(prj_dict['PROJECTION'])
  srs.SetProjParm('False_Easting', float(prj_dict['PARAMETER_FE']))
  srs.SetProjParm('False_Northing', float(prj_dict['PARAMETER_FN']))
  srs.SetProjParm('Central_Meridian', float(prj_dict['PARAMETER_CM']))
  srs.SetProjParm('Standard_Parallel_1', float(prj_dict['PARAMETER_SP1']))
  srs.SetProjParm('Standard_Parallel_2', float(prj_dict['PARAMETER_SP2']))
  srs.SetProjParm('Scale_Factor', float(prj_dict['PARAMETER_SF']))
  srs.SetProjParm('Latitude_Of_Origin', float(prj_dict['PARAMETER_LO']))

  # write prj file
  handle = open(filename, "w")
  handle.write(srs.ExportToWkt())
  handle.close

  print('prj data written to ' + filename)
  return 1
