# *********************************************************************************************
# *********************************************************************************************
# PRJ parser1.py
# author - Troy Hicks, troy.hicks@alaska.gov   
#
# **** Most recent revision done on Feb 15,2019
#     

#
# The PURPOSE of this python script is to parse a PRJ file.
#

#
# important variables
# -----------------------
# a         Semi-major axis, (Equatorial Radius) in meters [6378137.0 for GRS80]
# inFl      Inverse flattening (1/f) unitless [298.257222100882 for GRS80]
# FE        False Easting, in meters
# FN        False Northing, in meters
# phiO      latitude of true projection origin [for our purposes a single parallel definition will always be used, also as origin]
# phi1,phi2 phi1 and phi2 are the standard parallels for a double, for a single they must share the same value
# lamO      Central Meridian, longitude at grid origin [for our purposes will always be used for both]
# kO        grid scale factor at Central Parallel
# skew      skew azimuth for hotine oblique mercator, lamO used for center, phiO used for center
# *********************************************************************************************

import ast
import os.path


# **** READ A PRJ FILE ****

# ************ hard code a specific file to open **************
#with open("PRJ\Cordova.prj") as file:  
#    prjdata = file.read()
    
# ************ hard code in a .prj definition to test out ******
#prjFairbanks = 'PROJCS["Fairbanks",GEOGCS["GCS_NAD_1983_2011",DATUM["D_NAD_1983_2011",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["False_Easting",800000],PARAMETER["False_Northing",200000],PARAMETER["Central_Meridian",-146.933333333333],PARAMETER["Standard_Parallel_1",64.85],PARAMETER["Standard_Parallel_2",64.85],PARAMETER["Scale_Factor",1.00003],PARAMETER["Latitude_Of_Origin",64.85],UNIT["Foot_US",0.304800609601219]]'
#prjdata = prjFairbanks

def prj_reader(filename):
# check if filename exists
  if not os.path.isfile(filename):
    print ("I can't find .prj file " + filename + " to read.")
    return 0
# read and parse .prj file
  handle = open(filename, "r")
  prjdata = handle.read().replace('\n', '')
  handle.close()

# check if prj is a PROJCS or might just be a GEOGCS
  if prjdata.find("PROJCS") != 0:
    print('The ' + filename + " file does not seem to be a Projected Coordinate System 'PROJCS':")
    if prjdata.find("GEOGCS") == 0:
      print('The ' + filename + " file appears to be a Geographic Coordinate System 'GEOGCS':")
    return 0

# Since it is in fact a PROJC, then lets make it
# become a literal expression in Python, a list of lists 
  temp = prjdata  
  temp = temp.replace("PROJCS", "")
  temp = temp.replace("GEOGCS", "")
  temp = temp.replace("DATUM", "")
  temp = temp.replace("SPHEROID", "")
  temp = temp.replace("PRIMEM", "")
  temp = temp.replace("UNIT", "")
  temp = temp.replace("PROJECTION", "")
  temp = temp.replace("PARAMETER", "")

# myCRS is a literal list of lists of the prj string. CRS = cordinate reference system
  myCRS = ast.literal_eval(temp)  



# check for allowable projections
  if myCRS[2][0] != 'Lambert_Conformal_Conic' and myCRS[2][0] != 'Hotine_Oblique_Mercator_Azimuth_Natural_Origin' and myCRS[2][0] != 'Transverse_Mercator':
    print('The ' + filename + ' file is not in an allowable projection:')
    print('  ' + str(myCRS[2]))
    print('Hey get with the program! Only CONFORMAL projections allowed. And only these specific types:')
    print('  Transverse_Mercator')
    print('  Hotine_Oblique_Mercator_Azimuth_Natural_Origin')
    print('  Lambert_Conformal_Conic')
    return 0

# return projection information
  PROJCS_Name    = myCRS[0]        # PROJCS = ["PROJCS_name", GEOGCS, PROJECTION, PARAMETER, ... , PARAMETER, UNIT]
  GEOGCS         = myCRS[1]        # GEOGCS = ["GEOGCS_type", DATUM, PRIMEM, UNIT]
  DATUM          = GEOGCS[1]       # DATUM = ["DATUM_type", SPHEROID]
  SPHEROID       = DATUM[1]        # SPHEROID = ["SPHEROID_name", 0.0, 0.0]
  PRIMEM         = GEOGCS[2]       # PRIMEM = ["PRIMEM_name", 0.0]
  a              = SPHEROID[1]     # semi-major axis, in meters
  inFl           = SPHEROID[2]     # inverse flattening
  PROJECTION     = myCRS[2]        # PROJECTION = ["PROJECTION_method"]

# if Lambert_Conformal_Conic
  if myCRS[2][0] == 'Lambert_Conformal_Conic':
    for i in range(7) :
      if myCRS[i + 3][0] == "False_Easting":
        FE             = myCRS[i + 3][1]    
      if myCRS[i + 3][0] == "False_Northing":
        FN             = myCRS[i + 3][1]      
      if myCRS[i + 3][0] == "Latitude_Of_Origin":
        phiO           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Standard_Parallel_1":
        phi1           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Standard_Parallel_2":
        phi2           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Central_Meridian":
        lamO           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Scale_Factor":
        kO             = myCRS[i + 3][1]
    UNIT           = myCRS[10]       # UNIT = ["UNIT_type", 0.0]
    return PROJCS_Name, GEOGCS, DATUM, SPHEROID, PRIMEM, a, inFl, PROJECTION, FE, FN, phiO, phi1, phi2, lamO, kO, UNIT

# if 'Transverse_Mercator'
  if myCRS[2][0] == 'Transverse_Mercator':
    for i in range(5) :
      if myCRS[i + 3][0] == "False_Easting":
        FE             = myCRS[i + 3][1]    
      if myCRS[i + 3][0] == "False_Northing":
        FN             = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Latitude_Of_Origin":
        phiO           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Central_Meridian":
        lamO           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Scale_Factor":
        kO             = myCRS[i + 3][1]
    UNIT           = myCRS[8]       # UNIT = ["UNIT_type", 0.0]
    return PROJCS_Name, GEOGCS, DATUM, SPHEROID, PRIMEM, a, inFl, PROJECTION, FE, FN, phiO, lamO, kO, UNIT

# if 'Hotine_Oblique_Mercator_Azimuth_Natural_Origin'
  if myCRS[2][0] == 'Hotine_Oblique_Mercator_Azimuth_Natural_Origin':
    for i in range(6) :
      if myCRS[i + 3][0] == "False_Easting":
        FE             = myCRS[i + 3][1]    
      if myCRS[i + 3][0] == "False_Northing":
        FN             = myCRS[i + 3][1]      
      if myCRS[i + 3][0] == "Scale_Factor":
        kO             = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Azimuth":
        skew           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Longitude_Of_Center":
        lamO           = myCRS[i + 3][1]
      if myCRS[i + 3][0] == "Latitude_Of_Center":
        phiO           = myCRS[i + 3][1]
    UNIT           = myCRS[9]       # UNIT = ["UNIT_type", 0.0]
    return PROJCS_Name, GEOGCS, DATUM, SPHEROID, PRIMEM, a, inFl, PROJECTION, FE, FN, phiO, lamO, kO, skew, UNIT 

 