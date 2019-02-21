#!/usr/bin/env python3 
# *********************************************************************************************
# *********************************************************************************************
# Parser_tester.py
# author - Troy Hicks, troy.hicks@alaska.gov   
#
# **** Most recent revision done on Feb 18,2019
#     

#
# The PURPOSE of this python script is to test PRJ parser1.
# PRJ parser1.py is supposed to parse PRJ files for the explicit purpose of
# reading them for LDP design. It is also supposed to do this without any other libraries
# being required. A version that uses GDAL is powerful but requires it added as a package.
# this is a possible way to keep it simpler.
#


# import local modules
import PRJ_parser1 as p
from tkinter import filedialog
from tkinter import *
 
# look up the prj file 
root = Tk()
root.filename =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("prj files","*.prj"),("all files","*.*")))
prj_data = p.prj_reader(root.filename)
if prj_data == 0:
    print('Error reading the prj file, see error above!')
    exit(1)


# Here is the writer, which should be moved to a method in the parser or something.
# print out values for proof it works. This is just a quick look for testing.
if prj_data[7][0] == 'Lambert_Conformal_Conic':
  prj_string = '[PROJCS["' + prj_data[0] + '"'
  prj_string = prj_string + ',GEOCS["' + prj_data[1][0] + '",'
  prj_string = prj_string + 'DATUM["' +prj_data[1][1][0] + '",'
  prj_string = prj_string + 'SPHEROID["' + prj_data[1][1][1][0] + '",' + str(prj_data[1][1][1][1]) +"," + str(prj_data[1][1][1][2]) + "]],"
  prj_string = prj_string + 'PRIMEM["' +  prj_data[1][2][0] + '",' + str(prj_data[1][2][1]) + '],'
  prj_string = prj_string + 'UNIT["' + prj_data[1][3][0] + '",' + str(prj_data[1][3][1]) + "]],"
  prj_string = prj_string + 'PROJECTION["' +prj_data[7][0] + '"],'
  prj_string = prj_string + 'PARAMETER["False_Easting",' + str(prj_data[8]) + '],'
  prj_string = prj_string + 'PARAMETER["False_Northing",' + str(prj_data[9]) + '],'
  prj_string = prj_string + 'PARAMETER["Central_Meridian",' + str(prj_data[10]) + '],'
  prj_string = prj_string + 'PARAMETER["Standard_Parallel_1",' + str(prj_data[11]) + '],'
  prj_string = prj_string + 'PARAMETER["Standard_Parallel_2",' + str(prj_data[12]) + '],'
  prj_string = prj_string + 'PARAMETER["Scale_Factor",' + str(prj_data[14]) + '],'
  prj_string = prj_string + 'PARAMETER["Latitude_Of_Origin",' + str(prj_data[13]) + '],'
  prj_string = prj_string + "UNIT[" + prj_data[15][0] + "," + str(prj_data[15][1]) + ']]'
  print(prj_string)
  
if prj_data[7][0] == 'Transverse_Mercator':
  prj_string = '[PROJCS["' + prj_data[0] + '"'
  prj_string = prj_string + ',GEOCS["' + prj_data[1][0] + '",'
  prj_string = prj_string + 'DATUM["' +prj_data[1][1][0] + '",'
  prj_string = prj_string + 'SPHEROID["' + prj_data[1][1][1][0] + '",' + str(prj_data[1][1][1][1]) +"," + str(prj_data[1][1][1][2]) + "]],"
  prj_string = prj_string + 'PRIMEM["' +  prj_data[1][2][0] + '",' + str(prj_data[1][2][1]) + '],'
  prj_string = prj_string + 'UNIT["' + prj_data[1][3][0] + '",' + str(prj_data[1][3][1]) + "]],"
  prj_string = prj_string + 'PROJECTION["' +prj_data[7][0] + '"],'
  prj_string = prj_string + 'PARAMETER["False_Easting",' + str(prj_data[8]) + '],'
  prj_string = prj_string + 'PARAMETER["False_Northing",' + str(prj_data[9]) + '],'
  prj_string = prj_string + 'PARAMETER["Central_Meridian",' + str(prj_data[10]) + '],'
  prj_string = prj_string + 'PARAMETER["Scale_Factor",' + str(prj_data[12]) + '],'
  prj_string = prj_string + 'PARAMETER["Latitude_Of_Origin",' + str(prj_data[11]) + '],'
  prj_string = prj_string + "UNIT[" + prj_data[13][0] + "," + str(prj_data[13][1]) + ']]'
  print(prj_string)

if prj_data[7][0] == 'Hotine_Oblique_Mercator_Azimuth_Natural_Origin':
  prj_string = '[PROJCS["' + prj_data[0] + '"'
  prj_string = prj_string + ',GEOCS["' + prj_data[1][0] + '",'
  prj_string = prj_string + 'DATUM["' +prj_data[1][1][0] + '",'
  prj_string = prj_string + 'SPHEROID["' + prj_data[1][1][1][0] + '",' + str(prj_data[1][1][1][1]) +"," + str(prj_data[1][1][1][2]) + "]],"
  prj_string = prj_string + 'PRIMEM["' +  prj_data[1][2][0] + '",' + str(prj_data[1][2][1]) + '],'
  prj_string = prj_string + 'UNIT["' + prj_data[1][3][0] + '",' + str(prj_data[1][3][1]) + "]],"
  prj_string = prj_string + 'PROJECTION["' +prj_data[7][0] + '"],'
  prj_string = prj_string + 'PARAMETER["False_Easting",' + str(prj_data[8]) + '],'
  prj_string = prj_string + 'PARAMETER["False_Northing",' + str(prj_data[9]) + '],'
  prj_string = prj_string + 'PARAMETER["Scale_Factor",' + str(prj_data[12]) + '],'
  prj_string = prj_string + 'PARAMETER["Azimuth",' + str(prj_data[13]) + '],'
  prj_string = prj_string + 'PARAMETER["Longitude_Of_Center",' + str(prj_data[11]) + '],'
  prj_string = prj_string + 'PARAMETER["Latitude_Of_Center",' + str(prj_data[10]) + '],'
  prj_string = prj_string + "UNIT[" + prj_data[14][0] + "," + str(prj_data[14][1]) + ']]'
  print(prj_string)
