# *********************************************************************************************
# *********************************************************************************************
# PRJ tester.py
# author - Troy Hicks,  Nov 2018  troy.hicks@alaska.gov   
#
# **** Most recent revision done on Feb 02,2019
#
# **** NOTES for future edits
#        I realize that I should use tuples of tuples instead of lists, or at least save the first list of list made as a tuple
#        because the orders of the list elements, is not mutable. It cannot be sorted, would ruin the structure, therefore really it is
#        a tuple.  lists use []. Tuples use (). Maybe read as string, trim, then save as list, then save right away as tuple. Maybe.
# 
# **** main edits:
#      Feb. 02, 2019 - minor edits to the documentation to clarify and remove typos, prepare to share with Pete Hickman as a draft
#      Dec. 01, 2018 - now read in a .prj file for testing. Testing on Cordova.prj
#      Nov.     2018 - successfully parsed a string of .prj contents, Fairbanks.prj, which is LCC, and viewed results.
# ****
#
# The PURPOSE of this python script is to help figure out how to code in a PRJ reader, that will parse. And eventually code a PRJ writer.
#
# NOTES about .prj
# .prj is from Esri, and uses well known text (WKT) to make a world file to define a (CRS) cordinate reference system.
# It is not exactly following WKT, but is nicely simplified and assumes unit types and so on. This is NOT truly a WKT 2 format.
# So why use it? Because it works and its easier to find software to read and write .prj now than pure WKT.
#
# Notes about .prj structure
# The CRS will either, only, be a geographic 'geodetic' coordinate system (GEOGCS) or a projected 'grid' coordinate system (PROJCS).
# If it is PROJCS then it requires both a GEOGCS and also a projection definition.
#
# GEOGCS is defined by DATUM and Prime Meridian and UNIT (usually degrees).
#  DATUM is defined by SPHEROID.
#   SPHEROID is defined by 'a' (equatorial radius of the reference ellipsoid) and 'f' (flattening of the ellipsoid)
#  PRIMEM
#  UNIT
#
# PROJCS is defined by a GEOGCS and PROJECTION (type) and PARAMETERS (# of being dependendent on projection type) and UNIT (usually meters)
#  GEOGCS
#  PROJECTION this is the projection method, name of the method
#  PARAMETER usually has a few of these, not one. 
#  PARAMETER
#  ...
#  UNIT

# I realized that the .prj file looked like a Python list of lists, but with keyword headers. Kind of like
# a Python dictionary of lists. I don't think that is a possible data type. So for now will reduce it to a list of lists.
#
# The .prj seems it could be reduced to these sub-components lists:
#PRIMEM = ["PRIMEM_name", 0.0]           list = [string, number]
#SPHEROID = ["SPHEROID_name", 0.0, 0.0]  list = [string, number, number]
#UNIT = ["UNIT_type", 0.0]               list = [string, number]
#DATUM = ["DATUM_type", SPHEROID]        list = [string, list] => [string, [string, number, number]]
#PROJECTION = ["PROJECTION_type"]        list = [string] the string is exactly the name of the projection type, such as "Lambert_Conformal_Conic"
#PARAMETER = ["PARAMETER_method", 0.0]   list = [string, number] the string is exactly a named parameter that some PROJECTION needs, such as "Central_Meridian"
#GEOGCS = ["GEOGCS_type", DATUM, PRIMEM, UNIT] list = [string, list, list, list]
#PROJCS = ["PROJCS_name", GEOGCS, PROJECTION, PARAMETER, ... , PARAMETER, UNIT] list = [string, list, list, list, ..., list] 
#
# Notice that PROJECTION is not a list of PARAMETER's, but rather PROJECTION and each PARAMETER and UNIT and GEOGCS are the parts
# of the PROJCS. A person might think of a projection definition as being its own thing, and therefore might ought to have been
# a list itself, such as [string, list] but that is NOT how this is done in .prj. The PROJCS is the projection definition as well
# as the specific GEOGCS combined to make a CRS. Really any projection definition is independent of the GEOGCS. So ... a person
# was probably correct to have expected the format to be different. A CRS that is grid (PROJCS) should be made up of three things
# such as [name, GEOGCS, PROJ, UNITS] where PROJ is the projection definition. PROJ might then be a list = [type, parameter, ..., parameter]
# Maybe not though, as where do you put the False Eastings, etc. As part of the PROJ, or in the PROJCS. Until you make a PROJCS some
# of the parameters dont matter really. So I can see why it is just all embedded directly in PROJCS. Because really PROJECTION is a method.
# Once you decide the false easting, false northing, origin, and the GEOGCS you then have a full CRS.
#
#
# *** # parameters needed by PROJECTION type. LCC = 7, TM = 5, OM = 6 ***
#PROJCS_LCC = ["PROJCS_LCC_name",GEOGCS,PROJECTION,PARAMETER,PARAMETER,PARAMETER,PARAMETER,PARAMETER,PARAMETER,PARAMETER,UNIT]
#PROJCS_TM = ["PROJCS_TM_name",GEOGCS,PROJECTION,PARAMETER,PARAMETER,PARAMETER,PARAMETER,PARAMETER,UNIT]
#PROJCS_OM = ["PROJCS_OM_name",GEOGCS,PROJECTION,PARAMETER,PARAMETER,PARAMETER,PARAMETER,PARAMETER,UNIT]
#
# Work Flow:
# Read in the contents of a typical .prj file, and save contents to a variable
# Then parse the results in order to get each main part of the .prj as sub-components.
#    In order to parse:
#       Seems that the list of lists, should always have a few items in the same index position.
#       list[0] should always be the PROJCS name.
#       list[1] should always be the GEOGCS list
#       list[2] should always be the PROJECTION type
#       list[3] ... list[should be PARAMETERS, how even many there are. Each is a list.
# Assign the sub-components to variables so they are available for use.
# *********************************************************************************************

# Bring in the Abstract Syntax Tree library for Python. This will let a literal evaluation of an expression occur as if it is python code.
# it will be used to read a .prj file contents in and strip a few parts so it is a python list of lists.
import ast





# ************ here try a result of reading a PRJ file **************
# assume a typical prj file is read into this variable as a string
#with open("PRJ\Cordova.prj") as file:  
#    data = file.read()
    
# ************ or here we just code in a .prj defintion to test out ******
prjFairbanks = 'PROJCS["Fairbanks",GEOGCS["GCS_NAD_1983_2011",DATUM["D_NAD_1983_2011",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["False_Easting",800000],PARAMETER["False_Northing",200000],PARAMETER["Central_Meridian",-146.933333333333],PARAMETER["Standard_Parallel_1",64.85],PARAMETER["Standard_Parallel_2",64.85],PARAMETER["Scale_Factor",1.00003],PARAMETER["Latitude_Of_Origin",64.85],UNIT["Foot_US",0.304800609601219]]'
data = prjFairbanks

# Now remove text to make it
# become a literal expression in Python, a list of lists! [ [],[[],[]],[] ]  --- for example
temp = data  
temp = temp.replace("PROJCS", "")
temp = temp.replace("GEOGCS", "")
temp = temp.replace("DATUM", "")
temp = temp.replace("SPHEROID", "")
temp = temp.replace("PRIMEM", "")
temp = temp.replace("UNIT", "")
temp = temp.replace("PROJECTION", "")
temp = temp.replace("PARAMETER", "")

# myCRS is a literal list of lists of the prj string. CRS = cordinate reference system
myCRS = ast.literal_eval(temp)  # basically this myCRS is the PROJCS, maybe I should just use this as variable name


# CAUTION: Maybe I should be using tuples. Tuples are immutable, cannot have order changed or values. Order matters!
# Here I dump the sub lists into list variables  to use later. This only works if order of PRJ is same.
# Probably should read myCRS[2] and find out what type of projection, and therefore which parameters to look for.
# Then can use if then and read thru the list elements for

# need to look at myCRS[2] first and see which type of projection it is, then based on the answer use one of
# three ways to parse. Because a Lambert, Transverse mercator, and a Oblique Mercator all exactly have a set number
# of PARAMETERS. below I just have it for the LAMBERT for quick testing.

# if myCRS[2] is Lambert then do this, else if other projection type do something else, so on.

PROJCS_Name    = myCRS[0]        # PROJCS = ["PROJCS_name", GEOGCS, PROJECTION, PARAMETER, ... , PARAMETER, UNIT]
GEOGCS         = myCRS[1]        # GEOGCS = ["GEOGCS_type", DATUM, PRIMEM, UNIT]
DATUM          = GEOGCS[1]       # DATUM = ["DATUM_type", SPHEROID]
SPHEROID       = DATUM[1]        # SPHEROID = ["SPHEROID_name", 0.0, 0.0]
PRIMEM         = GEOGCS[2]       # PRIMEM = ["PRIMEM_name", 0.0]
PROJECTION     = myCRS[2]        # PROJECTION = ["PROJECTION_method"]
PARAMETER_FE   = myCRS[3]        # PARAMETER = ["PARAMETER_type", 0.0]
PARAMETER_FN   = myCRS[4]        # PARAMETER = ["PARAMETER_type", 0.0]
PARAMETER_CM   = myCRS[5]        # PARAMETER = ["PARAMETER_type", 0.0]
PARAMETER_SP_1 = myCRS[6]        # PARAMETER = ["PARAMETER_type", 0.0]
PARAMETER_SP_2 = myCRS[7]        # PARAMETER = ["PARAMETER_type", 0.0]
PARAMETER_SF   = myCRS[8]        # PARAMETER = ["PARAMETER_type", 0.0]
PARAMETER_LO   = myCRS[9]        # PARAMETER = ["PARAMETER_type", 0.0]
UNIT           = myCRS[10]       # UNIT = ["UNIT_type", 0.0]
# only grabbed one of the UNITS, because the one inside the DATUM has not yet been needed.

# print out values for proof it works. would dump into variables to pass to other scripts.
# print(myCRS)
print(myCRS[0])
print(myCRS[1])
print(myCRS[2])
print(myCRS[3])
print(myCRS[4])
print(myCRS[5])
print(myCRS[6])
print(myCRS[7])
print(myCRS[8])
print(myCRS[9])
print(myCRS[10])


