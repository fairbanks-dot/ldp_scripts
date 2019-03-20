# *********************************************************************************************
# *********************************************************************************************
# TransverseMercator.py
# current author    - Troy Hicks,  Dec 2017 - Mar 2018  troy.hicks@alaska.gov
#
# last edit: Mar. 16, 2019 
#
#
# still needs tested to make sure it works.
#
# The purpose of this code is to perform Transverse Mercator direct computations such that when provided
# a geodetic position in decimal degrees to return a x,y (grid values), k (grid scale factor at a point),
# and mercon (meridian convergence angle). Angle is in radians. To go from grid values backward to geodetic requires an inverse
# computation. Not included here.
#
# To use this:
# import TransverseMercator as TM
# get values for the variables needed by TM, phi, phiO, lam, lamO, FN, FE, kO, a, inFl.
# an hard code them in or open a prj file and read them in.
# results = TM.compute_TM(phi, phiO, lam, lamO, FN, FE, kO, a, inFl)
# this will create a list [x,y,k,merCon]
#
# *********************************************************************************************

import numpy as np # this is importing numerical tools. numpy needs installed for this to work.

   
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# Transverse Mercator (TM) supporting functions and definitions
#
# variables and constants, from Coordinate Reference System definition
# -------------------------------------------------------------------
# phi         geodetic latitude, positive North, in decimal degree
# lam         geodetic longitude, positive East, in decimal degree
# phiO        Latitude of grid origin
# lamO        Central Meridian, and meridian of grid origin
# FN, FE      false northing and false easting, in meters
# kO          grid scale factor at Central Meridian, lamO
# k           grid scale factor at a general point
# a           semi-major axis, from ellipsoid definition
# inFl        = 1/f, where f is from ellipsoid definition

# -- computed variables --
# b           semi-minor axis, from ellipsoid definition
# fiEcSq      ellipsoid first eccentricity squared (e**2)
# seEcSq      ellipsoid second eccentricity squared (e^2 prime), or (e')^2 
# N,T,C,A,M   intermediate variables and equations from Snyder
# C1,C3,C5,t  intermediate variables from Stem page 34, for convergence angle computation
# x           computed Easting (grid), in meters
# y           computed Northing (grid), in meters
# mercon      convergence angle

# math is adapted from page 61,       "Map Projections - A working Manual", Snyder, 1997 4th printing
#                                     "State Plane Coordinate System of 1983" Stem, James, mar 1990. NOAA Manual NOS NGS 5
#                                     "Ground Truth: Optimized Design of Low Distortion Projections", Dennis, Michael, Mar 2017
#                                     "IOGP Publication 373-7-2 –  Geomatics Guidance Note number 7, part 2" October 2018 www.epsg.org


# ----------------------------------------------------------------------------------------------
# main function
def compute_TM(phi, phiO, lam, lamO, FN, FE, kO, a, inFl): # generic, requires phi, lam, lamO in decimal degree
   phi = np.deg2rad(phi)
   phiO = np.deg2rad(phiO)
   lam = np.deg2rad(lam)
   lamO = np.deg2rad(lamO)
   fiEcSq = 2 / inFl - inFl ** -2 # ellipsoid first eccentricity squared
   seEcSq = fiEcSq/ (1 - fiEcSq) # ellipsoid second eccentricity squared
   N = a / np.sqrt(1 - fiEcSq * np.sin(phi) ** 2) # N, page 61 Snyder
   T = np.tan(phi) ** 2 # T, page 61 Snyder, also same 't'**0.5 page 33 Stem  note that Stem is tan(phi), not squared as is here
   C = seEcSq * np.cos(phi) ** 2 # C, page 61 Snyder, same as (eta)**2 on page 33 Stem
   A = (lam - lamO) * np.cos(phi) # A, page 61 Snyder - *** I changed the order, to change sign.  Also L page 34 Stem
   # next line is M and Mo on page 61 Snyder, adjusted since fiExSq is square of e, and not actually e.
   M = a  * ((1 - fiEcSq / 4 - 3 * fiEcSq ** 2 / 64 - 5 * fiEcSq ** 3 / 256) * phi - (3 * fiEcSq / 8 + 3 * fiEcSq ** 2 / 32 + 45 * fiEcSq ** 3 / 1024) * np.sin(2 * phi) + (15 * fiEcSq ** 2 / 256 + 45 * fiEcSq ** 3 / 1024) * np.sin(4 * phi) - (35 * fiEcSq ** 3 / 3072) * np.sin(6 * phi))
   # again, but use phiO instead of phi
   MO = a  * ((1 - fiEcSq / 4 - 3 * fiEcSq ** 2 / 64 - 5 * fiEcSq ** 3 / 256) * phi - (3 * fiEcSq / 8 + 3 * fiEcSq ** 2 / 32 + 45 * fiEcSq ** 3 / 1024) * np.sin(2 * phiO) + (15 * fiEcSq ** 2 / 256 + 45 * fiEcSq ** 3 / 1024) * np.sin(4 * phiO) - (35 * fiEcSq ** 3 / 3072) * np.sin(6 * phiO))
   x = FE + kO * N * (A + (1 - T + C) * A ** 3 / 6 + (5 - 18 * T + T ** 2 + 72 * C - 58 * seEcSq) * A ** 5 / 120)
   y = FN + kO * (M - MO + N * np.tan(phi) * (A ** 2 / 2 + (5 - T + 9 * C + 4 * C ** 2) * A ** 4 / 24 + (61 - 58 * T + T ** 2 + 600 * C - 330 * seEcSq) * A ** 6 / 720))
   k = kO * (1 + (1 + C) * A ** 2 / 2 + (5 - 4 * T + 42 * C + 13* C ** 2 - 28 * seEcSq) * A ** 4 / 24 + (61 - 148* T + 16 * T **2) * A ** 6 /720)
   # next five lines are adapted from page 34 Stem, to get convergence. Not shown in Snyder.
   t  = np.sqrt(T)
   C1 = -1 * t
   C3 = 1/3 * (1 + 3*C + 2*C**2)
   C5 = 1/15 * (2 - t**2)
   merCon = C1 * A * (1 + A**2 * (C3 + C5 * A**2))  # here 'A' is same as 'L' in page 34 Stem
   return x,y,k,merCon


