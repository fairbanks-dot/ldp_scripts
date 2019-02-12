# *********************************************************************************************
# *********************************************************************************************
# LambertConformalConic.py
# current author    - Troy Hicks,  Dec 2017 - Mar 2018  troy.hicks@alaska.gov
#
# last edit: Feb. 04, 2019
#
# The purpose of this code is to perform Lambert Conformal Conic direct computations such that when provided
# a geodetic position in decimal degrees to return a x,y (grid values), k (grid scale factor at a point),
# and mercon (convergence angle). Angle is in radians. To go from grid values backward to geodetic requires an inverse
# computation. Not included here.
#
# This is a revision of past code, to allow for double parallel Lambert definitions
# where first efforts assumed single. But I notice that .prj often refers to single as a double
# and just has both be the same value. I wanted to be able to compute on a real double so that
# I could see what happens when someone misuses a single and accidentily enters it wrong in their
# survey software.
#
# To use this:
# import LambertConformalConic as LCC
# # get values for the variables needed by LCC, phi, phiO, phi1, phi2, lam, lamO, FN, FE, kO, a, inFl.
# # can hard code them in or open a prj file and read them in.
# # for a single parallel Lambert the phi1 and phi2 are the same
# results = LCC.generic_LCC(phi, phiO, phi1, phi2, lam, lamO, FN, FE, kO, a, inFl)
#  # this will create a list [x,y,k,merCon]
#
# *********************************************************************************************

import numpy as np # this is importing numerical tools. numpy needs installed for this to work

# ----------------------------------------------------------
# misc ellipsoid and other functions or methods that might be needed, if these are not used they should be removed
# and maybe just added to a datum script

#def min2Dec(deg,mi,se):
#   h = np.sign(deg)
#   deg = abs(deg)
#   return h * (deg + mi / 60.0 + se / 3600.0)

# phibar is average between two points, for computations. phi is latitude
#def phiBar(phi1,phi2):  
#   return (phi1 + phi2) / 2

# semi-minor axis (b) inFl = 1/f, where f is from ellipsoid definition
#def seMiAx(a,inFl): 
#   return a * (1 - 1 / inFl)

# ellipsoid first eccentricity squared (e**2)
# inFl = 1/f, where f is from ellipsoid definition
def fiEcSq(inFl): 
    return 2 / inFl - inFl ** -2

# ellipsoid second Eccentricity Squared
# e^2 prime, or (e')^2 = fiEcSq(inFl)/ (1 - fiEcSq(inFl))
# inFl = 1/f, where f is from ellipsoid definition
def seEcSq(inFl): 
   return fiEcSq(inFl)/ (1 - fiEcSq(inFl))

# earth radius (Rsube)
#def eaRa(phi): 
#   return seMaAx * (1 - np.sin(np.deg2rad(phi)) ** 2 / inFl)

# meridian radius (Rsubm)
def meRa(phi, a, inFl): 
    return a * (1 - fiEcSq(inFl)) / (1 - fiEcSq(inFl) * np.sin(np.deg2rad(phi)) ** 2) ** 1.5

# prime vertical radius (Rsubn)
def prVeRa(phi, a, inFl): 
    return a / np.sqrt(1 - fiEcSq(inFl) * np.sin(np.deg2rad(phi)) ** 2)

# radius of curvature by azimuth (Rsubalpha)
#def raCur(alpha,phiBar): 
#   return meRa(phiBar) * prVeRa(phiBar) / ((meRa(phiBar) * np.sin(np.deg2rad(alpha)) ** 2) + (prVeRa(phiBar) * np.cos(np.deg2rad(alpha)) ** 2))

# geometric mean ellipsoid radius of curvature (Rg)  
# when h = height above ellipsoid then Rg/(Rg + h) = elevation factor. h is in meters
def geMeRaCu(phi, a, inFl): 
    return np.sqrt(meRa(phi, a, inFl) * prVeRa(phi, a, inFl))



   
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# Lambert Conformal Conic (LCC) supporting functions and definitions
#
# variables and constants
# -----------------------
# phi       geodetic latitude, positive North, in decimal degree
# lam       geodetic longitude, positive East, in decimal degree
# phiO      Central/Standard Parallel - latitude of true projection origin
# phiB      latitude at grid origin, often phiB = phiO 
# phi1,phi2 note: double parallel, phi1 and phi2 should be same for a single
# phiC      phiO = phiC, otherwise PhiC is computed from Phi1 and Phi2. Not used here.
# lamO      Central Meridian, longitude at grid origin
# kO        grid scale factor at Central Parallel
# k         grid scale factor at a general point
# rho       mapping radius for a given phi
# inFl      inverse flatenning, 1/f where f is defined by ellipsoid used
# fiEcSq    first eccentricity squared, from ellipsoid definition
# seEcSq    second eccentricity squared, from ellipsoid definition
# mercon    convergence angle
# m,n,t,F   these are intermediate variables and equations
#
# math is adapted from pages 107-108, "Map Projections - A working Manual", Snyder, 1997 4th printing
#                                     "State Plane Coordinate System of 1983" Stem, James, mar 1990. NOAA Manual NOS NGS 5
#                                     "Ground Truth: Optimized Design of Low Distortion Projections", Dennis, Michael, Mar 2017
#                                     "IOGP Publication 373-7-2 –  Geomatics Guidance Note number 7, part 2" October 2018 www.epsg.org
#
# compute the 'map constants' first, m1, m2, tO, t1, t2, n, F, rhoO then,
# compute the Phi specific m, t, rho, mercon, k, x, and y.
# note that I had to multiply kO to both rho and rhoO, which is not shown in Snyder p.108


# ****** Helper Methods ******
def helperLCC(phi, phiO, phi1, phi2, lam, lamO, kO, a, inFl): # returns rhoO, rho, merCon, and k
    # compute the map constants first. 
    e = np.sqrt(fiEcSq(inFl)) # e, from (e**2)
    m1 = np.cos(phi1) / (1 - fiEcSq(inFl) * np.sin(phi1) ** 2) ** .5
    m2 = np.cos(phi2) / (1 - fiEcSq(inFl) * np.sin(phi2) ** 2) ** .5
    tO = np.tan(np.pi / 4 - phiO / 2) / ((1 - e * np.sin(phiO)) / (1 + e * np.sin(phiO))) ** (e / 2) 
    t1 = np.tan(np.pi / 4 - phi1 / 2) / ((1 - e * np.sin(phi1)) / (1 + e * np.sin(phi1))) ** (e / 2)
    t2 = np.tan(np.pi / 4 - phi2 / 2) / ((1 - e * np.sin(phi2)) / (1 + e * np.sin(phi2))) ** (e / 2)
    if phi1 == phi2 :
        n = np.sin(phiO)     # this works when phi0 is single parallel latitude.
    else :
        n = (np.log(m1) - np.log(m2)) / (np.log(t1) - np.log(t2))
    F = m1/(n * t1 ** n)
    rhoO = a * F * kO * tO ** n   #kO added, not shown in Snyder p.108 (15-7a) 
    # now calculate the specific for phi
    m = np.cos(phi) / (1 - fiEcSq(inFl) * np.sin(phi) ** 2) ** .5
    t = np.tan(np.pi / 4 - phi / 2) / ((1 - e * np.sin(phi)) / (1 + e * np.sin(phi))) ** (e / 2)
    rho = a * F * kO * t ** n     #kO added, not shown in Snyder p.108 (15-7)
    merCon = n* (lam - lamO)      # convergence angle
    k = rho * n / (a * m)
    return rhoO, rho, merCon, k

def generic_LCC(phi, phiO, phi1, phi2, lam, lamO, FN, FE, kO, a, inFl):  #  requires decimal degree
    phi = np.deg2rad(phi)
    phiO = np.deg2rad(phiO)
    phi1 = np.deg2rad(phi1)
    phi2 = np.deg2rad(phi2)
    lam = np.deg2rad(lam)
    lamO = np.deg2rad(lamO)
    helper = helperLCC(phi, phiO, phi1, phi2, lam, lamO, kO, a, inFl) 
    rhoO = helper[0]   # mapping radius at latitude of grid origin
    rho = helper[1]    # mapping radius at phi, R
    merCon = helper[2] # convergence angle, in radians
    x = FE + rho * np.sin(merCon)
    y = FN + rhoO - rho * np.cos(merCon)
    k = helper[3]      # grid scale factor, k, at a point
    return x,y,k,merCon


