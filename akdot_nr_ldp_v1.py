# *********************************************************************************************
# *********************************************************************************************
# akdot_nr_ldp_v1.py
# original author   - Mitch Slife, Oct 2017
# current author    - Troy Hicks,  Dec 2017 - Mar 2018  troy.hicks@alaska.gov
#
# **** Most recent revision done on Mar 13,2018
# 
# **** Major edit notes:
#      edited to add content and documentation:
#      added functionality to compute k, grid scale factor at a point. For use in computing linear distortion at a point.
#      some minor edits to improve comments and documentation, and to make easier to read and reference the source of
#      computational methods.
#
# ****** edits: 3/16/2017
#        - added generic LCC and TM projections to use for computing optimized parameters to fit design points
# ******
#
# The PURPOSE of this python script is to perform geodetic computations on various projections
#     in order to gain Grid x and y, grid point scale factor, linear distortion, 
#     and other useful tools to aid in surveying and mapping.
#
# Projection computations and formulas come from the following sources:
#     "Map Projections - A Working Manual", Snyder, John, 1997 printing
#     "State Plane Coordinate System of 1983" Stem, James, mar 1990. NOAA Manual NOS NGS 5
#     "Ground Truth: Optimized Design of Low Distortion Projections", Dennis, Michael, Mar 2017
#     EPSG Guidance Note 7, Coordinate Conversions and Transformations including Formulas. An EPSG circular, http://www.ihsenergy.com/epsg/guid7.pdf
#
# notes: Generally, the actual computations for projections are from Snyder. See comments. Also there is some inverse functions
#        using various methods. A noteworthy one is the Vincenty method.
# *********************************************************************************************

import numpy as np # this is importing numerical tools. numpy needs installed for this to work

# ----------------------------------------------------------
# ----------------------------------------------------------
# Global variables and constants
# ------------------------------
ift = 1 / 0.3048        # International feet
usft = 3937.0 / 1200.0  # 1 meter = 3937/1200 US Survey feet - legal unit in Alaska
# maybe THESE SHOULD BE PART OF A GRS80 DEFINITION under DATUMs, but ok as is
seMaAx = 6378137 #semi-major axis (a) GRS-80, meters
inFl = 298.257222100882 #inverse flattening (1/f) unitless
# ----------------------------------------------------------

def min2Dec(deg,mi,se): # a simple function to convert DMS to decimal
   h = np.sign(deg)
   deg = abs(deg)
   return h * (deg + mi / 60.0 + se / 3600.0)

def phiBar(phi1,phi2):  # phi is latitude, phibar is average between two points, for computations
   return (phi1 + phi2) / 2

def seMiAx(a,f): #semi-minor axis (b)  f= inFl so this function asks for 1/f, not f
   return a * (1 - 1 / f)

def fiEcSq(f): #ellipsoid first eccentricity squared (e**2) f= inFl so this asks for 1/f, not f
   return 2 / f - f ** -2

def eaRa(phi): #earth radius (Rsube)
   return seMaAx * (1 - np.sin(np.deg2rad(phi)) ** 2 / inFl)

def meRa(phi): #meridian radius (Rsubm)
   return seMaAx * (1 - fiEcSq(inFl)) / (1 - fiEcSq(inFl) * np.sin(np.deg2rad(phi)) ** 2) ** 1.5

def prVeRa(phi): #prime vertical radius (Rsubn)
   return seMaAx / np.sqrt(1 - fiEcSq(inFl) * np.sin(np.deg2rad(phi)) ** 2)

def raCur(alpha,phiBar): #radius of curvature by azimuth (Rsubalpha)
   return meRa(phiBar) * prVeRa(phiBar) / ((meRa(phiBar) * np.sin(np.deg2rad(alpha)) ** 2) + (prVeRa(phiBar) * np.cos(np.deg2rad(alpha)) ** 2))

def geMeRaCu(phi): # geometric mean ellipsoid radius of curvature (Rg)  when h = height above ellipsoid then Rg/(Rg + h) = elevation factor. h is in meters
   return np.sqrt(meRa(phi) * prVeRa(phi))


# ------------------------------------------------------------
# compute geodetic coordinates from ECEF coordinates, returns degrees and meters
# this is finicky as it depends on the quadrant - to get the appropiate angle value
# edit the 'adjust for quadrant' section below to get it to work right for your area
def geoECEF(X,Y,Z):
   a   = seMaAx 
   b   = seMiAx(a,inFl)
   e = np.sqrt(fiEcSq(inFl)) # e, from (e**2)
   ee  = e ** 2 / (1 - e **2)
   p   = (X **2 + Y **2) ** 0.5
   q   = np.arctan(Z * a / (p * b))
   phi = np.arctan((Z + ee * b * (np.sin(q)) **3) / (p - (e ** 2) * a * (np.cos(q)) ** 3))
   v   = prVeRa(np.rad2deg(phi)) # prime vertical radius, convert phi to degrees first
   lam = np.arctan (Y / X)
   # now attempt to adjust for the quadrant, and west = negative, east positive
   lam = lam - np.pi # this just work locally ... need to fix this
   h   = p / np.cos(phi) - v
   return  (np.rad2deg(phi), np.rad2deg(lam), h)

# ------------------------------------------------------------
# compute ECEF cartesian cordinates from geodetic
def ecefX(phi,lam,h):
   return (prVeRa(phi) + h) * np.cos(np.deg2rad(phi)) * np.cos(np.deg2rad(lam)) 

def ecefY(phi,lam,h):
   return (prVeRa(phi) + h) * np.cos(np.deg2rad(phi)) * np.sin(np.deg2rad(lam))

def ecefZ(phi,h):
   return (prVeRa(phi) * (1 - fiEcSq(inFl)) + h) * np.sin(np.deg2rad(phi))

# compute inverse using ECEF values.
def ecefDelta(x1,y1,z1,h1,x2,y2,z2,h2): #Cartesian difference (H)    (This is the XYZh INVERSE METHOD and gets about same answer as VINCENTY method (after scaling Vincenty)
   return np.sqrt((x1 - x2)**2 + (y1 - y2)**2 +(z1 - z2)**2 - (h1 - h2)**2)

def ecefCorr(R,H): #local cartesian correction factor, corrects for curvature. only about 0.2ppm until over 20 miles - typically.
   return 2 * R * np.sinh(H / (2 * R * H))

# ------------------------------------------------------------
# misc. geodetic computations
def cenAng(phi1,phi2,lam1,lam2): #central angle between two points on a sphere (psi)
   return np.arccos(np.sin(np.deg2rad(phi1)) * np.sin(np.deg2rad(phi2)) + np.cos(np.deg2rad(phi1)) * np.cos(np.deg2rad(phi2)) * np.cos(np.deg2rad(lam1) - np.deg2rad(lam2))) * 180 / np.pi

def forGeoAzi(phi1,phi2,lam1,lam2): #Approximate forward geodetic azimuth (alpha)
   return np.arctan(np.cos(np.deg2rad(phi2)) * (lam2 - lam1) / (phi2 - phi1)) * 180 / np.pi 

def geoInv(phi1,lam1,phi2,lam2): #Approximate geodetic inverse based on spherical angle (psi)
   return raCur(forGeoAzi(phi1,phi2,lam1,lam2),phiBar(phi1,phi2)) * cenAng(phi1,phi2,lam1,lam2) * np.pi / 180


# -------------------------------------------------------------
# Vincenty method
# computes s, the ellipsoid distance between two geodetic points. Can later scale up to topograpghy to get ground distance.
# This algorithym typically is accurate to about 0.1mm. There are other numerical methods that are good to a 0.001mm.
# THis was initially being used to get ground distance to compare to grid distances. To check projection performance.
# However has been replaced by simply computing linear distortion at a point directly. (CSF-1)*1,000,000 as parts per million.

# Vincenty definitions
def phiAux(phi): #U1,2
   return (np.arctan((1 - 1 / inFl) * np.tan(np.deg2rad(phi))))

def sinSigma(U1,U2,lam):
   return np.sqrt((np.cos(U2) * np.sin(lam)) ** 2 + (np.cos(U1) * np.sin(U2) - np.sin(U1) * np.cos(U2) * np.cos(lam)) ** 2)

def cosSigma(U1,U2,lam):
   return np.sin((U1)) * np.sin((U2)) + np.cos((U1)) * np.cos((U2)) * np.cos(lam)

def sigma(sinSigma,cosSigma):
   return np.arctan(sinSigma / cosSigma)

def sinAlpha(U1,U2,lam,sinSigma):
   return np.cos((U1)) * np.cos((U2)) * np.sin(lam) / sinSigma

def cosSqrAlpha(sinAlpha):
   return 1 - (sinAlpha) ** 2

def cosTwoSigma(cosSigma,U1,U2,cosSqrAlpha):
   return cosSigma - 2 * np.sin((U1)) * np.sin((U2)) / cosSqrAlpha

def C(cosSqrAlpha):
   return cosSqrAlpha * (4 + (4 - 3 * cosSqrAlpha) / inFl) / (16 * inFl)

def lambd(L,C,sinAlpha,sigma,cosTwoSigma):
   return L + (1 - C) * sinAlpha * (sigma + C * np.sin(sigma) * (cosTwoSigma + C * np.cos(sigma) * (2 * cosTwoSigma ** 2 - 1))) / inFl

# Vincenty iterative process
def vincenty(phi1,lam1,phi2,lam2):   
   U1 = phiAux(phi1)
   U2 = phiAux(phi2)
   L = np.deg2rad(lam2 - lam1)
   lam = L
   delLam = 1
   while delLam > float("1.0e-12"):
      sinSig = sinSigma(U1,U2,lam)
      cosSig = cosSigma(U1,U2,lam)
      sig = sigma(sinSig,cosSig)
      sinA = sinAlpha(U1,U2,lam,sinSig)
      cos2A = cosSqrAlpha(sinA)
      cos2Sig = cosTwoSigma(cosSig,U1,U2,cos2A)
      c = C(cos2A)
      lamNew = lambd(L,c,sinA,sig,cos2Sig)
      delLam = abs(lam - lamNew)
      lam = lamNew
   uSqr = cos2A * (seMaAx ** 2 - seMiAx(seMaAx,inFl) ** 2) / seMiAx(seMaAx,inFl) ** 2
   A = 1 + uSqr * (4096 + uSqr * (-768 + uSqr * (320 - 175 * uSqr))) / 16384
   B = uSqr * (256 + uSqr * (-128 + uSqr * (74 - 47 * uSqr))) / 1024
   delSig = B * sinSig * (cos2Sig + B * (cosSig * (2 * cos2Sig ** 2 - 1) - B * cos2Sig * (4 * sinSig ** 2 - 3) * (4 * cos2Sig ** 2 - 3)) / 4)
   s = seMiAx(seMaAx,inFl) * A * (sig - delSig)
   return s

def ellCor(phi1,h1,phi2,h2): #Ellipsoidal distance correction in meters
   phiBar = (phi1 + phi2) / 2
   hBar = (h1 + h2) / 2
   Rg = geMeRaCu(phiBar)
   return (Rg + hBar) / Rg
   
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# Lambert Conformal Conic (LCC) supporting functions and definitions
#           - these are used by each projection function to compute x (easting), y (northing)
#
# variables and constants
# -----------------------
# phi       geodetic latitude, positive North
# lam       geodetic longitude, positive West (is this right??)  : we use positive East, this may alrady work fine.
# phiO      Central/Standard Parallel - latitude of true projection origin
# phiB      latitude at grid origin, often phiB = phi0 
#           ***************************************************************************************
# phi1,phi2 note: double parallel defintions need reduced to single, otherwise phi1 and phi2 should
# phiC      be made and used to compute phiC, then let phi0 = phiC
#           ***************************************************************************************
# lamO      Central Meridian, longitude at grid origin
# kO        grid scale factor at Standard Parallel
# k         grid scale factor at a general point
# rhoLCC    function that returns mapping radius for a given phi, R (called Rho by Mitch Slife)
# merConLCC function that returns convergence angle, often called gamma, for a given lam

def rhoLCC(phi,phiO,k): # returns rho, mapping radius at latitude phi
   # seMaAx = semi-major axis (a) GRS-80, meters      : global variable
   # inFl   = inverse flattening (1/f) GRS-80, meters : global variable
   # fiEcSq = function that computes ellipsoid first eccentricity squared (e**2)
   a = seMaAx
   n = np.sin(phiO)     # this works because phi0 is single parallel latitude, pre-compute phi0 otherwise.
   e = np.sqrt(fiEcSq(inFl)) # e, from (e**2)
   t = np.tan(np.pi / 4 - phi / 2) / ((1 - e * np.sin(phi)) / (1 + e * np.sin(phi))) ** (e / 2)
   tO = np.tan(np.pi / 4 - phiO / 2) / ((1 - e * np.sin(phiO)) / (1 + e * np.sin(phiO))) ** (e / 2)
   m = np.cos(phiO) / (1 - fiEcSq(inFl) * np.sin(phiO) ** 2) ** .5 
   F = m/(n * tO ** n)
   return a * F * k * t ** n

def merConLCC(lam,lamO,phiO): # returns convergence angle, needed to compute Northing/Easting
   return np.sin(phiO) * (lam - lamO)

def xLCC(FE,rho,lam,lamO,phiO): # returns Easting grid value, us survey ft
   return FE + rho * np.sin(merConLCC(lam,lamO,phiO))

def yLCC(FN,rho,rhoO,lam,lamO,phiO): # returns Northing grid value, us survey ft *** should be using rhoB, not rho0
   return FN + rhoO - rho * np.cos(merConLCC(lam,lamO,phiO))

def kLCC(phi, phiO, rho): # returns grid scale factor, k, at a point.
   # fiEcSq = function that computes ellipsoid first eccentricity squared (e**2)
   e = np.sqrt(fiEcSq(inFl))  # e, from (e**2) **** yes this might be silly to do it this way ***
   a = seMaAx  
   return np.sqrt(1-e**2 * (np.sin(phi)**2)) * (rho * np.sin(phiO) / (a * np.cos(phi)))

# ----------------------------------------------------------------------------------------------
# Projection defintions that use LCC. most of these are Low Distortion Projections developed by AK DOT
# this area can also have published LCC projection defintions such as State Plane Zones

def generic_LCC(phi,phiO,kO):  # generic, requires phi in decimal degree
   phi = np.deg2rad(phi)
   #lam = np.deg2rad(lam)
   phiO = np.deg2rad(phiO)
   #lamO = np.deg2rad(min2Dec(-156,44,0.0))
   #FN = 0.0 
   #FE = 0.0 
   #kO = 0.99999
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   #x = xLCC(FE,rho,lam,lamO,phiO)
   #y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return k

def computed_LCC(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(64.71371417) # standard parallel
#   phiO = np.deg2rad(min2Dec(64,35,3.774)) # standard parallel
   lamO = np.deg2rad(min2Dec(-146,56,0.0))
   FN = 200000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 800000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000025
   rho = rhoLCC(phi,phiO,kO) # mapping radius at phi, R
   rhoO = rhoLCC(phiO,phiO,kO) # mapping radius at latitude of grid origin, R0, assuming phiB = phi0
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def computed_LCC2(phi,lam):   # Fairbanks Baseline
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(64,50,56.36374)) # standard parallel
   lamO = np.deg2rad(min2Dec(-147,48,0.0))
   FN = 200000.0    # convert this value, which is in usft, to meters.
   FE = 800000.0    # convert this value, which is in usft, to meters.
   kO = 1.000022329
   rho = rhoLCC(phi,phiO,kO) # mapping radius at phi, R
   rhoO = rhoLCC(phiO,phiO,kO) # mapping radius at latitude of grid origin, R0, assuming phiB = phi0
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Fairbanks(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(64,51,0.0)) # standard parallel
#  phiB = np.deg2rad(64.85) # would need to have the phi, phiB, for grid origin to compute rhoB - called rho0 in this case, when phi0 not equal to phiB
   lamO = np.deg2rad(min2Dec(-146,56,0.0))
   FN = 200000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 800000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.00003
   rho = rhoLCC(phi,phiO,kO) # mapping radius at phi, R
   rhoO = rhoLCC(phiO,phiO,kO) # mapping radius at latitude of grid origin, R0, assuming phiB = phi0
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def PCCS_4(phi,lam):    # Pima County Coordinate System 4, Mt Lemmon.  
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(30,30,00.0))
   lamO = np.deg2rad(min2Dec(-110,45,0.0))
   FN = -620000.0 ** ift ** -1   # convert this value, which is in ift, to meters.
   FE = 30000 * ift ** -1        # convert this value, which is in ift, to meters.
   kO = 0.999800
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Barrow(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(71,2,00.0))
   lamO = np.deg2rad(min2Dec(-156,44,0.0))
   FN = 0.0 ** usft ** -1   # convert this value, which is in usft, to meters.
   FE = 30000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 0.99999
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Cordova(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(60,30,0.0))
   lamO = np.deg2rad(min2Dec(-145,15,0.0))
   FN = 100000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 300000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000004
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Dalton_Hwy_Zone_2(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(66,56,0.0))
   lamO = np.deg2rad(min2Dec(-149,38,0.0))
   FN = 100000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 50000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 0.999952
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Denali_Hwy(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(63,25,0.0))
   lamO = np.deg2rad(min2Dec(-149,16,0.0))
   FN = 160000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 0.0 * usft * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000122
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Galena(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(64,44,0.0))
   lamO = np.deg2rad(min2Dec(-158,19,0.0))
   FN = 116000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000007
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Glenn_Hwy_Zone_1(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(62,18,0.0))
   lamO = np.deg2rad(min2Dec(-146,59,23.1189))
   FN = 150000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000089
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Kotzebue(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(66,53,0.0))
   lamO = np.deg2rad(min2Dec(-162,39,0.0))
   FN = 92000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.0
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Nome(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(64,47,0.0))
   lamO = np.deg2rad(min2Dec(-166,26,0.0))
   FN = 127000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 0.999992
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Parks_Hwy_Zone_2(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(64,2,00.0))
   lamO = np.deg2rad(min2Dec(-149,25,0.0))
   FN = 300000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 600000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.00006
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Parks_Hwy_Zone_3(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(62,42,0.0))
   lamO = np.deg2rad(min2Dec(-150,27,0.0))
   FN = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 500000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.00006
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def St_Lawrence_Is(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(63,23,0.0))
   lamO = np.deg2rad(min2Dec(-171,39,0.0))
   FN = 600000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 50000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 0.999987
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

def Tanana_Manley(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(65,9,00.0))
   lamO = np.deg2rad(min2Dec(-152,35,0.0))
   FN = 62000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000023
   rho = rhoLCC(phi,phiO,kO)
   rhoO = rhoLCC(phiO,phiO,kO)
   x = xLCC(FE,rho,lam,lamO,phiO)
   y = yLCC(FN,rho,rhoO,lam,lamO,phiO)
   k = kLCC(phi, phiO, rho)
   return x,y,k

# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# Transverse Mercator (TM) supporting functions and definitions

# fiEcSq(inFl)/ (1 - fiEcSq(inFl)) = e^2 prime, (e')^2 : this is used below.

def nTMerc(phi): # N, page 61 Snyder
   return seMaAx / np.sqrt(1 - fiEcSq(inFl) * np.sin(phi) ** 2)

def tTMerc(phi): # T, page 61 Snyder
   return np.tan(phi) ** 2

def cTMerc(phi): # C, page 61 Snder, fiEcSq(inFl)/ (1 - fiEcSq(inFl)) = e^2 prime, (e')^2
   return fiEcSq(inFl)/ (1 - fiEcSq(inFl)) * np.cos(phi) ** 2 

def aTMerc(phi,lam,lamO): # A, page 61 Snyder - *** I changed the order, to change sign.
   return (lam - lamO) * np.cos(phi) 

def mTMerc(phi): # not same, but looks similar, to M on page 61 Snyder
   e = fiEcSq(inFl) # this is e^2, not actually e. eccentricity squared. therefore next equation reflects e as actually being e^2
   return seMaAx  * ((1 - e / 4 - 3 * e ** 2 / 64 - 5 * e ** 3 / 256) * phi - (3 * e / 8 + 3 * e ** 2 / 32 + 45 * e ** 3 / 1024) * np.sin(2 * phi) + (15 * e ** 2 / 256 + 45 * e ** 3 / 1024) * np.sin(4 * phi) - (35 * e ** 3 / 3072) * np.sin(6 * phi)) 

def xTMerc(FE,kO,N,A,T,C):
   ep = fiEcSq(inFl) / (1 - fiEcSq(inFl))
   return FE + kO * N * (A + (1 - T + C) * A ** 3 / 6 + (5 - 18 * T + T ** 2 + 72 * C - 58 * ep) * A ** 5 / 120)

def yTMerc(FN,kO,M,MO,N,A,T,C,phi):
   ep = fiEcSq(inFl) / (1 - fiEcSq(inFl))
   return FN + kO * (M - MO + N * np.tan(phi) * (A ** 2 / 2 + (5 - T + 9 * C + 4 * C ** 2) * A ** 4 / 24 + (61 - 58 * T + T ** 2 + 600 * C - 330 * ep) * A ** 6 / 720)) 

def kTMerc(kO,A,T,C): # compute k, grid scale factor at a point
   return kO * (1 + (1 + C) * A ** 2 / 2 + (5 - 4 * T + 42 * C + 13* C ** 2 - 28 * (fiEcSq(inFl)/ (1 - fiEcSq(inFl)))) * A ** 4 / 24 + (61 - 148* T + 16 * T **2) * A ** 6 /720)

# ----------------------------------------------------------------------------------------------
# Projection defintions that use TM. most of these are Low Distortion Projections developed by AK DOT
# this area can also have published TM projection defintions such as State Plane Zones

def generic_TM(phi, lam, lamO, kO): # generic, requires phi, lam, lamO in decimal degree
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
#   phiO = np.deg2rad(min2Dec(54,0,0))
   lamO = np.deg2rad(lamO)
#   FE = 500000.0 # value is in meters
#   FN = 0.0 # value is in meters
#   kO = 0.9999
#   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
#   M = mTMerc(phi)
#   MO = mTMerc(phiO)
#   x = xTMerc(FE,kO,N,A,T,C)
#   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return k

# *****************************************************************************
def computed_TM(phi,lam): # temporary test TM, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(31,30,00.0))
   lamO = np.deg2rad(-113.31338723267214)
   FN = 300000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 600000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000039
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k
# *****************************************************************************

def mod_banner(phi,lam): # Mod Banner Creek for a test, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(54,00,00.0))
   lamO = np.deg2rad(-146.0)
   FN = -533.0215 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 1640424.6519 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.00004173
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def mod_AKSPZ3(phi,lam): # Fairbanks Baseline test TM, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(54,00,00.0))
   lamO = np.deg2rad(-146.0)
   FN = -40.4946
   FE = 500002.8516
   kO = 0.99993347111
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def PCCS_3W(phi,lam): # PCCS 3W, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(31,30,0))
   lamO = np.deg2rad(min2Dec(-113,10,0))
   FE = 600000.0 * ift **-1 # convert value to meters
   FN = 0.0 # value is in meters
   kO = 1.000055
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def PCCS_2C(phi,lam): # PCCS 2C, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(31,15,0))
   lamO = np.deg2rad(min2Dec(-112,10,0))
   FE = 1800000.0 * ift **-1 # convert value to meters
   FN = 1000000.0 # value is in meters
   kO = 1.000090
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def AKSPZ2(phi,lam): # Alaska State Plane Zone 2, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(54,0,0))
   lamO = np.deg2rad(min2Dec(-142,0,0))
   FE = 500000 # value is in meters
   FN = 0 # value is in meters
   kO = 0.9999
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k


def AKSPZ3(phi,lam): # Alaska State Plane Zone 3, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(54,0,0))
   lamO = np.deg2rad(min2Dec(-146,0,0))
   FE = 500000 # value is in meters
   FN = 0 # value is in meters
   kO = 0.9999
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def AKSPZ4(phi,lam): # Alaska State Plane Zone 4, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(54,0,0))
   lamO = np.deg2rad(min2Dec(-150,0,0))
   FE = 500000.0 # value is in meters
   FN = 0.0 # value is in meters
   kO = 0.9999
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def AKSPZ7(phi,lam): # Alaska State Plane Zone 7, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(54,0,0))
   lamO = np.deg2rad(min2Dec(-162,0,0))
   FE = 500000.0 # value is in meters
   FN = 0.0 # value is in meters
   kO = 0.9999
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def UTM6N(phi,lam): # UTM Zone 6 North, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(0,0,0))
   lamO = np.deg2rad(min2Dec(-146,0,0))
   FE = 500000.0 # value is in meters
   FN = 0.0 # value is in meters
   kO = 0.9996
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def UTM4N(phi,lam): # UTM Zone 4 North, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(0,0,0))
   lamO = np.deg2rad(min2Dec(-159,0,0))
   FE = 500000.0 # value is in meters
   FN = 0.0 # value is in meters
   kO = 0.9996
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def rich3(phi,lam): # Rich Zone 3 LDP, returns values in meters
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(61,7,0))
   lamO = np.deg2rad(min2Dec(-144,46,0))
   FE = 1000000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FN = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000071
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def Anchorage(phi,lam):  # Anchorage 2015 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(61,30,0.0))
   lamO = np.deg2rad(min2Dec(-149,35,0.0))
   FN = 550000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 250000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000013
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def Dalton_Hwy_Zone_1(phi,lam): # Dalton Zone 1 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(68,8,00.0))
   lamO = np.deg2rad(min2Dec(-145,47,0.0))
   FN = 300000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 600000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 0.999872   
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k


def Dalton_Hwy_Zone_3(phi,lam): # Dalton Zone 3 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(66,4,00.0))
   lamO = np.deg2rad(min2Dec(-150,14,0.0))
   FN = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 900000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.00006
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def Parks_Hwy_Zone_1(phi,lam): # Parks Zone 1 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(64,11,0.0))
   lamO = np.deg2rad(min2Dec(-148,37,0.0))
   FN = 450000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 300000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.00002
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k

def Richardson_Hwy_Zone_3(phi,lam): # Rich Zone 3 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(61,7,00.0))
   lamO = np.deg2rad(min2Dec(-144,46,0.0))
   FN = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 1000000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000071
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k


def Valdez(phi,lam): # Valdez LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiO = np.deg2rad(min2Dec(60,54,00.0))
   lamO = np.deg2rad(min2Dec(-147,30,0.0))
   FN = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 0.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 0.99995
   N = nTMerc(phi)
   T = tTMerc(phi)
   C = cTMerc(phi)
   A = aTMerc(phi,lam,lamO)
   M = mTMerc(phi)
   MO = mTMerc(phiO)
   x = xTMerc(FE,kO,N,A,T,C)
   y = yTMerc(FN,kO,M,MO,N,A,T,C,phi)
   k = kTMerc(kO,A,T,C)                                           
   return x,y,k


# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# Oblique Mercator (OM) - Rectified Skew Orthomorphic - functions and definitions

# suporting functions
#
# (fiEcSq(inFl)) is a global function used to comute e**2, or e^2. Maybe 'e' is a reserved global variable, I dont think so. This should just be computed at beginning of this script.
# the e being used is just coming from GRS80. In the future the ellipsoid being used can be changed, would not afect this area of code.
# seMaAx is a global variable for 'a', from GRS80. Again, this can be changed above for a different ellipsoid in the future.

def bOMerc(phiC): # B, page 74 Snyder 
   return np.sqrt(1 + (fiEcSq(inFl)) * np.cos(phiC) ** 4 / (1 - fiEcSq(inFl)))

def aOMerc(kO,B,phiC): # A, page 74 Snyder
   return  seMaAx * B * kO * np.sqrt(1 - (fiEcSq(inFl))) / (1 - (fiEcSq(inFl)) * np.sin(phiC) ** 2) 

def tOMerc(phiC): # tO, page 74 Snyder later we run similar function, but it uses phi and not phiC
   return np.tan(np.pi / 4 - phiC / 2) / ((1 - np.sqrt(fiEcSq(inFl)) * np.sin(phiC)) / (1 + np.sqrt(fiEcSq(inFl)) * np.sin(phiC))) ** (np.sqrt(fiEcSq(inFl)) / 2)

def dOMerc(B,phiC): # D, page 74 Snyder
   return B * np.sqrt(1 - (fiEcSq(inFl))) / (np.cos(phiC) * np.sqrt(1 - fiEcSq(inFl) * np.sin(phiC) ** 2))

def fOMerc(D): # F, page 74 Snyder 
   return D + np.sqrt(D ** 2 - 1)  # would be D - (minus) if phiC was negative. Alaska is positive latitude

def eOMerc(F,tO,B): # E, page 74 Snyder     
   return F * tO ** B

def gOMerc(F): # G, page 74 Snyder
   return (F - 1 / F) / 2

def gamOMerc(skew,D): # gammaO, page 74 Snyder
   return np.arcsin(np.sin(skew) / D)

def lamOMerc(lamC,G,gamO,B): # computes lamO, using lamC, the two are not the same. page 74 Snyder
   return lamC - (np.arcsin(G * np.tan(gamO))) / B 

def uCOMerc(A,B,D,skew):
   return (A / B) * np.arctan(np.sqrt(D ** 2 - 1) / np.cos(skew))

def tMerc(phi): # t, page 72 Snyder eariler we an similar function, but it uses phi and not phi
   return np.tan(np.pi / 4 - phi / 2) / ((1 - np.sqrt(fiEcSq(inFl)) * np.sin(phi)) / (1 + np.sqrt(fiEcSq(inFl)) * np.sin(phi))) ** (np.sqrt(fiEcSq(inFl)) / 2)

def qOMerc(E,t,B): # Q, page 72 Snyder, (equation 9-25)
   return E / t ** B

def sOMerc(Q): # S, page 72 Snyder, (equation 9-26)
   return (Q - 1 / Q) / 2

def TOMerc(Q): # T, page 72 Snyder, (equation 9-27)
   return (Q + 1 / Q) / 2

def VOMerc(B,lam,lamO): # V, page 72 Snyder, (equation 9-28)
   return np.sin(B * (lam - lamO))

def UOMerc(V,gamO,S,T): # U, page 72 Snyder, (equation 9-29)
   return (-V * np.cos(gamO) + S * np.sin(gamO)) / T

def vOMerc(A,U,B): # v, page 72 Snyder, (equation 9-30)
   return A * np.log((1 - U) / (1 + U)) / (2 * B) # in Python, log is the natural log. log10 would be the 10 base, etc. 

def uOMerc(A,S,gamO,B,V,lam,lamO): # u, page 72 Snyder, (equation 9-31)
   return A * np.arctan((S * np.cos(gamO) + V * np.sin(gamO)) / (np.cos(B * (lam - lamO)))) / B

def xOMerc(v,u,skew,FE):
   return v * np.cos(skew) + u * np.sin(skew) + FE

def yOMerc(v,u,skew,FN):
   return u * np.cos(skew) - v * np.sin(skew) + FN

def kOMerc(A,B,u,phi,lam,lamO): # k, page 72 Snyder, (equation 9-32). k = the scale factor at a point. Can be used with elevation factor to compute linear distortion at a point.
   e = np.sqrt(fiEcSq(inFl))  # e, from (e**2) 
   a = seMaAx # in meters
   k = ((A * np.cos(B * u / A)) * (1 - e**2 * np.sin(phi)**2)**0.5 / ((a * np.cos(phi)) * (np.cos(B * (lam - lamO)))))
   return k
  

# ----------------------------------------------------------------------------------------------
# Projection definitions that use OM. most of these are Low Distortion Projections developed by AK DOT
# this area can also have published OM projection defintions such as State Plane Zones
#
# phiC, lamC = lat and long of the selected center of the map, falling on the central line. Comes from
#              the projection definition. NOTE: phiC is called phiO in page 74 Snyder
# skew = angle of azimuth east of north, for the central line. Note: called 'alpha' page 74 Snyder
#
# each projection defintion below is given as a function, and calls the other supporting supporting functions as needed.
# 

def generic_OM(phi, PhiC, lam, LamC, Skew, kO): # generic, requires phi, PhiC, lam, LamC, Skew in decimal degree, and kO
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(PhiC)  # this is also phiO, in this case phiC = phiO both are same thing.
   lamC = np.deg2rad(LamC) # lamC is the longitude at the origin, part of definition 
   skew = np.deg2rad(Skew) 
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
#   x = xOMerc(v,u,skew,FE)
#   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return k

def computed_OM(phi,lam): # test computed OM, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(33.49541751951218)  # this is also phiO, in this case phiC = phiO both are same thing.
   lamC = np.deg2rad(-109.34292579999995) # lamC is the longitude at the origin, part of definition
   FN = -22400000.0  * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 5700000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000108 
   skew = np.deg2rad(55) 
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k

def PCCS_1E(phi,lam): # Dalton Zone 4 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(32,15,0.0))
   lamC = np.deg2rad(min2Dec(-111,24,0.0))
   kO = 1.000110
   FE = 160000.0 * ift ** -1   # convert this value, which is in ift, to meters.
   FN = 800000.0 * ift ** -1   # convert this value, which is in ift, to meters.
   skew = np.deg2rad(45.0)
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k

def Richardson_Hwy_Zone_2(phi,lam): # Richardson Zone 2 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(62,59,0.0))
   lamC = np.deg2rad(min2Dec(-146,44,0.0))
   kO = 1.000084
   FE = 5700000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FN = -22400000.0 * usft ** -1   # convert this value, which is in usft, to meters..
   skew = np.deg2rad(-13)
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k

def dalt4(phi,lam): # Dalton Zone 4 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(65,45,0.0))
   lamC = np.deg2rad(min2Dec(-149,27,0.0))
   kO = 1.000045
   FE = 20300000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FN = -17900000.0 * usft ** -1   # convert this value, which is in usft, to meters..
   skew = np.deg2rad(min2Dec(-47,0,0))
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k


def Alaska_Hwy(phi,lam): # Alaska LDP, returns values in meters. This is for the Alaska Canadian Highway.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(63,2,00.0))
   lamC = np.deg2rad(min2Dec(-143,52,0.0))
   FN = -16200000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 22000000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000046
   skew = np.deg2rad(min2Dec(-52,0,00.0))
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k


def Dalton_Hwy_Zone_4(phi,lam): # Dalton Zone 4 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(65,45,0.0))
   lamC = np.deg2rad(min2Dec(-149,27,0.0))
   FN = -17900000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 20300000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000045
   skew = np.deg2rad(min2Dec(-47,0,00.0))
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   t = tOMerc(phi)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)
   uC = uCOMerc(A,B,D,skew)
   Q = qOMerc(E,t,B)
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam, lamC)
   return x,y,k

def Glenn_Hwy_Zone_2(phi,lam): # Glenn Zone 2 LDP, returns values in meters.
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(62,57,0.0))
   lamC = np.deg2rad(min2Dec(-148,49,0.0))
   FN = -18165000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = -17479000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 0.999904
   skew = np.deg2rad(min2Dec(43,0,00.0))
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k


def Richardson_Hwy_Zone_1(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(64,18,0.0))
   lamC = np.deg2rad(min2Dec(-145,43,0.0))
   FN = -13700000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = 25100000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000047
   skew = np.deg2rad(min2Dec(-60,0,00.0))
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k

  
def Taylor_Hwy(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(63,57,0.0))
   lamC = np.deg2rad(min2Dec(-142,22,0.0))
   FN = -20000000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = -8000000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.00012
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k

def Tok_Cutoff(phi,lam):
   phi = np.deg2rad(phi)
   lam = np.deg2rad(lam)
   phiC = np.deg2rad(min2Dec(63,5,00.0))
   lamC = np.deg2rad(min2Dec(-143,59,0.0))
   FN = -18500000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   FE = -16700000.0 * usft ** -1   # convert this value, which is in usft, to meters.
   kO = 1.000091
   skew = np.deg2rad(min2Dec(42,0,00.0))
   B = bOMerc(phiC)
   A = aOMerc(kO,B,phiC)
   tO = tOMerc(phiC)
   D = dOMerc(B,phiC)
   F = fOMerc(D)
   E = eOMerc(F,tO,B)
   G = gOMerc(F)
   gamO = gamOMerc(skew,D)
   lamO = lamOMerc(lamC,G,gamO,B)   # this is NOT same as lamC, lamC is longitude at center of zone definition.
   uC = uCOMerc(A,B,D,skew)
   t = tMerc(phi)
   Q = qOMerc(E,t,B) 
   S = sOMerc(Q)
   T = TOMerc(Q)
   V = VOMerc(B,lam,lamO)
   U = UOMerc(V,gamO,S,T)
   v = vOMerc(A,U,B)
   u = uOMerc(A,S,gamO,B,V,lam,lamO)
   x = xOMerc(v,u,skew,FE)
   y = yOMerc(v,u,skew,FN)
   k = kOMerc(A,B,u,phi,lam,lamO)
   return x,y,k

# End of specific projection definitions
# ------------------------------------------------------------------------------------------------------------

