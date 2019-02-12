# Distortion writer
#
# reads a point file, then writes out distortion results to another file
# last two lines are standard deviation and then average
#
# *** this is meant to use Defined Projections in imported LDP script ***
#
# Author: Troy Hicks
# last revision: Mar 16, 2018
#
# notes:
#   edit this later on to just get called by another script
#   that will send over the projection type, and the parameters, and even which file to use
#   for points, and what file to write to.
#

import numpy as np
import akdot_nr_ldp_v1 as ldp


# test points, in decimal degrees, Lat, long, then Ellipsoid Height in meters.
points = np.genfromtxt("Fairbanks gridded ellipsoid_A2.csv", delimiter=',')

distortion_results = np.zeros([points.shape[0] + 2]) # last two slots are stdDev and Average
######### write out results of exact LDP already designed ##################
outfile = open("Fairbanks gridded ellipsoid_A2_computed_OM.csv", "w")
print("Lat" ",", "Long", ",", "Ellipsoid", ",", "DISTORTION", file = outfile)
for i in range(points.shape[0]) :
    Rg = ldp.geMeRaCu(points[i][2])
    EF = Rg/(Rg + points[i][2])
#    k = ldp.Fairbanks(points[i][0],points[i][1])  # needs a (phi and lam)
#    k = ldp.computed_LCC(points[i][0],points[i][1])  # needs a (phi and lam)
    k = ldp.generic_OM(points[i][0], 64.94111057676989, points[i][1], -150.00409533805313, -81, 1.000027) # needs phi, PhiC, lam, LamC, Skew, and kO
#    k = ldp.generic_LCC(points[i][0], 64.7137141700027, 1.000025)  # needs a (phi,phiO,kO)
#    d = (EF * k[2] - 1) * 1000000  # d = linear distortion in units of ppm  USE THIS if calling normal function
    d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm  USE THIS if calling generic function
    distortion_results[i] = d # build the sample array of d values   
#    print(distortion_results[i], file = outfile)
    print(points[i][0], ",", points[i][1], ",", points[i][2], ",", distortion_results[i], file = outfile)
outfile.close()
##########################################################################

#for i in range(distortion_results.shape[0]) :
#   print("%.2f" % distortion_results[i])



