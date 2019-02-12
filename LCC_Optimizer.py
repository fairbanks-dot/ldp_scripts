import numpy as np
import akdot_nr_ldp_v1 as ldp

# This is file for finding parameters for projectionsto best fit test points.
#
# Troy Hicks
# last revision: mar 13, 2018
#
# notes: added some points to test points. Also started to write kO optimization routine
#


# global constant for converting meters to us survey feet
usft = 3937.0 / 1200.0  # 1 meter = 3937/1200 US Survey feet

# test points, in decimal degrees, Lat, long, then Ellipsoid Height in meters.
points = np.genfromtxt("Fairbanks design points.csv", delimiter=',') 

print ("Running Lambert Conformal Conic Optimizer on your design points. Standby ...")
################################################################################
################################################################################
#
# Lambert Conformal Conic Projection parameter OPTIMIZER
# ------------------------------------------------------
#     Uses brute force to find the best values for PhiO and kO
#     against supplied test points. Best fit to the points.
#     The linear distortion at each point is computed based on various values of
#     either PhiO or kO. For PhiO, the variance of distortion values for
#     the test points is minimized. For kO the average of the distortion values
#     is minimized. 
#
################################################################################
# Optimize the PhiO parameter
#
#      accomplished by finding the PhiO that has the minimum standard deviation
#      of linear distortion per test point
#
################################################################################
min_PhiO = (np.amin(points[:,0])) - 4.5
max_PhiO = (np.amax(points[:,0])) + 4.5
step_PhiO = 15 # seconds to step for each test
step_PhiO = step_PhiO / 3600 # convert to decimal
tries_PhiO = np.int((max_PhiO-min_PhiO) / step_PhiO +0.5)
# build results array of rows = tries and 2 columns: [test_Phio, stdDev]
results_PhiO = np.zeros([tries_PhiO,2]) 
test_PhiO = min_PhiO   # start the test with the minimum value
best_PhiO = min_PhiO   #
best_stdDev = 99999999 # 
for c in range(tries_PhiO) :
   # for each overall PhiO try, sample the distortion for the test points
   try_PhiO_sample = np.zeros([1,points.shape[0]]) 
   for i in range(points.shape[0]) :
      Rg = ldp.geMeRaCu(points[i][2])
      EF = Rg/(Rg + points[i][2])
      k = ldp.generic_LCC(points[i][0],test_PhiO,1.0)  # needs a (phi, phiO, kO)
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_PhiO_sample[0][i] = d # build the sample array of d values   
   stdDev = np.std(try_PhiO_sample[0,:])
   results_PhiO[c][1] = stdDev  
   results_PhiO[c][0] = test_PhiO  
   if np.abs(best_stdDev) > np.abs(stdDev) :
       best_PhiO = test_PhiO
       best_stdDev = stdDev
   test_PhiO = test_PhiO + step_PhiO   # increase test value to next
print ("PhiO", best_PhiO, "stdDev", best_stdDev)

################################################################################
# Optimize the kO parameter
#
#      accomplished by finding the kO that has the minimum average
#      of linear distortion per test point
#
################################################################################
min_kO = 0.9999
max_kO = 1.0001
step_kO = 0.000001 # step for each test
tries_kO = np.int((max_kO-min_kO) / step_kO +0.5)
test_kO = min_kO  # start the test with the minimum value
best_kO = max_kO  # reset best 
best_avg = 999999999999 # start with high number so can look for lower value
for c in range(tries_kO) :
   # for each overall kO try, sample the distortion for the test points
   try_kO_sample = np.zeros([points.shape[0]]) 
   for i in range(points.shape[0]) :
      Rg = ldp.geMeRaCu(points[i][2])
      EF = Rg/(Rg + points[i][2])
      k = ldp.generic_LCC(points[i][0], best_PhiO, test_kO)
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_kO_sample[i] = d # build the sample array of d values
   avg = np.average(try_kO_sample)
   if np.abs(best_avg) > np.abs(avg) :
       best_kO = test_kO
       best_avg = avg
   test_kO = test_kO + step_kO   # increase test value to next

print ("kO:",best_kO, "  average:", best_avg)



######### write the kO results  - to graph later ##################################
#outfile = open("results_kO_LCC.csv", "w")
#for i in range(results_kO.shape[0]) :
#   print(results_kO[i][0],",", results_kO[i][1],file = outfile)
#outfile.close()
####################################################################################


######### write the LamO results  - to graph later #################################
#outfile = open("results_PhiO_LCC.csv", "w")
#for i in range(results_PhiO.shape[0]) :
#   print(results_PhiO[i][0],",", results_PhiO[i][1],file = outfile)
#outfile.close()
####################################################################################

