import numpy as np
import akdot_nr_ldp_v1 as ldp

# This is file for finding parameters for projectionsto best fit test points.
#
# Troy Hicks
# last revision: mar 16, 2018
#
# notes: pulling points from a file
#


# global constant for converting meters to us survey feet
usft = 3937.0 / 1200.0  # 1 meter = 3937/1200 US Survey feet

# test points, in decimal degrees, Lat, long, then Ellipsoid Height in meters.
points = np.genfromtxt("Fairbanks design points.csv", delimiter=',')                                                   


print ("Running Transverse Mercator Optimizer on your design points. Standby ...")
################################################################################
################################################################################
#
# Transverse Mercator Projection parameter OPTIMIZER
# ------------------------------------------------------
#     Uses brute force to find the best values for lamO and kO
#     against supplied test points. Best fit to the points.
#     The linear distortion at each point is computed based on various values of
#     either LamO or kO. For LamO, the variance of distortion values for
#     the test points is minimized. For kO the average of the distortion values
#     is minimized. 
#
################################################################################
# Optimize the LamO parameter
#
#      accomplished by finding the LamO that has the minimum standard deviation
#      of linear distortion per test point
#
################################################################################
min_LamO = (np.amin(points[:,1])) -6
max_LamO = (np.amax(points[:,1])) +6
step_LamO = 15 # seconds to step for each test
step_LamO = step_LamO / 3600 # convert to decimal
tries_LamO = np.int((max_LamO-min_LamO) / step_LamO +1)
# build results array of rows = tries and 2 columns: [test_Phio, stdDev]
results_LamO = np.zeros([tries_LamO,2]) 
test_LamO = min_LamO  # start the test with the minimum value
best_LamO = min_LamO  # reset best 
best_stdDev = 999999999999 # start with high number so can look for lower value
for c in range(tries_LamO) :
   # for each overall PhiO try, sample the distortion for the test points
   try_LamO_sample = np.zeros([points.shape[0]]) 
   for i in range(points.shape[0]) :
      Rg = ldp.geMeRaCu(points[i][2])
      EF = Rg/(Rg + points[i][2])
      k = ldp.generic_TM(points[i][0], points[i][1], test_LamO,1.0)  # needs a (phi, lam, lamO, kO)
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_LamO_sample[i] = d # build the sample array of d values   
   stdDev = np.std(try_LamO_sample)
   results_LamO[c][1] = stdDev  
   results_LamO[c][0] = test_LamO 
   if np.abs(best_stdDev) > np.abs(stdDev) :
       best_LamO = test_LamO
       best_stdDev = stdDev
   test_LamO = test_LamO + step_LamO   # increase test value to next

print ("LamO", best_LamO, "stdDev", best_stdDev)


################################################################################
# Optimize the kO parameter
#
#      accomplished by finding the kO that has the minimum average
#      of linear distortion per test point
#
################################################################################
min_kO = 0.999
max_kO = 1.001
step_kO = 0.000001 # step for each test
tries_kO = np.int((max_kO-min_kO) / step_kO +1)
# build results array of rows = tries and 2 columns: [test_kO, avg]
results_kO = np.zeros([tries_kO,2]) 
test_kO = min_kO  # start the test with the minimum value
best_kO = max_kO  # reset best 
best_avg = 999999999999 # start with high number so can look for lower value
for c in range(tries_kO) :
   # for each overall kO try, sample the distortion for the test points
   try_kO_sample = np.zeros([points.shape[0]]) 
   for i in range(points.shape[0]) :
      Rg = ldp.geMeRaCu(points[i][2])
      EF = Rg/(Rg + points[i][2])
      k = ldp.generic_TM(points[i][0], points[i][1], best_LamO, test_kO)   # using best_LamO try a kO
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_kO_sample[i] = d # build the sample array of d values
   avg = np.average(try_kO_sample)
   stdDev2 = np.std(try_kO_sample)

 #  print (try_kO_sample, "avg:", avg, "test_kO", test_kO, "best_LamO", best_LamO, "std", stdDev2)
   results_kO[c][1] = avg  
   results_kO[c][0] = test_kO  
   if np.abs(best_avg) > np.abs(avg) :
       best_kO = test_kO
       best_avg = avg
   test_kO = test_kO + step_kO   # increase test value to next

print ("kO:",best_kO, "  average:", best_avg)


######### write the kO results  - to graph later ##################################
outfile = open("results_kO_TM.csv", "w")
for i in range(results_kO.shape[0]) :
   print(results_kO[i][0],",", results_kO[i][1],file = outfile)
outfile.close()
####################################################################################


######### write the LamO results  - to graph later #################################
outfile = open("results_LamO_TM.csv", "w")
for i in range(results_LamO.shape[0]) :
   print(results_LamO[i][0],",", results_LamO[i][1],file = outfile)
outfile.close()
####################################################################################

















