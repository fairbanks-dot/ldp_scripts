import numpy as np
import TransverseMercator as tm
import time
start = time.time()

# This is file for finding parameters for projectionsto best fit test points.
#
# Troy Hicks
# last revision: mar 18, 2019
#
# edited to use TransverseMercator.py
# added range and quartile as additional ways to optimize.
#       unknown if one or other is better method. Might remove this later.
#


# global constants
usft = 3937.0 / 1200.0  # 1 meter = 3937/1200 US Survey feet
a = 6378137 #semi-major axis (a) GRS-80, meters
inFl = 298.257222100882 #inverse flattening (1/f) unitless
fiEcSq = 2 / inFl - inFl ** -2 # ellipsoid first eccentricity squared (e**2)

# A few needed ellipsoid functions here in order to compute elevation factor
# ***********************************************************************************
# geometric mean ellipsoid radius of curvature (Rg)  
# when h = height above ellipsoid then Rg/(Rg + h) = elevation factor. h is in meters
def geMeRaCu(phi, a, inFl): 
   return np.sqrt(meRa(phi, a, inFl) * prVeRa(phi, a, inFl))

# meridian radius 
def meRa(phi, a, inFl): 
   return a * (1 - fiEcSq) / (1 - fiEcSq * np.sin(np.deg2rad(phi)) ** 2) ** 1.5                  

# prime vertical radius
def prVeRa(phi, a, inFl): 
   return a / np.sqrt(1 - fiEcSq * np.sin(np.deg2rad(phi)) ** 2)

# Here pull in design points, this should be done by selecting via tkinter
# for now hard code it in.
# test points, in decimal degrees, Lat, long, then Ellipsoid Height in meters.
points = np.genfromtxt("PCCS_3W.csv", delimiter=',')                                                   


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
min_LamO = (np.amin(points[:,1])) - 4.5
max_LamO = (np.amax(points[:,1])) + 4.5
step_LamO = 15 # seconds to step for each test
step_LamO = step_LamO / 3600 # convert to decimal
tries_LamO = np.int((max_LamO-min_LamO) / step_LamO +1)
# build results array of rows = tries and 2 columns: [test_Phio, stdDev]
results_LamO = np.zeros([tries_LamO,2]) 
test_LamO = min_LamO  # start the test with the minimum value
best_LamO = min_LamO  # reset best 
best_stdDev = 999999999999 # start with high number so can look for lower value
best_range = 99999999
best_quartile = 99999999 # inner 50%, (75%-25%) only portion of full range
for c in range(tries_LamO) :
   # for each overall PhiO try, sample the distortion for the test points
   try_LamO_sample = np.zeros([points.shape[0]]) 
   for i in range(points.shape[0]) :
      Rg = geMeRaCu(points[i][2], a, inFl)
      EF = Rg/(Rg + points[i][2])
      TM = tm.compute_TM(points[i][0], points[i][0], points[i][1], test_LamO, 0.0, 0.0, 1.0, a, inFl)        
      k = TM[2]
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_LamO_sample[i] = d # build the sample array of d values   
   k_stdDev = np.std(try_LamO_sample)
   k_range = np.max(try_LamO_sample) - np.min(try_LamO_sample)
   k_quartile = np.percentile(try_LamO_sample, 75) - np.percentile(try_LamO_sample,25)
   if np.abs(best_stdDev) > np.abs(k_stdDev) :
       best_LamO = test_LamO
       best_stdDev = k_stdDev
   if np.abs(best_range) > np.abs(k_range) :
       best_range_LamO = test_LamO
       best_range = k_range
   if np.abs(best_quartile) > np.abs(k_quartile) :
       best_quartile_LamO = test_LamO
       best_quartile = k_quartile
   test_LamO = test_LamO + step_LamO   # increase test value to next

print(" ")
print ("best LamO (stdDev method):  ", "%.8f" % best_LamO,          "stdDev:   ", "%.1f" % best_stdDev)
print ("best LamO (range method):   ", "%.8f" % best_range_LamO,    "range:    ", "%.1f" % best_range)
print ("best LamO (quartile method):", "%.8f" % best_quartile_LamO, "50% range:", "%.1f" % best_quartile)

################################################################################
# Optimize the kO parameter
#
#      accomplished by finding the kO that has the minimum average
#      of linear distortion per test point
#
################################################################################
# check for best scale
#
# start with best using stdDev method for LamO
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
      Rg = geMeRaCu(points[i][2], a, inFl)
      EF = Rg/(Rg + points[i][2])
      TM = tm.compute_TM(points[i][0], points[i][0], points[i][1], best_LamO, 0.0, 0.0, test_kO, a, inFl)            
      k = TM[2]
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_kO_sample[i] = d # build the sample array of d values
   avg = np.average(try_kO_sample) 
   if np.abs(best_avg) > np.abs(avg) :
       best_kO = test_kO
       best_avg = avg
   test_kO = test_kO + step_kO   # increase test value to next
print ("kO (stdDev method):  ","%.6f" % best_kO, "  average:","%.1f" % best_avg)

# now best using range method for LamO
min_kO =best_kO - 0.0001
max_kO = best_kO + 0.0001
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
      Rg = geMeRaCu(points[i][2], a, inFl)
      EF = Rg/(Rg + points[i][2])
      TM = tm.compute_TM(points[i][0], points[i][0], points[i][1], best_range_LamO, 0.0, 0.0, test_kO, a, inFl)            
      k = TM[2]
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_kO_sample[i] = d # build the sample array of d values
   avg = np.average(try_kO_sample) 
   if np.abs(best_avg) > np.abs(avg) :
       best_kO = test_kO
       best_avg = avg
   test_kO = test_kO + step_kO   # increase test value to next
print ("kO (range method):   ","%.6f" % best_kO, "  average:","%.1f" % best_avg)

# now best using quartile method for LamO
min_kO =best_kO - 0.0001
max_kO = best_kO + 0.0001
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
      Rg = geMeRaCu(points[i][2], a, inFl)
      EF = Rg/(Rg + points[i][2])
      TM = tm.compute_TM(points[i][0], points[i][0], points[i][1], best_quartile_LamO, 0.0, 0.0, test_kO, a, inFl)            
      k = TM[2]
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_kO_sample[i] = d # build the sample array of d values
   avg = np.average(try_kO_sample) 
   if np.abs(best_avg) > np.abs(avg) :
       best_kO = test_kO
       best_avg = avg
   test_kO = test_kO + step_kO   # increase test value to next
print ("kO (quartile method):","%.6f" % best_kO, "  average:","%.1f" % best_avg)

# print out some stats on computation time
end = time.time()
total_mins2 = (end - start)/60
num_points = points.shape[0]
seconds_per_point = (end - start) / num_points
print(" ")
print(num_points, " design points used.")
print("%.1f" % total_mins2, " mins total")
print("%.1f" % seconds_per_point, " seconds per design point.")


















