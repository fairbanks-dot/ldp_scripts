import numpy as np
import akdot_nr_ldp_v1 as ldp
import time
start = time.time()
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



################################################################################
################################################################################
#
# Oblique Mercator Projection parameter OPTIMIZER
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
#min_PhiO = (np.amin(points[:,0])) -3
#max_PhiO = (np.amax(points[:,0])) +3
step_minutes = 6  # number minutes to step each time

range_PhiO = 3.5
avg_PhiO = (np.average(points[:,0]))
step_PhiO = step_minutes*60 # seconds to step for each test
step_PhiO = step_PhiO / 3600.0 # convert to decimal

range_LamO = 3.5
avg_LamO = (np.average(points[:,1]))
step_LamO = step_minutes*60 # seconds to step for each test
step_LamO = step_LamO / 3600 # convert to decimal

step_Skew = 2          # set step for Skew (in degrees), initial. later will look each degree
best_LamO = avg_LamO   # reset best 
best_PhiO = avg_PhiO   # reset best
best_Skew = -90        # reset best
best_stdDev = 99999999 # start with high number so can look for lower value
phase1 = False  # phase1 is not complete. During phase1 we use a very coarse search, and stop once it seems no longer gets better.


for c in range (2) :
   test_PhiO = avg_PhiO - range_PhiO  # reset test_PhiO
   while test_PhiO <= avg_PhiO + range_PhiO :     
      test_LamO = avg_LamO - range_LamO #reset test_LamO
      while test_LamO <= avg_LamO + range_LamO:  
         test_Skew = -89              # reset test_Skew
         while test_Skew <= 89 :      
            distortion_sample = np.zeros([points.shape[0]])  # reset the distortion sample for points
            for i in range(points.shape[0]) :
               Rg = ldp.geMeRaCu(points[i][2])
               EF = Rg/(Rg + points[i][2])
               k = ldp.generic_OM(points[i][0], test_PhiO, points[i][1], test_LamO, test_Skew, 1.0)  # needs a (phi, phiC,Lam, LamC, Skew, kO)
               d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
               distortion_sample[i] = d # build the sample array of d values   
            stdDev = np.std(distortion_sample)
            if np.abs(best_stdDev) > np.abs(stdDev) :
               best_LamO = test_LamO
               best_PhiO = test_PhiO
               best_Skew = test_Skew
               best_stdDev = stdDev
               print(" ")
               print("********************************************************************************************************************")
               print("******** WINNER! ************ PhiO", best_PhiO, "LamO", best_LamO, "Skew", best_Skew, "stdDev", "%.5f" % best_stdDev)
               print("********************************************************************************************************************")
               print(" ")
            # this is a check to see if solutions are simply getting progressively worse and therefore lets stop searching.
            if stdDev > 20*best_stdDev and test_Skew == best_Skew and phase1 != True :
               print(" ")
               print("phase 1 complete, begin phase 2")
               phase1 = True
               avg_LamO = best_LamO
               range_LamO = step_minutes/60/2
               avg_PhiO = best_PhiO
               range_PhiO = step_minutes/60/2
               step_LamO = step_LamO / step_minutes  # to now look at each minute interval
               step_PhiO = step_PhiO / step_minutes  # to now look at each minute interval
               step_Skew = 1
               test_PhiO = avg_PhiO + range_PhiO + 1 # stop checking PhiO for now
               test_LamO = avg_LamO + range_LamO + 1 # stop checking LamO for now
               test_Skew = 91                        # stop checking Skew for now
            if stdDev > 1.4*best_stdDev and test_Skew == best_Skew and phase1 == True and c == 1:
               print(" ")
               print("phase 2 may not find a better answer. ", stdDev, "PhiO spread", avg_PhiO-test_PhiO, " LamO spread", avg_LamO-test_LamO)              
            test_Skew = test_Skew + step_Skew
         print(".", end="")
         test_LamO = test_LamO + step_LamO
      print(" ")
      #percent_complete = (avg_PhiO - range_PhiO - test_PhiO)/(2 * range_PhiO) * 100
      current = time.time()
      time_so_far = np.int((current-start)*100)
      time_so_far = time_so_far/100
      current_mins = (time_so_far)/60
      #total_mins = ((percent_complete/100) **-1) * (time_so_far)/60
      print ("%.1f" % current_mins, "minutes so far")
      test_PhiO = test_PhiO + step_PhiO
 
print ("PhiO:", best_PhiO, "LamO", best_LamO, "Skew:", best_Skew)
print ("stdDev", best_stdDev)


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
      k = ldp.generic_OM(points[i][0], best_PhiO, points[i][1], best_LamO, best_Skew, test_kO)
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_kO_sample[i] = d # build the sample array of d values
   avg = np.average(try_kO_sample)
   if np.abs(best_avg) > np.abs(avg) :
       best_kO = test_kO
       best_avg = avg
   test_kO = test_kO + step_kO   # increase test value to next

print ("kO:",best_kO, "  average:", best_avg)


end = time.time()
total_mins2 = (end - start)/60

print("%.1f" % total_mins2, " mins total")
















