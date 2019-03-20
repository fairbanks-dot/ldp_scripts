import numpy as np
import akdot_nr_ldp_v2 as ldp
import time
start = time.time()
# This is file for finding parameters for projectionsto best fit test points.
#
# Troy Hicks
# last revision: mar 18, 2019
#
# edits: improved display of time to compute, figured out good default values
# notes: pulling points from a file
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
# test points, in decimal degrees, Lat, long, then Ellipsoid Height in meters.
points = np.genfromtxt("PCCS_1E_weighted.csv", delimiter=',')                                                   



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
# Optimize the PhiO, LamO, and skew parameters
#
#      accomplished by finding the phiO, lamO, and skew that has the minimum
#      standard deviation of linear distortion per test point.
#
################################################################################
# search_Range = number of degrees to look from design points list center.
# risk is if it is too small you might miss a good solution. 1.5 to 3.5 suggested.
# step_Minutes = number of minutes to iterate. 6 minutes suggested.
# step_Skew = number of degrees to iterate. 6 degrees suggested.
# inital iterations = search_Range*2*60/step_Minutes*2*60/step_Skew
# 2.5*2*2*60/6*2*180/6 = 6000  for first loop
# 2.5*2*2*60/3*2*180/2 = 36000 for first loop

search_Range = 2.5                   # inital range box, in degrees around average value of design points
step_Minutes = 15                     # initial number minutes to step each time, later will iterate per 1 minute
step_Skew = 8                        # initial number of degrees to iterate per skew, later will iterate per 1 degree
range_PhiO = search_Range            # set initial boxed in range for phiO from center
range_LamO = search_Range            # set initial boxed in range for lamO from center
avg_PhiO = (np.average(points[:,0])) # initially centering the range on average phiO
avg_LamO = (np.average(points[:,1])) # initially centering the range on average lamO
center_PhiO = avg_PhiO
center_LamO = avg_LamO
step_PhiO = step_Minutes/60          # set inital step as decimal of degree (minutes converted to degree)
step_LamO = step_Minutes/60          # set inital step as decimal of degree (minutes converted to degree)
best_LamO = avg_LamO   # reset best lamO
best_PhiO = avg_PhiO   # reset best phiO
best_Skew = -90        # reset best skew
best_stdDev = 999999   # start with high number so can look for lower value
phase1 = False         # phase1 is not complete. During phase1 we use a very coarse search, and stop once it seems no longer gets better.

total_Iterations = search_Range*2*search_Range*2*60/step_Minutes*60/step_Minutes*180/step_Skew*points.shape[0]  # I think this is total number of first phase iterations

for c in range (3) :   # execute three phases if possible, if not then it looks all the way through and never gets far enough from best_stdDev   
   test_PhiO = center_PhiO - range_PhiO              # reset test_PhiO
   while test_PhiO <= center_PhiO + range_PhiO :     
      test_LamO = center_LamO - range_LamO #reset test_LamO
      while test_LamO <= center_LamO + range_LamO:  
         test_Skew = -89              # reset test_Skew
         while test_Skew <= 89 :            
            distortion_sample = np.zeros([points.shape[0]])  # reset the distortion sample for points
            for i in range(points.shape[0]) :
               Rg = geMeRaCu(points[i][2], a, inFl)
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
               print("******** WINNER! ************ PhiO","%.7f" % best_PhiO, "LamO", "%.7f" %best_LamO, "Skew", "%.1f" %best_Skew, "stdDev", "%.2f" % best_stdDev)
               print(" ")
            test_Skew = test_Skew + step_Skew
         print(".", end="")
         test_LamO = test_LamO + step_LamO
      print(" ")
      current = time.time()
      time_so_far = np.int((current-start)*100)
      time_so_far = time_so_far/100
      current_mins = (time_so_far)/60
      phase = c + 1
      print(" ")
      print ("%.1f" % current_mins, "minutes so far.                        phase = ", phase)
      print("phiO/lamO range: +/-", range_PhiO, "  phiO/lamO step: ", "%.4f" % step_PhiO, "degrees", "  skew step:", step_Skew)
      test_PhiO = test_PhiO + step_PhiO
   # finished a phase, no now lets focus in tighter
   center_PhiO = best_PhiO               # re-centering the range on current best phiO
   center_LamO = best_LamO               # re-centering the range on current best lamO
   range_LamO = range_LamO / step_Minutes*3         # 
   range_PhiO = range_PhiO / step_Minutes*3         # 
   step_LamO = step_LamO / (60 / step_Minutes)          # 
   step_PhiO = step_PhiO / (60 / step_Minutes)          # 
   step_Skew = step_Skew / 2                                  # 
print ("PhiO:","%.8f" % best_PhiO, "LamO","%.8f" % best_LamO, "Skew:", best_Skew)
print ("stdDev","%.2f" % best_stdDev)


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
tries_kO = np.int((max_kO-min_kO) / step_kO +0.5)
test_kO = min_kO  # start the test with the minimum value
best_kO = max_kO  # reset best 
best_avg = 99999999 # start with high number so can look for lower value
for c in range(tries_kO) :
   # for each overall kO try, sample the distortion for the test points
   try_kO_sample = np.zeros([points.shape[0]]) 
   for i in range(points.shape[0]) :
      Rg = geMeRaCu(points[i][2], a, inFl)
      EF = Rg/(Rg + points[i][2])
      k = ldp.generic_OM(points[i][0], best_PhiO, points[i][1], best_LamO, best_Skew, test_kO)
      d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm
      try_kO_sample[i] = d # build the sample array of d values
   avg = np.average(try_kO_sample)
   if np.abs(best_avg) > np.abs(avg) :
       best_kO = test_kO
       best_avg = avg
   test_kO = test_kO + step_kO   # increase test value to next

print ("kO:","%.7f" % best_kO, "  average:","%.2f" % best_avg)


# print out some stats on computation time
end = time.time()
total_mins2 = (end - start)/60
num_points = points.shape[0]
seconds_per_point = (end - start) / num_points
print(" ")
print(num_points, " design points used.")
print("%.1f" % total_mins2, " mins total")
print("%.1f" % seconds_per_point, " seconds per design point.")















