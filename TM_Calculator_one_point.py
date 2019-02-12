import numpy as np
import akdot_nr_ldp_v1 as ldp




# global constant for converting meters to us survey feet
usft = 3937.0 / 1200.0  # 1 meter = 3937/1200 US Survey feet

# test points, in decimal degrees, Lat, long, then Ellipsoid Height in meters.
points = np.genfromtxt("one point.csv", delimiter=',')                                                   



# Elevation Factor is based on ellipsoid to ground, nothing to do with any particular projection.
Rg = ldp.geMeRaCu(points[0,2])
EF = Rg/(Rg + points[0,2])

k = ldp.generic_TM(points[0,0], points[0,1], points[0,1], 1/EF)

d = (EF * k - 1) * 1000000  # d = linear distortion in units of ppm  USE THIS if calling generic function
   
start = ldp.generic_TM1(points[0,0], points[0,1], points[0,1], 1/EF)
end = ldp.generic_TM1(64.83701077, -147.5656682, points[0,1], 1/EF)
point1 = [start[0], start[1], 148.8012192]
point2 = [end[0], end[1], 148.8012192]

print ("Using the provided position and the selected projection, your results are ...")
print ("-------------------------------------------------------------------------------------------------------")
print (  "Grid scale factor (k):", k, "    (EF):", EF, "    CSF:", "%.12f" % (k* EF), "   d(ppm) =", "%.2f" % d)

print ("-------------------------------------------------------------------------------------------------------")
print ("-------------------------------------------------------------------------------------------------------")
print (point1, point2)

inv = (( (point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 ) **0.5)

print ("%.4f" % inv)

d2 = (EF * end[2] - 1) * 1000000
print (d2)




