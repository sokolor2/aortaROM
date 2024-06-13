from utilities import utilities
import os
import glob
import numpy

##directoryRuns = "/home/sokolor2/Data"
##numberOfIterations = 120
##
##utilities = utilities(directoryRuns)
###utilities.plottingTest(directoryRuns, numberOfIterations)
##utilities.writeOutput(directoryRuns, 102, 0, 0)
###flow = utilities.createFlows(0.000264, 0.75)
###print(flow[0])

for i in range(143,197):
    os.rename("/mnt/e/constC/%s/RCR_%s.dat"%(i,i+1),"/mnt/e/constC/%s/RCR_%s.dat"%(i,i))
