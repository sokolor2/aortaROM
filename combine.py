from subprocess import call 
from postProcessing import postProcessing
from preProcessing import preProcessing
from utilities import utilities
from SALib.sample import saltelli
from concurrent import futures
from datetime import datetime
import glob
from scipy.stats import qmc
import numpy as np
import os
from numpy.fft import fft, ifft
import cmath
import shutil

directory     = "/home/sokolor2/Geoms"
directoryRuns = "/home/sokolor2/Data"
constC = True
constR = False

# EVERYTHING NEEDS TO BE IN SI UNITS (kg,m,s)
                    
def function_m( params):
    counter, (R, C, lamda) = params
    identity = "%s"%(counter)
    
    preProcessing.writeInput_parallel( R, C, lamda, counter)
    try:
        call(["./oneDBio.vserial",  "bifurcation_%s.in"%(counter)])
        files =  glob.glob("*bifurcation_%s.his"%counter)
        filesOut =  glob.glob("*bifurcation_%s.out"%counter)
        filesTex =  glob.glob("*bifurcation_%s.tex"%counter)
        filesTxt =  glob.glob("*bifurcation_%s.txt"%counter)
        filesIn =  glob.glob("*bifurcation_%s.in"%counter)
        postProcess = postProcessing(identity, filesOut, filesTex, filesTxt, filesIn, lamda=lamda)
        postProcess.readOutput_Parallel(R, C, lamda, directory, directoryRuns, counter, files, plot=True, save=True)
    except Exception as e:
        print(e)

if __name__ == "__main__":
    utilities = utilities(directoryRuns)
    utilities.clear(directory)
    length = []
    avgRadius = []

    lengthDat = open('length.dat','r') # aorta length data file (mm)
    for line in lengthDat:
        length.append(float(line))
    lengthDat.close()

    avgRadiusDat = open('avgRadius.dat','r') # avg radius data file (mm^2)
    for line in avgRadiusDat:
        x = line.split()
        avgRadius.append(float(x[0]))
    avgRadiusDat.close()


    # Generate sobo sample before iterating through vessels
    numberOfSamples = 4
    problem = {
        'num_vars': 3,
        'names'   : ['R', 'C'],
        'bounds'  : [[0.71E+8,2.98E+8],
                    [3.3E-9,26.8E-9],
                    [0.9,1.1]]
        }
    param_values = saltelli.sample(problem, numberOfSamples, calc_second_order=False)
    ite = param_values.shape[0]

    try:
        preProcessing = preProcessing(numberOfSamples)
    except Exception as e:
        print(e)

    for i in range(198,200):
        # iterate through each aorta in the list
        # add functionality to get length and radius values into this code itself
        # add functionality to generate sum of sines and cosines for new input waveforms
        # currently done seperately
        ################################################################################
        period = 0.75
        threshold = 0.0000025
        radiusData = open('radiusDat/geomData' + str(i) + '.dat')
        lst = []
        radsMM = []
        rads = []
        areaFunc = "3.14*("
        for line in radiusData:
            lst += [line.split()]

        for b in range(1,len(lst)-1):
            radsMM.append(float(lst[b][3]))

        
        radsMM = radsMM[10:len(radsMM)-50]
        radsMM.reverse()
        
        for c in range(0, len(radsMM)):
            rads.append(radsMM[c]/1000)
            
        xRange = np.arange(0,length[i],length[i]/len(rads))
        fit = np.polyfit(xRange,rads,20)
        for b in range(0,len(fit)):
            if b < len(fit)-1:
               areaFunc = areaFunc + str(fit[b]) + '*x^' + str(len(fit)-b-1) + ' + '
            else:
                areaFunc = areaFunc + str(fit[b]) + ')^2'
        calcArea = np.square(avgRadius[i]/1000)*np.pi

        betaFunc = utilities.betaFunction(fit, xRange)

        flowData = utilities.createFlows(period, threshold)
        waveData = open("wave.dat","w")
        for item in flowData:
            waveData.write(item)

        
        # Generate flow waveforms from FFT
                
        for m in range(0, len(flowData)):
            print(flowData[m])
            preProcessing.getVals(length[i]/1000, areaFunc, flowData[m], constC, constR, calcArea, betaFunc)
            datasetCreation = True
            if datasetCreation:
                with futures.ProcessPoolExecutor(max_workers=10) as executor:
                    executor.map(function_m,[(counter,(R, C, lamda)) \
                                          for counter,(R, C, lamda) in enumerate(param_values)])
            print("AORTA #%s, FLOW #%d"%(i,m))
            utilities.writeOutput(directoryRuns, 102 , i , m, constC, constR)
            label = os.listdir(directoryRuns)

            # Save and Clear .npz files for next wave
            os.chdir("/home")
            if constC:
                os.chdir("/mnt/e/constC/%s"%i)
            elif constR:
                os.chdir("/mnt/e/constC/")
            else:
                os.chdir("/mnt/e/Runs/")
            
            if not os.path.exists("Data_%s_%d"%(i,m)):
                os.makedirs("Data_%s_%d"%(i,m), mode=0o777)

            os.chdir("/home/sokolor2")
            
            for k in range(0,len(label)):
                if label[k].endswith(".npz"):
                    os.chmod(os.path.join(directoryRuns,label[k]),0o777)
                    origin = os.path.join(directoryRuns,label[k])
                    target = "/mnt/e/constC/%s/Data_%s_%d/"%(i,i,m) + label[k]
                    shutil.copyfile(origin,target)
                    os.remove(os.path.join(directoryRuns,label[k]))
            

        print("AORTA #%s COMPLETE"%i)
