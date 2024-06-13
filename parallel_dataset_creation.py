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

#   Change the directory name to the folder that the analysis will take place

directory     = "/home/sokolor2/Geoms"
directoryRuns = "/home/sokolor2/Data"

# directory     = "/home1/g/gkissas/OperatorLearningVascular/BifurcationGeom"
# directoryRuns = "/scratch/gkissas/OperatorVascular/Runs5000Geom/sanityRuns"
    
def function_m( params):  
    
#   Change the directory name to the folder that the analysis will take place      
        counter, (R, C, l1, lamda) = params
        identity = "%s"%(counter)

        preProcessing.writeInput_parallel( R, C, l1, lamda, counter)
        
        # try:
        call(["./oneDBio.vserial",  "bifurcation_%s.in"%(counter)])
        files =  glob.glob("*bifurcation_%s.his"%counter)
        filesOut =  glob.glob("*bifurcation_%s.out"%counter)
        filesTex =  glob.glob("*bifurcation_%s.tex"%counter)
        filesTxt =  glob.glob("*bifurcation_%s.txt"%counter)
        filesIn =  glob.glob("*bifurcation_%s.in"%counter)
        postProcess = postProcessing(identity, filesOut, filesTex, filesTxt, filesIn, lamda=lamda)
        postProcess.readOutput_Parallel(R, C, l1, lamda, directory, directoryRuns, counter, files, plot=True, save=True)
        # except:
        #     print("Except")
        # utilities.clear_parallel(directory, counter)
        
def serial_debugging_test(param_values):
    counter = 0
    R2, R4, C2, C4,l1, l2, l3, l4, lamda1, lamda2, lamda3 = param_values[0,:]
    function_m((counter,(R2, R4, C2, C4, l1, l2, l3, l4, lamda1, lamda2, lamda3)))       
    print("All good")
    
if __name__ == "__main__": 
    utilities = utilities(directoryRuns)
    utilities.clear(directory)
    # utilities.purge()
    

    startTime = datetime.now()

    sobolSample     = True
    debuggingMode   = False
    datasetCreation = True
    test_sample     = False
    sanityCheck     = False

    if sobolSample:
        numberOfSamples = 16
        problem = {
            'num_vars': 4,
            'names'   : ['R', 'C','l1', 'lamda'],
            'bounds'  : [[8.3E+7,2.0E+8],
                        [1.5E-9,1.05E-8],
                        [0.3, 0.7],
                        [0.9,1.1]]
            }
        param_values = saltelli.sample(problem, numberOfSamples, calc_second_order=False)
        ite = param_values.shape[0]
    
    preProcessing = preProcessing(numberOfSamples)
    
    
    if debuggingMode:
        serial_debugging_test(param_values)

    if datasetCreation:
        with futures.ProcessPoolExecutor(max_workers=4) as executor:
            executor.map(function_m,[(counter,(R, C, l1, lamda)) \
                                  for counter,(R, C, l1, lamda) in enumerate(param_values)])    

    if test_sample:
        counter = 1000000
        # R2 = 6.0E+09
        # R4 = 4.0E+08
        # C2 = 6.0E-10
        # C4 = 30.0E-09
        N = 0.1

        R2 = 2.102840E+09
        R4 = 1.667377E+08
        C2 = 2.538206E-10
        C4 = 9.002950E-09
        l1 = 0.04964
        l2 = 0.0532
        l3 = 0.08866
        l4 = 0.03225
        lamda1 = 1.0
        lamda2 = 1.0
        lamda3 = 1.0

        function_m((counter,(R2, R4, C2, C4, l1, l2, l3, l4, lamda1, lamda2, lamda3)))       

    if sanityCheck:
        import numpy as np
        filename = "data_ol_test.npy"
        d = np.load(filename)

        R2 = d[:,0][:,None]*1e+10
        R4 = d[:,1][:,None]*1e+09
        C2 = d[:,2][:,None]*1e-09
        C4 = d[:,3][:,None]*1e-07
        l1 = d[:,4][:,None]
        l2 = d[:,5][:,None]
        l3 = d[:,6][:,None]
        l4 = d[:,7][:,None]
        lamda0 =  8.683831e-05
        lamda1 = np.ones_like(R2)
        lamda2 = lamda0*np.ones_like(R2)
        lamda3 = np.ones_like(R2)

        numberOfSamples = R2.shape[0]

        param_values = np.concatenate((R2,R4,C2,C4,l1,l2,l3,l4,lamda1,lamda2,lamda3),axis=-1)

        with futures.ProcessPoolExecutor(max_workers=20) as executor:
            executor.map(function_m,[(counter,(R2, R4, C2, C4, l1, l2, l3, l4, lamda1, lamda2, lamda3)) \
                                  for counter,(R2, R4, C2, C4, l1, l2, l3, l4, lamda1, lamda2, lamda3) in enumerate(param_values)])    
        
        # for i in range(R2.shape[0]):
        #     print("%e"%R2[i],"%e"%R4[i], "%e"%C2[i], "%e"%C4[i], l1[i], l2[i], l3[i], l4[i], lamda1[i], lamda2[i], lamda3[i])
        #     function_m([i,(R2[i], R4[i], C2[i], C4[i], l1[i], l2[i], l3[i], l4[i], lamda1[i], lamda2[i], lamda3[i])])

    # utilities.plotResultsVelocity(ite, lthreshold=1, identityVessel=5, diagnostics=False)
    # utilities.plotResultsFlow(ite, lthreshold=1, diagnostics=True)
