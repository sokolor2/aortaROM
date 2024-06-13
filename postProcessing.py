     # -*- coding: utf-8 -*-
import numpy as np
from subprocess import call
import os,  os.path
import matplotlib as mlp
mlp.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from glob import glob

class postProcessing():
   def __init__(self, identity, filesOut, filesTex, filesTxt, filesIn, num_points = 1, lamda=1):
    
        self.num_points = num_points
        self.identity = identity
        self.col_number = 18
        self.color=iter(cm.rainbow(np.linspace(0,1,7)))
        self.T = 0.75
        self.convergence_index = []
        for filename in filesOut:
            call(["rm", "%s"%filename])

        call(["rm", "%s"%filesTex[0]])
        call(["rm", "%s"%filesTxt[0]])
        #call(["rm", "%s"%filesIn[0]])
        
   def readOutput_Parallel(self, R, C, lamda, directory, directoryRuns, counter, files, \
                  plot=True, save = True, keepPeriod = True, diagnostics = True):
        counterF = 0
        counterS = 0
        for filename in files:
            f = open(filename,"r")
            data = f.readlines()
            f.close()

            num_header_lines =self.num_points + 3
            num_no_header_lines = num_header_lines 
            data_no_header = data[num_no_header_lines:]
            N_d = int(len(data_no_header)/self.num_points)
            
            data_no_header = "".join(data_no_header)
            data_no_header = "".join(data_no_header.split(" \n"))
            data_no_header = list(map(float, data_no_header.split()))
            data_no_header = np.asarray(data_no_header)[:,None]

            temporary = data_no_header.reshape((N_d, self.num_points, self.col_number))
            self.values = np.moveaxis(temporary,0,1)  
            name = os.path.splitext(filename)[0]
            flag, counter = self.check_periodicity(self.T, name)
            print(flag)
            if flag == 1:
                if keepPeriod:
                    i_min, i_max = self.find_index(self.T, self.values, counter)
                    self.values = self.values[:,i_min:i_max,:]
                print("Pre-save")
                if save:
                    print("Saving")
                    self.save(R, C, lamda, directoryRuns, self.values, self.identity, name)       
                counterS += 1
                #call(["rm", "%s"%filename])
                print("Post-save")
            else:
                counterF += 1
                #call(["rm", "%s"%filename])
        # if glob.glob("%s_*.his"%counter) ==False :
        #     print("The code did not converge!")
        #     counterF += 1
            
        if diagnostics:
            f = self.run_diagnostics(self.identity, counterS, counterF)

   def define_filename(self, name, identity):
        vessel_id = os.path.splitext(name)[0].split("_")[3]   
        
        if identity < 10:
            counter = 0
        elif identity < 1000 and identity > 10:
            counter = -1
        elif identity <10000 and identity >= 1000:
            counter = -2
        elif identity <10000 and identity >= 1000:
            counter = -3
        else:
            print("Something is wrong with filename in parallel save")
        
        name2 = '_'.join(name.split('_'))[:-3 + counter]
        filename = "%s"%identity + name2 + "%s"%vessel_id
        
        return filename

       
   def save(self, R, C, lamda, directoryRuns, values, identity, name):
        vessel_id = os.path.splitext(name)[0].split("_")[1]
        filename = "%s"%identity + "bifurcation_"
        #filename = self.define_filename(name, identity)
        np.savez_compressed(directoryRuns +"/%s"%filename, R=R, C=C, quantities=values, lamda=lamda)
            
   def plotter(self,  values):
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        c=next(self.color)
        ax.plot(values[:,:,0], values[:,:,1], 'ro',linewidth=1, markersize=0.5, c=c)

       
   def find_index(self, T, t_r, i):
        T_min = i*T*np.ones((t_r.shape[0]))[:,None]
        T_max = (i+1)*T*np.ones((t_r.shape[0]))[:,None]
        i_min = self.find_nearest_index(t_r, T_min)
        i_max = self.find_nearest_index(t_r, T_max)
        return i_min, i_max
    
   def find_nearest_index(self, array, value):
        array = np.asarray(array).T
        idx = (np.abs(array - value)).argmin()
        return idx

   def find_nearest_value(self, array, value):
        array = np.asarray(array).T
        idx = (np.abs(array - value)).argmin()
        return array[idx]
    
   def check_periodicity(self, T, name):
        t_r = self.values[:,:,0]
        number_of_periods = int(t_r[0,-1]/T)
        counter = 0
        vessel_id = os.path.splitext(name)[0].split("_")[1]
        print("Iteration: %s Has %d Periods"%(vessel_id,number_of_periods))
        #vess = os.path.splitext(name)[0].split("_")[2]
        for i in range(number_of_periods -1,number_of_periods):
            i_min, i_max  = self.find_index(T, t_r, i)     
            discrepancy = np.abs(self.values[:,i_min,1]-self.values[:,i_max,1])
            if discrepancy < 500 and discrepancy !=0:
                #print("The solution of iteration %s, vessel %s is periodic after period %d"%(vessel_id, vess, i))
                print("The solution of iteration %s is periodic after period %d"%(vessel_id, i))
                counter = i
                self.convergence_index.append(i)
                return 1, counter
            else:
                #print("The solution of iteration %s, vessel %s is aperiodic"%(vessel_id,vess))
                print("The solution of iteration %s is aperiodic"%(vessel_id))
       
        counter = -1
        return 0, -1
      
   def update_dict(self, identity, AllArea, AllPressure, AllVelocity, AllFlow):
        AllPressure.update({identity: self.values[:,:,1].T})
        AllVelocity.update({identity: self.values[:,:,2].T})
        AllArea.update({identity: self.values[:,:,3].T})
        AllFlow.update({identity: self.values[:,:,4].T})
        return AllArea, AllPressure, AllVelocity, AllFlow
        
   def run_diagnostics(self, iteration, counterS, counterF):
       dis = counterS - counterF  
       if dis == counterS:
           print("The solution algorithm for iteration %s converged successfully for all vessels."%iteration)
           return True
       else:
           print("The number of successfully converged solutions for iteration %s are %d."%(iteration, counterS))
           print("The number of unsuccessfully converged solutions for iteration %s are %d."%(iteration, counterF))
           return False
          
