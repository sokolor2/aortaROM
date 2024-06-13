import numpy as np
from vtk import vtkUnstructuredGridReader, \
                vtkDataSetMapper, vtkActor, \
                vtkRenderer, vtkRenderWindow,\
                vtkRenderWindowInteractor
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import glob
import os
from os import path
from subprocess import call 
from SALib.analyze import sobol, delta
import plotting
import re
import math
from numpy.fft import fft, ifft
import cmath

class utilities:
   def __init__(self, directoryRuns):
       self.directoryRuns = directoryRuns
       self.files = []

   def visualizeGeometry(self, file_name):
        
        # Read the source file.
        reader = vtkUnstructuredGridReader()
        reader.SetFileName(file_name)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update() # Needed because of GetScalarRange
        output = reader.GetOutput()
        scalar_range = output.GetScalarRange()
        
        # Create the mapper that corresponds the objects of the vtk file
        # into graphics elements
        mapper = vtkDataSetMapper()
        mapper.SetInputData(output)
        mapper.SetScalarRange(scalar_range)
        
        # Create the Actor
        actor = vtkActor()
        actor.SetMapper(mapper)
        
        # Create the Renderer
        renderer = vtkRenderer()
        renderer.AddActor(actor)
        renderer.SetBackground(1, 1, 1) # Set background to white
        
        # Create the RendererWindow
        renderer_window = vtkRenderWindow()
        renderer_window.AddRenderer(renderer)
        
        # Create the RendererWindowInteractor and display the vtk_file
        interactor = vtkRenderWindowInteractor()
        interactor.SetRenderWindow(renderer_window)
        interactor.Initialize()
        interactor.Start()
                
   def plotResultsPressure(self, numberOfIterations, identityVessel=4, \
                         uthreshold=160, lthreshold=60, diagnostics=False):
       
       # GOAL: Create a loop over the number of iterations for the chosen vessel and plot
       # At each iteration step you need to:
       # 1) Load the file with the iteration id and read numpy array stored in .npz file
       # 2) You need to have a control (if statement) because the stored pressure waveforms might have length equyal to zero
       # 3) Check if number of point in the steady state cardiac cycle solution is equal to a predifend number i.e. P =39 -> pressure=array["Pressure"][:P]
       # chose p equal to the least one in your dataset
       # 4) You need to plot those arrays
       fig1 = plt.figure(1,figsize=(15, 10), dpi=300, facecolor='w', frameon = False)
       fig1.clf()
       ax1 = fig1.add_subplot(111)  


       filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
                            %(0,2)
       data_t = np.load(filename, allow_pickle=True).item()            
       n = data_t["Pressure"][-1,:].shape[0]
       for i in range(1,numberOfIterations):
           if glob.glob(self.directoryRuns+"%s_bifurcation_%d.npy"\
           %(i,identityVessel)):
               filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
               %(i,identityVessel)
               data = np.load(filename, allow_pickle=True).item()
               if  np.shape(data["Pressure"]/133.3)==(1,0):
                   continue
               else:
                   p = data["Pressure"][-1,:n-1]/133.3
#                   print(p.shape)
                   t = data["Time"][-1,:n-1]
                   A1 = data["A1"]
                   C1 = data["C1"]
                   R1 = data["R1"]
                   
                   A4 = data["A4"]
                   C4 = data["C4"]
                   R4 = data["R4"]  
                   
                   A5 = data["A5"]
                   C5 = data["C5"]
                   R5 = data["R5"]
                                     
                   A7 = data["A7"]
                   C7 = data["C7"]
                   R7 = data["R7"]
                                    
                   
                   t = t - t.min()
                   ax1.plot(t,p)
    
                   try:
                    if p.max() < uthreshold and p.max() > lthreshold:
                        if diagnostics == True:
                               A1 = data["A1"]
                               C1 = data["C1"]
                               R1 = data["R1"]
                               
                               A4 = data["A4"]
                               C4 = data["C4"]
                               R4 = data["R4"]  
                               
                               A5 = data["A5"]
                               C5 = data["C5"]
                               R5 = data["R5"]
                                                 
                               A7 = data["A7"]
                               C7 = data["C7"]
                               R7 = data["R7"]
                               print(A1, A4, A5, A7, R1, R4, R5, R7,C1, C4,C5, C7)
                        t = t - t.min()
                        ax1.plot(t,p)
                   except ValueError:
                    pass
           else:
               continue
       ax1.set_xlabel('Time in $s$')
       ax1.set_ylabel('Pressure in $mmHg$')
       plt.savefig("pressure_for_threshold_%s.png"%lthreshold)

   def plotResultsVelocity(self, numberOfIterations, identityVessel=7, \
                         uthreshold=160, lthreshold= 60, diagnostics=False):
       
       fig2 = plt.figure(2,figsize=(15, 10), dpi=300, facecolor='w', frameon = False)
       fig2.clf()
       ax2 = fig2.add_subplot(111)  
       filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
                            %(0,2)       
       n = 48
       for i in range(1,numberOfIterations):
           if glob.glob(self.directoryRuns+"%s_bifurcation_%d.npy"\
           %(i,identityVessel)):
               filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
               %(i,identityVessel)
               data = np.load(filename, allow_pickle=True).item()
               p = data["Pressure"][-1,:n]/133.3
               u = data["Velocity"][-1,:n]
               t = data["Time"][-1,:n]
               A1 = data["A1"]
               C1 = data["C1"]
               R1 = data["R1"]
                
               A4 = data["A4"]
               C4 = data["C4"]
               R4 = data["R4"]  
                
               A5 = data["A5"]
               C5 = data["C5"]
               R5 = data["R5"]
                                  
               A7 = data["A7"]
               C7 = data["C7"]
               R7 = data["R7"]
               # try:
               if len(p)!=0 and len(t)!=0:                  
                   if p.max() < uthreshold and p.max() > lthreshold:
                        if diagnostics == True:
                                   A1 = data["A1"]
                                   C1 = data["C1"]
                                   R1 = data["R1"]
                                   
                                   A4 = data["A4"]
                                   C4 = data["C4"]
                                   R4 = data["R4"]  
                                   
                                   A5 = data["A5"]
                                   C5 = data["C5"]
                                   R5 = data["R5"]
                                                     
                                   A7 = data["A7"]
                                   C7 = data["C7"]
                                   R7 = data["R7"]
                                   print(A1, A4, A5, A7, R1, R4, R5, R7,C1, C4,C5, C7)
                        t = t - t.min()
                        # print(u.shape)
                        ax2.plot(t,u)
               # except ValueError:
               #  pass
               # print(u.shape)

           else:
               continue
       ax2.set_xlabel('Time in $s$',fontsize=20)
       ax2.set_ylabel('Velocity in $m/s$',fontsize=20)
       plt.savefig("velocity_for_threshold_%s.png"%lthreshold)


   def plotResultsFlow(self, numberOfIterations, identityVessel=7, \
                         uthreshold=160, lthreshold=60, diagnostics=False):
       
       fig3 = plt.figure(3,figsize=(15, 10), dpi=300, facecolor='w', frameon = False)
       fig3.clf()
       ax3 = fig3.add_subplot(111)  
       filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
                            %(0,2)
       data_t = np.load(filename, allow_pickle=True).item()            
       n = data_t["Pressure"][-1,:].shape[0]
       for i in range(1,numberOfIterations):
           if glob.glob(self.directoryRuns+"%s_bifurcation_%d.npy"\
           %(i,identityVessel)):
               filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
               %(i,identityVessel)
               data = np.load(filename, allow_pickle=True).item()
               p = data["Pressure"][-1,:n-1]/133.3
               q = data["Flow"][-1,:n-1]
               t = data["Time"][-1,:n-1]
               A1 = data["A1"]
               C1 = data["C1"]
               R1 = data["R1"]
                
               A4 = data["A4"]
               C4 = data["C4"]
               R4 = data["R4"]  
                
               A5 = data["A5"]
               C5 = data["C5"]
               R5 = data["R5"]
                                  
               A7 = data["A7"]
               C7 = data["C7"]
               R7 = data["R7"]
               try:
                if p.max() < uthreshold and p.max() > lthreshold:
                    if diagnostics == True:
                               A1 = data["A1"]
                               C1 = data["C1"]
                               R1 = data["R1"]
                               
                               A4 = data["A4"]
                               C4 = data["C4"]
                               R4 = data["R4"]  
                               
                               A5 = data["A5"]
                               C5 = data["C5"]
                               R5 = data["R5"]
                                                 
                               A7 = data["A7"]
                               C7 = data["C7"]
                               R7 = data["R7"]
                               print(A1, A4, A5, A7, R1, R4, R5, R7,C1, C4,C5, C7)
                    t = t - t.min()
                    ax3.plot(t,q)
               except ValueError:
                pass
           else:
               continue
       ax3.set_xlabel('Time in $s$',fontsize=25)
       ax3.set_ylabel('Flow in $L/s$',fontsize=25)
       plt.savefig("flow_for_threshold_%s.png"%lthreshold)

               
   def purge(self):
       for filename in glob.glob(os.path.join(self.directoryRuns, '*')):
                call(["rm",filename])
    
               
   def clear(self, directory):
        for counter in range(1,8):
            for filename in os.listdir(directory):
                if filename.endswith("%s.his"%counter) or filename.endswith("%s.out"%counter):
                   call(["rm", filename])
                if filename.endswith("%s.tex"%counter) or filename.endswith("%s.txt"%counter):
                   call(["rm", filename])
                if filename.endswith("_%s.in"%counter):  
                   call(["rm", filename])
                   
   def perform_SASobol(self, problem , ite, vesselId, \
                       uthreshold=190, lthreshold=10, QoI="Velocity"):
       
        # num_vars = problem["num_vars"]
        # samples = 2*num_vars + 2
        self.S1A1= []
        self.S1R1 = []
        self.S1C1 = []

        self.S1A4= []
        self.S1R4 = []
        self.S1C4 = []

        self.S1A5= []
        self.S1R5 = []
        self.S1C5 = []

        self.S1A7= []
        self.S1R7 = []
        self.S1C7 = []        
        
        self.STA1= []
        self.STR1 = []
        self.STC1 = []

        self.STA4= []
        self.STR4 = []
        self.STC4 = []

        self.STA5= []
        self.STR5 = []
        self.STC5 = []

        self.STA7= []
        self.STR7 = []
        self.STC7 = []   
        
        self.time = []
        
        self.S2betaR = []
        self.S2betaC = []

        self.S2Rbeta = []
        self.S2RC = []

        self.S2Cbeta = []
        self.S2CR = []
        
        self.S1_confE = []
        self.S2_confE = []
        
        filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
                             %(12,vesselId)
        data = np.load(filename, allow_pickle=True).item()

        n = data[QoI][-1,:].shape[0]
        self.find_nonFailed(ite,vesselId)
        successfulFiles = self.find_nonFailedNumber(12) 
        print(successfulFiles)
        for j in range(0,n-1):      
            Y= []
            for filename in self.files:
                    data = np.load(filename, allow_pickle=True).item()
                    if QoI == "Pressure":
                        p = data[QoI][-1,j]/133.3
                    else:
                        p = data[QoI][-1,j]
                    Y.append(p)
            self.YY= Y    
            Y = np.asarray(Y)
            Si = sobol.analyze(problem,Y, calc_second_order=False, \
                               print_to_console=False, parallel=True, n_processors=8) 
            self.Si = Si
            
            self.S1A1.append(Si["S1"][0])
            self.S1A4.append(Si["S1"][1])
            self.S1A5.append(Si["S1"][2])
            self.S1A7.append(Si["S1"][3])

            self.S1R1.append(Si["S1"][4])
            self.S1R4.append(Si["S1"][5])
            self.S1R5.append(Si["S1"][6])
            self.S1R7.append(Si["S1"][7])
            
            self.S1C1.append(Si["S1"][8])
            self.S1C4.append(Si["S1"][9])
            self.S1C5.append(Si["S1"][10])
            self.S1C7.append(Si["S1"][11])
        

            self.STA1.append(Si["ST"][0])
            self.STA4.append(Si["ST"][1])
            self.STA5.append(Si["ST"][2])
            self.STA7.append(Si["ST"][3])
            
            self.STR1.append(Si["ST"][4])
            self.STR4.append(Si["ST"][5])
            self.STR5.append(Si["ST"][6])
            self.STR7.append(Si["ST"][7])
            
            self.STC1.append(Si["ST"][8])
            self.STC4.append(Si["ST"][9])
            self.STC5.append(Si["ST"][10])
            self.STC7.append(Si["ST"][11])
                                  
            
            self.S1_confE.append(Si["S1_conf"][2])
            
            calc_second_order=False
            if calc_second_order:
                self.S2betaR.append(Si["S2"][0,1])
                self.S2betaC.append(Si["S2"][0,2])

                self.S2Rbeta.append(Si["S2"][1,0])
                self.S2RC.append(Si["S2"][1,2])

                self.S2Cbeta.append(Si["S2"][2,0])
                self.S2CR.append(Si["S2"][2,1])

                self.S2_confE.append(Si["S2_conf"][2])

        data = np.load(self.files[0], allow_pickle=True).item()
        tmin = data["Time"].min()
        t = data["Time"] - tmin
        self.time.append(t) 

   def perform_SADelta(self, problem , ite, vesselId = 4):
        self.S1beta = []
        self.S1R = []
        self.S1C = []
        
        self.time = []
        
        filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
                             %(0,vesselId)
        data = np.load(filename, allow_pickle=True).item()
        for j in range(data["Pressure"][-1,:].shape[0]):      
            Y= []
            X = []
            for i in range(0,ite):
                
                filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
                             %(i+1,vesselId)
                data = np.load(filename, allow_pickle=True).item()
                p = data["Pressure"][-1,j]/133.3
                beta = data["A0"]
                C = data["C0"]
                R = data["R0"]
                Y.append(p)
                X.append((beta,R,C))
            Y = np.asarray(Y)
            X = np.asarray(X)
            Si = delta.analyze(problem, X, Y, \
                               print_to_console=False)
            self.S1beta.append(Si["S1"][0])
            self.S1R.append(Si["S1"][1])
            self.S1C.append(Si["S1"][2])
            
        tmin = data["Time"].min()
        t = data["Time"] - tmin
        self.time.append(t) 
        
   def plot_SA(self,vesselId):
       
       fig4 = plt.figure(4,figsize=(15, 10))
       fig4.clf()
       ax4 = fig4.add_subplot(111)  
       
       self.S1A1 = np.asarray(self.S1A1)
       self.S1R1 = np.asarray(self.S1R1)
       self.S1C1 = np.asarray(self.S1C1)

       self.S1A4 = np.asarray(self.S1A4)
       self.S1R4 = np.asarray(self.S1R4)
       self.S1C4 = np.asarray(self.S1C4)

       self.S1A5 = np.asarray(self.S1A5)
       self.S1R5 = np.asarray(self.S1R5)
       self.S1C5 = np.asarray(self.S1C5)

       self.S1A7 = np.asarray(self.S1A7)
       self.S1R7 = np.asarray(self.S1R7)
       self.S1C7 = np.asarray(self.S1C7)       
       
       self.time = np.asarray(self.time)
       
       m1 = self.time.flatten().shape[0]
       m2 = self.S1A1.flatten().shape[0]
       
       if m1 > m2:
           n = m2
       else:
           n = m1
      
       ax4.plot(self.time.flatten()[:n],  self.S1A1.flatten()[:n], 'ro-',linewidth=2.5, markersize=6.5, label='A1')
       ax4.plot(self.time.flatten()[:n],  self.S1R1.flatten()[:n], 'mo-',linewidth=2.5, markersize=6.5, label='R1')
       ax4.plot(self.time.flatten()[:n],  self.S1C1.flatten()[:n], 'ko-',linewidth=2.5, markersize=6.5, label='C1')

       ax4.plot(self.time.flatten()[:n],  self.S1A4.flatten()[:n], 'r-',linewidth=2.5, markersize=6.5, label='A4')
       ax4.plot(self.time.flatten()[:n],  self.S1R4.flatten()[:n], 'm-',linewidth=2.5, markersize=6.5, label='R4')
       ax4.plot(self.time.flatten()[:n],  self.S1C4.flatten()[:n], 'k-',linewidth=2.5, markersize=6.5, label='C4')

       ax4.plot(self.time.flatten()[:n],  self.S1A5.flatten()[:n], 'r>-',linewidth=2.5, markersize=6.5, label='A5')
       ax4.plot(self.time.flatten()[:n],  self.S1R5.flatten()[:n], 'm>-',linewidth=2.5, markersize=6.5, label='R5')
       ax4.plot(self.time.flatten()[:n],  self.S1C5.flatten()[:n], 'k>-',linewidth=2.5, markersize=6.5, label='C5')

       ax4.plot(self.time.flatten()[:n],  self.S1A7.flatten()[:n], 'r--',linewidth=2.5, markersize=6.5, label='A7')
       ax4.plot(self.time.flatten()[:n],  self.S1R7.flatten()[:n], 'm--',linewidth=2.5, markersize=6.5, label='R7')
       ax4.plot(self.time.flatten()[:n],  self.S1C7.flatten()[:n], 'k--',linewidth=2.5, markersize=6.5, label='C7')              
       
       ax4.legend(loc=1, prop={'size': 22})  
       ax4.set_xlabel("Time in s")
       ax4.set_ylabel("First order Sobol indices")
       plt.savefig("First_order_Sobol_%d.png"%vesselId)
       
       calc_second_order=False
       if calc_second_order:
           fig5 = plt.figure(5,figsize=(15, 10))
           fig5.clf()
           ax5 = fig5.add_subplot(111)  
           

           self.S2betaR = np.asarray(self.S2betaR)
           self.S2betaC = np.asarray(self.S2betaC)
    
           self.S2Rbeta = np.asarray(self.S2Rbeta)
           self.S2RC = np.asarray(self.S2RC)
           
           self.S2Cbeta = np.asarray(self.S2Cbeta)
           self.S2CR = np.asarray(self.S2CR)
          
           ax5.plot(self.time.flatten()[:n],  self.S2betaR.flatten()[:n],linewidth=1, markersize=0.5, label='A-h')
           ax5.plot(self.time.flatten()[:n],  self.S2betaC.flatten()[:n],linewidth=1, markersize=0.5, label='A-E')

           ax5.plot(self.time.flatten()[:n],  self.S2Rbeta.flatten()[:n], linewidth=1, markersize=0.5, label='h-A')
           ax5.plot(self.time.flatten()[:n],  self.S2RC.flatten()[:n], linewidth=1, markersize=0.5, label='h-E')

           ax5.plot(self.time.flatten()[:n],  self.S2Cbeta.flatten()[:n], linewidth=1, markersize=0.5, label='E-A')
           ax5.plot(self.time.flatten()[:n],  self.S2CR.flatten()[:n], linewidth=1, markersize=0.5, label='E-h')
           
           ax5.legend(loc=1, prop={'size': 19})
           ax5.set_xlabel("Time in s")
           ax5.set_ylabel("Second order Sobol indices")
           
           plt.savefig("Second_order_Sobol_%d.png"%vesselId)
           dictionary = {}
           
           dictionary  = { "Time"  : self.time.flatten(), 
                        "S1beta"   :  self.S1beta.flatten(),
                        "S1R4"     : self.S1R4.flatten(),
                        "S1R7"     : self.S1R7.flatten(),
                        "S2betaR"  : self.S2betaR.flatten(),
                        "S2betaC"  : self.S2betaC.flatten(),
                        "S2Rbeta"  : self.S2Rbeta.flatten(),
                        "S2RC"     : self.S2RC.flatten(),
                        "S2Cbeta"  : self.S2Cbeta.flatten(),
                        "S2CR"     : self.S2CR.flatten()}
           
           np.save("Sobol_indices",dictionary)
       
   def find_nonFailed(self, totalIterations, vesselId):
       counter = 53
       for i in range(0, totalIterations):
           if glob.glob(self.directoryRuns +"%s_*_%s.npy"%(i,vesselId)):
               filename = self.directoryRuns+"%s_bifurcation_%d.npy"\
                                 %(i,vesselId)
               data = np.load(filename, allow_pickle=True).item()
               p = data["Pressure"][-1,:]

               if len(p)>0:
                   self.files.append(filename)
                   continue
               else:
                  counter = counter + 1
           else:
               counter = counter + 1
               print(" Number of failed attempts %d"%counter)
   
   def find_nonFailedNumber(self, params) : 
    # Find the quotient 
    m = params + 2
    n = len(self.files)
    q = np.floor(n / m) 
      
    # 1st possible closest number 
    n1 = m * q 
    # 2nd possible closest number 
    if((n * m) > 0) : 
        n2 = (m * (q + 1))  
    else : 
        n2 = (m * (q - 1)) 
      
    # if true, then n1 is the required closest number 
    if n1<n2 : 
        return int(n1 )
      
    # else n2 is the required closest number  
    return int(n2 )
  
   def plot_SA_total(self,vesselId):
           
           fig4 = plt.figure(4,figsize=(15, 10))
           fig4.clf()
           ax4 = fig4.add_subplot(111)  
           
           self.STA1 = np.asarray(self.STA1)
           self.STR1 = np.asarray(self.STR1)
           self.STC1 = np.asarray(self.STC1)
           
           self.STA4 = np.asarray(self.STA4)
           self.STR4 = np.asarray(self.STR4)
           self.STC4 = np.asarray(self.STC4)
    
           self.STA5 = np.asarray(self.STA5)
           self.STR5 = np.asarray(self.STR5)
           self.STC5 = np.asarray(self.STC5)
    
           self.STA7 = np.asarray(self.STA7)
           self.STR7 = np.asarray(self.STR7)
           self.STC7 = np.asarray(self.STC7)       
           
           self.time = np.asarray(self.time)
           
           m1 = self.time.flatten().shape[0]
           m2 = self.S1A1.flatten().shape[0]
           
           if m1 > m2:
               n = m2
           else:
               n = m1
          
           ax4.plot(self.time.flatten()[:n],  self.STA1.flatten()[:n], 'ro-',linewidth=1, markersize=4.5, label='A1')
           ax4.plot(self.time.flatten()[:n],  self.STR1.flatten()[:n], 'mo-',linewidth=1, markersize=4.5, label='R1')
           ax4.plot(self.time.flatten()[:n],  self.STC1.flatten()[:n], 'ko-',linewidth=1, markersize=4.5, label='C1')
                                                    
           ax4.plot(self.time.flatten()[:n],  self.STA4.flatten()[:n], 'r-',linewidth=1, markersize=4.5, label='A4')
           ax4.plot(self.time.flatten()[:n],  self.STR4.flatten()[:n], 'm-',linewidth=1, markersize=4.5, label='R4')
           ax4.plot(self.time.flatten()[:n],  self.STC4.flatten()[:n], 'k-',linewidth=1, markersize=4.5, label='C4')
                                                    
           ax4.plot(self.time.flatten()[:n],  self.STA5.flatten()[:n], 'r*-',linewidth=1, markersize=4.5, label='A5')
           ax4.plot(self.time.flatten()[:n],  self.STR5.flatten()[:n], 'm*-',linewidth=1, markersize=4.5, label='R5')
           ax4.plot(self.time.flatten()[:n],  self.STC5.flatten()[:n], 'k*-',linewidth=1, markersize=4.5, label='C5')
                                                    
           ax4.plot(self.time.flatten()[:n],  self.STA7.flatten()[:n], 'r--',linewidth=1, markersize=4.5, label='A7')
           ax4.plot(self.time.flatten()[:n],  self.STR7.flatten()[:n], 'm--',linewidth=1, markersize=4.5, label='R7')
           ax4.plot(self.time.flatten()[:n],  self.STC7.flatten()[:n], 'k--',linewidth=1, markersize=4.5, label='C7')              
           
           ax4.legend(loc=1, prop={'size': 22})  
           ax4.set_xlabel("Time in s")
           ax4.set_ylabel("Total Sobol indices")
           plt.savefig("Total_order_Sobol_%d.png"%vesselId)

   def writeOutput(self, directoryRuns, numberOfIterations, identity, sub, constC, constR):
      resistance = []
      compliance = []
      compFinal = []
      pres = []
      resFinal = []
      presFinal = []
      vel = []
      velFinal = []
      for i in range (0,numberOfIterations):
         exists = path.exists(directoryRuns + "/" + "%sbifurcation_.npz"%(i))
         #print(exists)
         if exists:
            filename = self.directoryRuns + "/" + "%sbifurcation_.npz"\
               %(i)
            data = np.load(filename, allow_pickle=True)
            p = data['quantities'][:,:,1]/133.3
            p = p.tolist()
            v = data['quantities'][:,:,2]
            v = v.tolist()
            velocity = v[0]
            pressure = p[0]
            #print(p[0])
            current = open("bifurcation_%s.in"%(i), "r")
            for line in current:
               if 'Resistance' in line:
                  resistance.append(float(line[2:15]))
               if 'Compliance' in line:
                  compliance.append(float(line[2:15]))
            current.close()
            pres.append(pressure)
            vel.append(velocity)
      if constC:
         pair = zip(resistance, pres, vel)
         sorted_pairs = sorted(pair)
         tuples = zip(*sorted_pairs)
         resistance, pres, vel = [list(tuple) for tuple in tuples]
      if constR:
         pair = zip(compliance, pres, vel)
         sorted_pairs = sorted(pair)
         tuples = zip(*sorted_pairs)
         compliance, pres, vel = [list(tuple) for tuple in tuples]
      seen = set()
      if constC:
         for item in resistance:
            if item not in seen:
               resFinal.append(item)
               presFinal.append(pres[resistance.index(item)])
               compFinal.append(compliance[resistance.index(item)])
               velFinal.append(vel[resistance.index(item)])
               seen.add(item)
      elif constR:
         for item in compliance:
            if item not in seen:
               compFinal.append(item)
               presFinal.append(pres[compliance.index(item)])
               resFinal.append(resistance[compliance.index(item)])
               velFina.append(vel[compliance.index(item)])
               seen.add(item)
      else:
         for item in resistance:
            resFinal.append(item)
            presFinal.append(pres[resistance.index(item)])
            compFinal.append(compliance[resistance.index(item)])
            velFinal.append(vel[resistance.index(item)])
            

      os.chdir("/home/sokolor2")
      if constC:
         os.chdir("/mnt/e/constC/")
      elif constR:
         os.chdir("/mnt/e/constR/")
      else:
         os.chdir("/mnt/e/Runs/")

      if not os.path.exists("%s"%(identity)):
         os.makedirs("%s"%(identity))

      os.chdir("%s"%(identity))   

      f = open("output_%s_%d.dat"%(identity,sub),"w")
      flow = open("flow_%s_%d.dat"%(identity,sub),"w")
      RCR = open("RCR_%s.dat"%(identity),"w")
      
      for item in presFinal:
         f.write(str(item))
         f.write("\n")
      for item in velFinal:
         flow.write(str(item))
         flow.write("\n")
      if constC:
         for item in resFinal:
            RCR.write(str(item) + " " + str(compFinal[resFinal.index(item)]))
      if constR:
         for item in compFinal:
            RCR.write(str(item) + " " + str(resFinal[compFinal.index(item)]))
      print(presFinal)
      f.close()
      RCR.close()
      os.chdir("/home/sokolor2")
            
   def createFlows(self, period, threshold):
      minFlow = 0.0002*0.9
      maxFlow = 0.0002*1.1
      flowRange = np.linspace(minFlow,maxFlow,5)
      tSys = (2/5)*period
      tSpan = np.linspace(0,0.75,750)
      allWave = []
      flowData = []
      for item in flowRange:
         func = []
         for i in range(0,len(tSpan)):
            if tSpan[i] < tSys:
               val = item*np.sin(np.pi*tSpan[i]/tSys)
               func.append(val)
            else:
               func.append(0)
         flowData.append(func)
      for i in range(0,len(flowData)):
          fft3 = np.fft.fft(flowData[i])
          x = np.arange(0, 0.75, 0.001)
          freqs = np.fft.fftfreq(len(x), 0.001)
          wave = ""
          recomb = np.zeros((len(x),))
          middle = len(x)//2 + 1
          for j in range(middle):
              if abs(fft3[j]) / len(x) > threshold:
                  if j == 0:
                      coeff = 2
                      #print(str(fft3[i].real/len(x)))
                      wave += str(fft3[j].real/len(x))
                  else:
                      coeff = 1
                      #print(str(fft3[i].real/(len(x)/2)) + "*cos(" + str(freqs[i]*2*np.pi) + "*t) - " + str(fft3[i].imag/(len(x)/2)) + "*sin(" + str(freqs[i]*2*np.pi) + "*t)")
                      wave += (" + " + (str(fft3[j].real/(len(x)/2)) + "*cos(" + str(freqs[j]*2*np.pi) + "*t) - "))
                      wave += (str(fft3[j].imag/(len(x)/2)) + "*sin(" + str(freqs[j]*2*np.pi) + "*t)")             
                  sinusoid = 1/(len(x)*coeff/2)*(abs(fft3[j])*np.cos(freqs[i]*2*np.pi*x+cmath.phase(fft3[j])))
                  recomb += sinusoid
                  plt.plot(x, sinusoid)

          allWave.append(wave)
      return(allWave)
               
                  
         
      return(flowData)

   def decompose_fft(data: list, threshold: float = 0.0):
    fft3 = np.fft.fft(data)
    x = np.arange(0, 0.75, 0.001)
    freqs = np.fft.fftfreq(len(x), 0.001)
    wave = ""
    recomb = np.zeros((len(x),))
    middle = len(x)//2 + 1
    for i in range(middle):
        if abs(fft3[i]) / len(x) > 0.0:
            if i == 0:
                coeff = 2
                #print(str(fft3[i].real/len(x)))
                wave += str(fft3[i].real/len(x))
            else:
                coeff = 1
                #print(str(fft3[i].real/(len(x)/2)) + "*cos(" + str(freqs[i]*2*np.pi) + "*t) - " + str(fft3[i].imag/(len(x)/2)) + "*sin(" + str(freqs[i]*2*np.pi) + "*t)")
                wave += (" + " + (str(fft3[i].real/(len(x)/2)) + "*cos(" + str(freqs[i]*2*np.pi) + "*t) - "))
                wave += (str(fft3[i].imag/(len(x)/2)) + "*sin(" + str(freqs[i]*2*np.pi) + "*t)")             
            sinusoid = 1/(len(x)*coeff/2)*(abs(fft3[i])*np.cos(freqs[i]*2*np.pi*x+cmath.phase(fft3[i])))
            recomb += sinusoid
            plt.plot(x, sinusoid)
    return wave

   def betaFunction(self, fit, xRange):
      betaFunc = ""
      poly = np.polyval(fit,xRange)
      areas = [np.pi * x**2 for x in poly]
      betaData = [self.calculatebeta(calcArea) for calcArea in areas]
      betaFit = np.polyfit(xRange, betaData, 20)
      for i in range(0,len(betaFit)):
         if i < len(fit)-1:
            betaFunc = betaFunc + str(betaFit[i]) + '*x^' + str(len(betaFit)-i-1) + ' + '
         else:
            betaFunc = betaFunc + str(betaFit[i])
      return betaFunc
            

   def calculatebeta(self, calcArea):
        try:
            a = np.sqrt(calcArea/np.pi)
            k1 = 2.0E+06
            k2 = -2253.00
            k3 = 86500.00
            v= 0.5
            numerator = (k1*np.exp(k2*a) + k3)
            denominator = (1 - v**2)*a*np.sqrt(np.pi) 
        except Exception as e:
            print(e)
        return numerator/denominator
   
   def plottingTest(self, directoryRuns, numberOfIterations, identityVessel=4, \
                         uthreshold=160, lthreshold=60, diagnostics=False):
      plt.axis('On')
      for i in range (0,numberOfIterations):
         exists = path.exists(directoryRuns + "/" + "%sbifurcation_.npz"%(i))
         print(str(exists))
         if exists:
            filename = self.directoryRuns + "/" + "%sbifurcation_.npz"\
               %(i)
            data = np.load(filename, allow_pickle=True)
            fig1 = plt.figure(1,figsize=(15, 15), dpi=300, facecolor='w', frameon = False)
            fig1.clf()
            fig1.set_figheight(5.5)
      
            ax1 = fig1.add_subplot(111)
            ax1.axes.clear()
            if  np.shape(data['quantities'][:,:,1]/133.3)==(1,0):
               #print(np.shape(data['quantities'][:,:,1]/133.3))
               continue
            else:
               p = data['quantities'][:,:,1]/133.3
               v = data['quantities'][:,:,2]
               v = v.tolist()
               p = p.tolist()
               pressure = p[0]
               #print(p[0])
               velocity = v[0]
               t = data['quantities'][:,:,0]
               t = t - t.min()
               t = t.tolist()
               time = t[0]
               ax1.plot(time,pressure)
               #ax1.set_ylabel('Pressure in $mmHg$')
               plt.title("Pressure Waveform (Model)")
               plt.xlabel("Time (s)")
               plt.ylabel("Pressure (mmHG)")
               plt.savefig("%s_pressure.png"%(i))
               #print(time)
               fig1.clf()
   

            
            

