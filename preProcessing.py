import numpy as np
from subprocess import call 
from pyDOE import lhs
import os
from vtk import *

class preProcessing:
   def __init__(self, N):
       v = 0.5
       self.v= (1 - v**2)
       self.filename = 'bifurcation.in'
       self.f = open(self.filename,'w')
       self.beta  = 1.
       if os.path.exists("Runs/")==False:
           call(["mkdir", "Runs/"])
       self.N = N
##       self.length = length
##       self.constR = constR
##       self.constC = constC
##       self.area = area
       
   def calculateBeta(self, A1, A4, A5, A7):
       self.beta1  = self.calculatebeta(A1)
       self.beta4  = self.calculatebeta(A4)
       self.beta5  = self.calculatebeta(A5)
       self.beta7  = self.calculatebeta(A7)

   def getVals(self, length, areaFunc, wave, constC, constR, calcArea, betaFunc):
       self.length = length
       self.areaFunc = areaFunc
       self.constC = constC
       self.constR = constR
       self.wave = wave
       self.calcArea = calcArea
       self.betaFunc = betaFunc

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
       

   def writeInput_parallel(self, R, C, lamda, counter, nsteps=2.8e+06):
        try:
           self.filename = "bifurcation_%s.in"%counter
           self.f = open(self.filename,'w')
           self.f.write("12 parameter list\n") 
           self.f.write("0                       EQTYPE\n") 
           self.f.write("2.0E-4                  DT\n") 
           self.f.write("10.5E5                  NSTEPS\n") 
           # self.f.write("2.0E-4                 DT\n") 
           # self.f.write("2.28E4                  NSTEPS\n") 
           self.f.write("8E12                    IOSTEP\n") 
           self.f.write("200                    HISSTEP\n") 
           self.f.write("2                       INTTYPE\n") 
           self.f.write("1060                    Rho (Kg/m3)\n")
           self.f.write("0.75                   T (s)\n") 
           self.f.write("1.1                      Alpha (velocity profile param in the momentum eq) (for Alpha=1 no wall viscous effects considered)\n") 
           self.f.write("3.5E-03                  Viscosity (Pa*s)\n") 
           self.f.write("1.0                    Bscal\n") 
           self.f.write("7265.0                  pinf\n") 
           self.f.write("Mesh -- expansion order -- quadrature order Ndomains = 1\n") 
           self.f.write("1       nel domain 1 Beta Area Pext\n") 
           # self.f.write("0.0      0.04964    4  4 # x_lower x_upper L q\n") 
           self.f.write("0.0      %f    8  8 # x_lower x_upper L q\n"%self.length) 
           self.f.write("Beta = Bscal*(" + self.betaFunc + ")\n") 
           self.f.write("Ao = %f\n"%self.calcArea)
           self.f.write("Pext = 10000\n")
           self.f.write("Boundary conditions\n") 
           self.f.write("u  0     # A lhs boundary Domain 1\n") 
           self.f.write("   A = " + str(self.areaFunc) + "\n")
           self.f.write("q  0     # U lhs boundary\n")
           self.f.write("   q = " + str(self.wave) + "\n")
        
           if self.constC:
              self.f.write("W %e    # Compliance\n"%6.0E-9)
           else:
              self.f.write("W  %e    # Compliance\n"%C)

           if self.constR:
              self.f.write("W  %e    # Resistance\n"%0.55E+8)
           else:
              self.f.write("W  %e    # Resistance\n"%R)
          
           self.f.write("Initial conditions \n") 
           self.f.write("3 Lines Follow\n") 
           self.f.write("Given\n") 
           self.f.write("a = Ao\n") 
           self.f.write("u = 0.0\n") 
           self.f.write("History Pts \n") 
           self.f.write("1   #Number of Domains with history points\n") 
           self.f.write("1 1  #Npts Domain id x[0], x[1], x[2],...\n")
           self.f.write("%f\n"%self.length) 
           self.f.write("\n") 
           self.f.close()
        except Exception as e:
           print(e)
