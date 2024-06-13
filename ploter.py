# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
    
directory = './'
    
number_of_domains = 7
plt_values= {}
fig = plt.figure()
ax = fig.add_subplot(111)
color=iter(cm.rainbow(np.linspace(0,1,7)))

for filename in os.listdir(directory):
    if filename.endswith("cylinder.his"):
        print(filename)
        col_number = 17
        f = open(filename,"r")
        data = f.readlines()
        f.close()
        
        data_header = data[:3]
        data_no_header = data[4:]
        N_d = len(data_no_header)
        
        data_no_header = "".join(data_no_header)
        data_no_header = "".join(data_no_header.split("1\n"))
        data_no_header = list(map(float, data_no_header.split()))
        data_no_header = np.asarray(data_no_header)[:,None]
        
        values = np.zeros((N_d,col_number))
        
        for i in range(0,col_number):
            for j in range(0,N_d):
                values[j,i] = data_no_header[j*col_number + i, 0]
        
#        dictionary = {}
#        dictionary = {"Time": values[:,0],
#                      "Pressure": values[:,1],
#                      "Velocity": values[:,2],
#                      "Area": values[:,3],
#                      "Flow": values[:,4]}
        
#           np.save(filename,dictionary)


#        ax.cla() # or ax.clear()
        c=next(color)
        ax.plot(values[:,0], values[:,4]*1e6, 'ro',linewidth=1, markersize=0.5, label= "Flow_%s"%filename, c=c)
        ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('Flow in mL/s')
ax.set_title('Flow plot')
fig.savefig('flow.png')        

