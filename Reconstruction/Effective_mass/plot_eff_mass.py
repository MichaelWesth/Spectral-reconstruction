import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from subprocess import *
import sys

outfile = 'results/Zero_temp.dat'   #Saves the data here
path = 'Gen2l_128x32.meson.g5.uu'   #Uses this data
test = 'therm'                      #Uses this kernel, can be therm, zero or avg
output=check_output('./eff '+path+' 'test'> '+outfile, shell=True) 

fin = open(outfile, 'r'); 

lines=fin.readlines(); 

n = 1

x=[]; y=[]; z=[];
#reads the output from c++ code, and stores it.
for l in lines:
    temp=l.split();
    x.append(float(temp[0]));
    y.append(float(temp[1]));
    z.append(float(temp[2]));


plt.grid()
#plt.semilogy(x, y)             #plot data on logplot
plt.errorbar(x, y, yerr=z)      #plot data with error
#plt.ylim(0, 0.1)   
plt.xlabel('$ t $')
plt.ylabel('$ E_0 $')
plt.show()

plt.show();
