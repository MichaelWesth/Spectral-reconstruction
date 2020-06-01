import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from subprocess import *
import sys

#cFile = 'python3 Plot_points.py '
cFile = './BG'
lamb = str(10**-6)
E0 = '0.00001'
E = '1' 
tmin = '1'
tmax = '90'
T = '320'
test = 'rho'   #defines the test being run
path = '../Pion_data/Gen2l_128x32.meson.g5.uu'
kern = 'therm'
min = E0
max = '0.25'
num = '100'

pmin = sys.argv[1]
pmax = sys.argv[2]
pnum = sys.argv[3]

outfile = "data1.dat"

pmin = float(pmin)
pmax = float(pmax)
pnum = int(pnum)

#-------------------------------parameters for LQCD data
ls = [16, 20, 24, 32, 40, 48, 56, 64, 128]; #tmax
#ls = [40, 48, 56, 64, 128]; #tmax
#ls = [64, 128]; #tmax



#------------------------------parameters for O(3)-model data
#ls = [30, 60, 90, 120, 150]; #tmax
#ls = [4,5,6,7,8]; #-log(lambda)


#---------------------------------Code for plotting 

if (test == 'rho'):
    A = []
    for i in range(pnum):
        #lamb = str(10**-ls[i])
        #lamb = str(ls[i])
        tmax = str(int(ls[i])/2)
        T = str(int(ls[i]))
        path = 'Pion_data/Gen2l_'+T+'x32.meson.g5.uu'
        outfile = 'Result/pion_Thalf_rho_data'+T+'.dat'
        #eps = str(ls[i])
        #output=check_output(cFile+' '+lamb+' '+E0+' '+E+' '+tmin+' '+tmax+' '+T+' '+test+' '+path+' '+kern+' '+min+' '+max+' '+num+' > '+outfile, shell=True) 
            
        fin = open(outfile, 'r'); 
        lines=fin.readlines(); 

        x=[]; y=[]; z=[];
        for l in lines:
            temp=l.split();
            x.append(float(temp[0]));
            y.append(float(temp[1]));
            z.append(float(temp[2]));

        #A.append(np.max(y)*ls[i])   
        plt.errorbar(x, y, z, fmt = '-', ms=2, label = '$tmax $ = ' + tmax)
        #plt.errorbar(x, y, z, fmt = '-', ms=2, label = '$\lambda $ = ' + lamb)
        plt.legend()
        #plt.show()

    #plt.plot(ls, A)
    plt.grid()
    plt.xlabel('$E $')  
    plt.ylabel('$ \hat{\\rho} (E) $')
    #plt.title('$\lambda $ = '+lamb)
    #plt.ylim(0,70)

elif (test == 'delta'):
    for i in range(pnum):
        #lamb = str(np.power(10,ls[i]))
        output=check_output(cFile+' '+lamb+' '+E0+' '+E+' '+tmin+' '+tmax+' '+T+' '+test+' '+path+' '+kern+' '+min+' '+max+' '+num+' > '+outfile, shell=True) 
            
        fin = open(outfile, 'r'); 
        lines=fin.readlines(); 

        x=[]; y=[]; z=[];
        for l in lines:
            temp=l.split();
            x.append(float(temp[0]));
            y.append(float(temp[1]));
            z.append(float(temp[2]));
        
        plt.plot(x, y, label = '$\lambda $ = ' + lamb)
        plt.legend()
        plt.plot(x, z, label = '$\lambda $ = ' + lamb)
        plt.legend()
        

    plt.xlabel('$E $')  
    plt.ylabel('$ delta (E)$')
    #plt.ylim(0,70)
    #plt.savefig('../Plots/'+test+'/eps'+eps+'.png')

plt.show()