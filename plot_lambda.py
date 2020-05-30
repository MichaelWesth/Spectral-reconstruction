import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from subprocess import *
import sys

#cFile = 'python3 Plot_points.py '
cFile = './therm_kern'
lamb = '0.0001'
eps = '0.05'
E0 = '0.00001'
E = '0.06'
tmin = '1'
tmax = '64'
T = '32'
test = 'delta'
#path = 'Test2_2pt_d0'
path = 'Pion_data/Gen2l_16x32.meson.g5.uu'
#path = 'Gen2l_24x32.meson.g5.us'
#path = 'Gen2l_128x32.meson.g5.us'
#path = '128_ss_pion'
kern = 'therm'
min = '0'
max = '0.1'
num = '50'

pmin = sys.argv[1]
pmax = sys.argv[2]
pnum = sys.argv[3]

outfile = "data.dat"

pmin = float(pmin)
pmax = float(pmax)
pnum = int(pnum)

ls = np.linspace(pmin, pmax, pnum)
#ls = [16, 20, 24, 32, 40, 48, 56, 64, 128];
#ls = [0.01, 0.005, 0.001, 0.0005];
#ls_eps = [0.0673173, 0.0415415, 0.0392893, 0.0327828, 0.0247748, 0.018018, 0.016016, 0.0132633, 0.00650651] #zero_kern opt_vals [0.078097, 0.0538889, 0.045753, 0.0384309, 0.0289159, 0.0216779, 0.0198134, 0.0161877, 0.00780252]
ls_eps = [0.0950951, 0.0788288, 0.0643143, 0.0485485, 0.0395395, 0.0325325, 0.0282783, 0.0247748, 0.0122623] #therm_kern  opt_vals [0.114245, 0.101128, 0.0821657, 0.0605941, 0.0480301, 0.0405909, 0.0362493, 0.0318793, 0.0157264]

if (test == 'eps'):
    for i in range(pnum):
        tmax = str(int(ls[i])/2)
        T = str(int(ls[i]))
        path = 'Pion_data/Gen2l_'+T+'x32.meson.g5.uu'
        outfile = 'Result/pion_opteps_data'+T+'.dat'

        #output=check_output(cFile+' '+lamb+' '+eps+' '+E0+' '+E+' '+tmin+' '+tmax+' '+T+' '+test+' '+path+' '+kern+' '+min+' '+max+' '+num+' > '+outfile, shell=True) 
        
        fin = open(outfile, 'r'); 

        lines=fin.readlines(); 

        x=[]; y=[]; 

        for l in lines:
            temp=l.split();
            x.append(float(temp[0]));
            y.append(float(temp[1]));
        
        print(np.min(y))
        plt.plot(x[:400], y[:400], label = '$t_{max} $ = '+tmax )
        plt.legend(); 
     
    plt.plot([0, 0.1],[0, 0.1], label='$\epsilon = \epsilon_{opt}$', color='black') 
    plt.legend()
    plt.grid()
    plt.xlabel('$\epsilon $')  
    plt.ylabel('$\epsilon _{opt} $')
    plt.title('$E $ = '+E+'$, E_0 $ = '+E0)
    plt.title('$\lambda $ = '+lamb)
    #plt.ylim(-25,70)
    #plt.savefig('../Plots/'+test+'/e'+E+'.png')
   
elif (test == 'W'):
    for i in range(pnum):
        #lamb = str(np.power(10,ls[i]))
        tmax = str(int(ls[i]))
        outfile = 'Result/W_data'+tmax+'.dat'
        #output=check_output(cFile+' '+lamb+' '+eps+' '+E0+' '+E+' '+tmin+' '+tmax+' '+T+' '+test+' '+path+' '+kern+' '+min+' '+max+' '+num+' > '+outfile, shell=True) 
        
        fin = open(outfile, 'r'); 

        lines=fin.readlines(); 

        x=[]; y=[]; 

        for l in lines:
            temp=l.split();
            #x.append(float(temp[0]));
            #y.append(float(temp[1]));
            x.append(float(temp[2]));
            y.append(float(temp[3]));
        
        #plt.plot([i*int(tmax) for i in x],[i*int(tmax) for i in y], label = '$t_{max} $ = '+tmax )
        plt.plot(x,y, label = '$t_{max} $ = '+tmax )
        plt.legend();

    plt.grid()    
    plt.xlabel('$\omega \cdot t_{max}$')  
    plt.ylabel('$\sigma \cdot t_{max}$')    
    #plt.xlabel('$\omega  $')  
    #plt.ylabel('$ W(\omega, E) $')
    #plt.title('$t_{max} $ = '+tmax+'$, E_0 $ = '+E0)
    #plt.ylim(0,70)    
    plt.savefig('../plots/'+test+'_therm_eps'+eps+'tmax'+tmax+'.png')
    
elif (test == 'rho'):
    A = []
    for i in range(pnum):
        #lamb = '10e-'+str(int(ls[i]))
        #lamb = str(ls[i])
        tmax = str(int(ls[i]/2))
        T = str(int(ls[i]))
        #path = 'Pion_data/Gen2l_'+T+'x32.meson.g5.uu'
        #path = 'Data/Gen2l_'+T+'x32.meson.g5.us'
        #outfile = 'Result/pion_rho_data'+T+'.dat'
        path = 'Pion_data/Gen2l_'+T+'x32.meson.g5.uu'
        outfile = 'Data/pion_rho_eps'+eps+'_data'+T+'_withOptEpsTherm_Thalf.dat'
        eps = str(ls_eps[i])
        #output=check_output(cFile+' '+lamb+' '+eps+' '+E0+' '+E+' '+tmin+' '+tmax+' '+T+' '+test+' '+path+' '+kern+' '+min+' '+max+' '+num+' > '+outfile, shell=True) 
            
        fin = open(outfile, 'r'); 
        lines=fin.readlines(); 

        x=[]; y=[]; z=[];
        for l in lines:
            temp=l.split();
            x.append(float(temp[0]));
            y.append(float(temp[1]));
            z.append(float(temp[2]));

        #A.append(np.max(y)*ls[i])   
        plt.errorbar(x, y, z, fmt = '-', ms=2, label = '$T $ = ' + T)
        #plt.hist(y, bins=30)
        #plt.errorbar(x, y, z, label = '$\lambda $ = ' + lamb)
        #plt.errorbar(x, y, z, label = '$\epsilon $ = ' + eps)
        #plt.plot(x, y, label = '$\lambda $ = ' + lamb)
        plt.legend()
        #plt.show()

    #plt.plot(ls, A)
    plt.grid()
    plt.xlabel('$E $')  
    plt.ylabel('$ rho (E) $')
    plt.title('$t_{max} $ = T/2$, \lambda $ = '+lamb+'$, \epsilon = \epsilon_{opt}$')
    #plt.ylim(0,70)
    #plt.savefig('../pion_plots/'+test+'_therm_eps'+eps+'tmax'+tmax+'with__error.png')

elif (test == 'delta'):
    for i in range(pnum):
        tmax = str(int(ls[i]))
        output=check_output(cFile+' '+lamb+' '+eps+' '+E0+' '+E+' '+tmin+' '+tmax+' '+T+' '+test+' '+path+' '+kern+' '+min+' '+max+' '+num+' > '+outfile, shell=True) 
            
        fin = open(outfile, 'r'); 
        lines=fin.readlines(); 

        x=[]; y=[]; z=[];
        for l in lines:
            temp=l.split();
            x.append(float(temp[0]));
            y.append(float(temp[1]));
            z.append(float(temp[2]));
        
        plt.plot(x, y, label = '$tmax $ = ' + tmax)
        plt.legend()
        
        #plt.plot(x, y, label = '$\lambda $ = ' + lamb)
        #plt.legend()
        #plt.plot(x, z, label = '$\lambda $ = ' + lamb)
        #plt.legend()
        

    plt.xlabel('$E $')  
    plt.ylabel('$ delta (E)$')
    plt.title('$t_{max} $ = '+tmax+'$, E_0 $ = '+E0)
    #plt.ylim(0,70)
    plt.savefig('../plots/'+test+'_therm_eps'+eps+'tmax'+tmax+'.png')

plt.show()