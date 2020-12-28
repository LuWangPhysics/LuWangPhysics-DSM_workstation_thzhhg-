#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 11:17:19 2020

@author: luwang1
"""
from __future__ import print_function
import matplotlib.pyplot as plt
from math import pi as pi
import numpy  as np

#remove the thread module in the sys and add pach
#in the later

import sys
if 'threading' in sys.modules:
    del sys.modules['threading']
sys.path.insert(0, '/Users/luwang1/Library/Python/3.7/lib/python/site-packages/')

import os
os.environ['MKL_NUM_THREADS'] = '1'
#get current working directory
os.listdir(os.getcwd())
#from gevent import monkey
#monkey.patch_all()
#monkey.patch_all(thread=False)



print('load and save from python')
plt.rcParams['text.usetex'] = False
#plt.rcParams['font.family'] = "/Users/luwang1/opt/anaconda3/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/helvetica"
plt.rcParams['font.size'] = 20
def plot1(a,b):
  
    omega_thz= np.loadtxt(a+"thz_omega.txt")
    spec_thz= np.loadtxt(b+"thz_spec.txt")
    my_fig=plt.figure()
    plt.plot(omega_thz/2/pi/1e12,spec_thz,linewidth=4)
    plt.xlabel('Frequency (THz)')
    plt.ylabel('THz spectrum')
    my_fig.savefig(b+'thz_spec.pdf', dpi=300,\
                orientation='portrait', bbox_inches='tight')
    plt.close(my_fig) 
    
    omega_ir1= np.loadtxt(a+"ir1_omega.txt")
    spec_ir1= np.loadtxt(b+"ir1_spec.txt")  
    my_fig=plt.figure()
    plt.plot(omega_ir1/2/pi/1e12,spec_ir1,linewidth=4)
    plt.xlabel('Frequency (THz)')
    plt.ylabel('IR$_1$ spectrum')
    my_fig.savefig(b+'ir1_spec.pdf', dpi=300,\
                orientation='portrait', bbox_inches='tight')
    plt.close(my_fig) 
        
        
    omega_ir2= np.loadtxt(a+"ir2_omega.txt")
    spec_ir2= np.loadtxt(b+"ir2_spec.txt")
    my_fig=plt.figure()
    plt.plot(omega_ir2/2/pi/1e12,spec_ir2,linewidth=4)
    plt.xlabel('Frequency (THz)')
    plt.ylabel('IR$_2$ spectrum')
 
    my_fig.savefig(b+'ir2_spec.pdf', dpi=300,\
                orientation='portrait', bbox_inches='tight')
    plt.close(my_fig) 
    
    
    z = np.loadtxt(a+'z.txt')
    eff= np.loadtxt(b+'eff.txt')
    my_fig=plt.figure()
    plt.plot(z*1e9,eff,linewidth=4)
    plt.xlabel('z (nm)')
    plt.ylabel('eff (%)')
    my_fig.savefig(b+'eff.pdf', dpi=300,\
                 orientation='portrait', bbox_inches='tight')
    plt.close(my_fig) 
    
    
    ## the thz eff 3hg plot
    # eff_3hg= np.loadtxt(b+'eff_3hg.txt')
    # my_fig=plt.figure()
    # plt.plot(z*1e9,eff_3hg,linewidth=4)
    # plt.xlabel('z (nm)')
    # plt.ylabel('eff_3hg (%)')
    # my_fig.savefig(b+'eff_3hg.pdf', dpi=300,\
    #              orientation='portrait', bbox_inches='tight')
    # plt.close(my_fig) 
    
    
    
    plt.close('all')
    return


def plot2(a,b):
    z = np.loadtxt(a+'z.txt')

        
    econs= np.loadtxt(a+'energy_cons.txt')
    my_fig=plt.figure()
    plt.plot(z*1e9,econs,linewidth=4)
    plt.ticklabel_format(useOffset=False)
    plt.xlabel('z (nm)')
    plt.ylabel('energy cons')

    my_fig.savefig(a+'energy_cons.pdf', dpi=300,\
                 orientation='portrait', bbox_inches='tight')
    plt.close(my_fig)   
    
    econs= np.loadtxt(a+'energy_pump.txt')
    my_fig=plt.figure()
    plt.plot(z*1e9,econs,linewidth=4)
    plt.ticklabel_format(useOffset=False)
    plt.xlabel('z (nm)')
    plt.ylabel('energy pump')
 
    my_fig.savefig(a+'energy_pump.pdf', dpi=300,\
                 orientation='portrait', bbox_inches='tight')
    plt.close(my_fig)
        
    E_thz=np.genfromtxt(b+'et_thz.txt',delimiter=',')
    t=np.loadtxt(a+'t.txt')
    dt=t[2]-t[1];
    my_fig=plt.figure()
    if np.size(np.shape(E_thz))==1:
            max_p = np.argmax(abs(E_thz),axis=0)
            if E_thz[max_p]<0:
                E_thz=-1*E_thz
            plt.plot(t*1e12,np.fft.fftshift(E_thz[:]),linewidth=4)
    else:
         for k in range(np.size(np.shape(b))):
                 plt.plot(t*1e15,np.fft.fftshift(E_thz[:,k]),linewidth=4)
                 plt.legend(('re', 'im'))
    plt.xlabel('t (ps)')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel('E$_{THz}$ (V/m)')
    plt.xlim(-5, 5)
    my_fig.savefig(b+'et_thz.pdf', dpi=300,\
                 orientation='portrait', bbox_inches='tight')
    plt.close(my_fig) 
    
    plt.close('all')
    return




def plot_general(save_path,a_name,b_name,c_name):
    a = np.loadtxt(save_path+a_name)
    b= np.genfromtxt(save_path+b_name,delimiter=',')
    my_fig=plt.figure()
    if np.size(np.shape(b))==1:
         plt.plot(a,(b[:]))
    else:
         for k in range(np.size(np.shape(b))):
             plt.plot(a,(b[:,k]))
             plt.legend(('re', 'im'))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    my_fig.savefig(save_path+c_name+'.pdf', dpi=300,\
                 orientation='portrait', bbox_inches='tight')  
    plt.close(my_fig)  
    plt.close('all')
    return 




