# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 22:45:00 2019

@author: Guangze Chen
"""

# This script draws eigen-wavefunction distribution on the lattice for 2D systems, and generates figs.3(c) and (d) in the article 

import kwant
import tinyarray as ta
from math import sin, cos, sqrt, pi

lat_u = kwant.lattice.square(a=1,name='u',norbs=1)
lat_d = kwant.lattice.square(a=1,name='d',norbs=1)

def make_system(W,theta,kz): #make a 2D cut of the original system with fixed kz
    syst = kwant.Builder()
    lat = kwant.lattice.square(a=1,name='u',norbs=2)
    sigma_x = ta.array([[0, 1], [1, 0]])
    sigma_y = ta.array([[0, -1j], [1j, 0]])
    sigma_z = ta.array([[1, 0], [0, -1]])

    x,y=W,W
    M=0.7 #same parameter as before
    k0=1
    
    for i in range(x):
        for j in range(y):
            syst[lat(i,j)] = M*(k0**2-4-kz**2)*sigma_z+(-sin(theta)/sqrt(2)*sin(kz)-1j*sin(theta)/sqrt(2)*sin(kz))*sigma_x
                        
            if i>0:
                syst[lat(i,j),lat(i-1,j)] = 1j/2*(cos(theta/2)**2)*sigma_x+1j/2*(sin(theta/2)**2)*sigma_y+M*sigma_z
                           
            if j>0:
                syst[lat(i,j),lat(i,j-1)] = 1j/2*(cos(theta/2)**2)*sigma_y+1j/2*(sin(theta/2)**2)*sigma_x+M*sigma_z
    
    syst1=syst.finalized()
    return syst1

#plot strength of wavefunction on each site
def plot_data(syst,n):
    import scipy.linalg as la
    ham = syst.hamiltonian_submatrix()
    evecs = la.eigh(ham)[1]
    psi=evecs[:, n] #the nth wavefunction of the system

    J_0 = kwant.operator.Density(syst) 
    density=J_0(psi)
    
    wf = abs(psi)**2

    kwant.plot(syst, site_size=density*(1/wf.max()), site_color=(0, 0, 1, 0.3),
               hop_lw=0.1) #plot the density over the maximum value of density on each site
    
def main():
    syst2=make_system(40,pi/2-pi/6,-0.4) #the system is 40*40 sites, which means there are 3200 bands
    for i in range(15,16):    
        plot_data(syst2,1600+i) #KWANT sort eigenvectors with eigenvalues increasing, due to particle-hole symmetry, the 1600th lies in the middle
    
    #fig3(c) and (d) are obtained with kz=-0.2, i=8 and kz=-0.4, i=15, respectively
	
if __name__ == '__main__':
    main()
