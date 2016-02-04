# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 13:21:55 2016

@author: matthewg
"""
from scipy.integrate import odeint, trapz
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def round2int(x):
    return int(x+0.5)
    
#e=1.6e-19
e=1
#m0=9.109e-31
hbar=1.05e-34
m0=hbar**2

class Structure(object):
    def __init__(self,layers,dx):
        self.layers = layers
        self.dx = dx
        self.thickness=sum(layer[0] for layer in layers)
        self.numpoints = round2int(self.thickness/self.dx)        
        self.meff=np.zeros(self.numpoints)
        self.V=np.zeros(self.numpoints)
        
        index = 0
        for layer in self.layers:
            layer_thickness = layer[0]
            layer_V = layer[1]*e
            layer_m = 1*m0
            layer_numpoints = round2int(layer_thickness/self.dx)
            self.meff[index:index+layer_numpoints] = layer_m
            self.V[index:index+layer_numpoints] = layer_V
            index+=layer_numpoints
        self.x = np.linspace(0,self.thickness,self.numpoints)
        
        self.V_interp=interp1d(self.x, self.V,bounds_error=False,fill_value=max(self.V))
        
        print "length x = {}".format(len(self.x))
        print "length V = {}".format(len(self.V))
        
        
    def plot_bandstructure(self):
        plt.figure()
        plt.plot(self.x,self.V)
 
def Wave_function(E, struct):
    """
    Calculates wave function psi for the given value
    of energy E and returns value at point b
    """
    
    def f(psi, x, params):
        """
        Returns derivatives for the 1D schrodinger eq.
        State0 is first derivative of the wave function psi, and state1 is its second derivative.
        """
        #print "x={}".format(x)
        E,V,m = params
        state0 = psi[1]
        state1 = 2.0*m/hbar**2*(V(x) - E)*psi[0]
        return np.array([state0, state1]) 
    
    psi0 = np.array([0,1])       # Wave function initial states
    params=[E, struct.V_interp, struct.meff[0]]
    psi = odeint(f, psi0, struct.x, args=(params,))
    return psi[:,0]
    
def Wave_function_inf(E, struct):
    """
    Wrapper function which just calculates the wavefunction at the end of the 
    simulation domain. This is checked against the 0 value. For bound states
    psi(inf)->0. For all other states it diverges
    """
    psi=Wave_function(E, struct)
    return psi[-1]
    
    
    
def find_bound_states(struct, do_plots=False):
    """
    finds the zeros of psi_inf vs. E.
    Bound states have psi_inf==0 while it diverges for all other energies
    THis funciton looks for a sign change in psi_inf and then calls 
      brentq over this small range to get an accurate value for the zero.
      This is the bound state energy
    """
    Es = np.linspace(min(struct.V), max(struct.V), 100)   # vector of energies where we look for the stable states
    
    psi_inf = []      # vector of wave function at x -> inf (ie, end of domain)
    for e1 in Es:
        psi_inf.append(Wave_function_inf(e1, struct))     # for each energy e1 find the the psi(x) at x = b    

    print psi_inf
    
    E_zeroes = []
    s = np.sign(psi_inf)
    for i in range(len(psi_inf)-1):
        if s[i]+s[i+1] == 0:
            zero = brentq(lambda energy: Wave_function_inf(energy, struct), Es[i], Es[i+1])
            E_zeroes.append(zero)
    
    
    if do_plots:
        plt.figure()
        plt.plot(Es,psi_inf)
        plt.title('Values of the $\Psi(\inf)$ vs. Energy')
        plt.xlabel('Energy, $E$')
        plt.ylabel('$\Psi(x -> inf)$', rotation='horizontal')
        for E in E_zeroes:
            plt.plot(E, [0], 'go')
            plt.annotate("E = %.2f"%E, xy = (E, 0), xytext=(E, 30))
            plt.grid()    
    
    return E_zeroes
    
    
def normalise(psi,x):
    """
    Normalises a wavefunction.
    Note that this assumes it is a bound state (ie. that the integral is finite)
    """
    integral=trapz(np.abs(psi),x)
    print "integral = {}".format(integral)
    return psi/integral
 
 
def main():
    # main program        
    #layers = [[20e-9,10],[10e-9,0],[20e-9,10]]
    #dx=0.1e-9
    layers = [[20,10],[1,0],[20,10]]
    dx=0.1
    struct = Structure(layers,dx)
    
    struct.plot_bandstructure()

    bound_states = find_bound_states(struct,True)
    
    print "Bound state energies:"
    for E in bound_states:
        print " E = {:.2g}".format(E)


 

 
    # Plot the wavefunctions for first 4 eigenstates
#    col=['r','g','b','c']
#    plt.figure(2)
#    plt.plot(x,[V(i) for i in x],'k',label='Potential')
#    for i,E in enumerate(E_zeroes[0:4]):
#        psi=Wave_function(E,V,x)
#        psi=normalise(psi,x)
#        plt.plot(x, E+psi, label="E = %.2f"%E, color=col[i])
#        plt.plot(x,E*np.ones(len(x)),'--',label=None,color=col[i])
#    plt.legend(loc="upper right")
#    plt.title('Wave function')
#    plt.xlabel('x, $x/L$')
#    plt.ylabel('$\Psi(x)$', rotation='horizontal', fontsize = 15)
#    plt.grid()
 
if __name__ == "__main__":
    main()