# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 13:21:55 2016

@author: matthewg
"""
from scipy.integrate import odeint, trapz
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import numpy as np




def square_well(x, L, V0):
    if  abs(x) > L:
        return V0
    else:
        return 0

 
def Wave_function(E, V, x):
    """
    Calculates wave function psi for the given value
    of energy E and returns value at point b
    """
    
    def f(psi, x, params):
        """
        Returns derivatives for the 1D schrodinger eq.
        State0 is first derivative of the wave function psi, and state1 is its second derivative.
        """
        E,V = params
        state0 = psi[1]
        state1 = 2.0*(V(x) - E)*psi[0]
        return np.array([state0, state1]) 
    
    psi0 = np.array([0,1])       # Wave function initial states
    params=[E,V]
    psi = odeint(f, psi0, x, args=(params,))
    return psi[:,0]
    
def Wave_function_inf(E, V, x):
    """
    Wrapper function which just calculates the wavefunction at the end of the 
    simulation domain. This is checked against the 0 value. For bound states
    psi(inf)->0. For all other states it diverges
    """
    psi=Wave_function(E, V, x)
    return psi[-1]
    
    
    
def find_bound_states(E,psi_b,V,x):
    """
    finds the zeros of psi_inf vs. E.
    Bound states have psi_inf==0 while it diverges for all other energies
    THis funciton looks for a sign change in psi_inf and then calls 
      brentq over this small range to get an accurate value for the zero.
      This is the bound state energy
    """
    all_zeroes = []
    s = np.sign(psi_b)
    for i in range(len(psi_b)-1):
        if s[i]+s[i+1] == 0:
            zero = brentq(lambda energy: Wave_function_inf(energy,V,x), E[i], E[i+1])
            all_zeroes.append(zero)
    return all_zeroes
    
    
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
    V0 = 10
    L = 1
    N = 1000                  # number of points to take
    b = 2
    x = np.linspace(-b, b, N)    # x-axis
    

    en = np.linspace(0, V0, 100)   # vector of energies where we look for the stable states
    V= lambda x: square_well(x,L,V0)
    
    psi_b = []      # vector of wave function at x = b for all of the energies in en
    for e1 in en:
        psi_b.append(Wave_function_inf(e1, V, x))     # for each energy e1 find the the psi(x) at x = b
    E_zeroes = find_bound_states(en, psi_b, V, x)   # now find the energies where psi(b) = 0 
 
    # Print energies for the bound states
    print "Energies for the bound states are: "
    for E in E_zeroes:
        print "%.2f" %E
 
    # Plot wave function values at b vs energy vector
    plt.figure()
    plt.plot(en/V0,psi_b)
    plt.title('Values of the $\Psi(b)$ vs. Energy')
    plt.xlabel('Energy, $E/V_0$')
    plt.ylabel('$\Psi(x = b)$', rotation='horizontal')
    for E in E_zeroes:
        plt.plot(E/V0, [0], 'go')
        plt.annotate("E = %.2f"%E, xy = (E/V0, 0), xytext=(E/V0, 30))
    plt.grid()
 
    # Plot the wavefunctions for first 4 eigenstates
    col=['r','g','b','c']
    plt.figure(2)
    plt.plot(x,[V(i) for i in x],'k',label='Potential')
    for i,E in enumerate(E_zeroes[0:4]):
        psi=Wave_function(E,V,x)
        psi=normalise(psi,x)
        plt.plot(x, E+psi, label="E = %.2f"%E, color=col[i])
        plt.plot(x,E*np.ones(len(x)),'--',label=None,color=col[i])
    plt.legend(loc="upper right")
    plt.title('Wave function')
    plt.xlabel('x, $x/L$')
    plt.ylabel('$\Psi(x)$', rotation='horizontal', fontsize = 15)
    plt.grid()
 
if __name__ == "__main__":
    main()