"""
numerical schrodinger-poisson solver
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#layers format [thickness nm, energy eV]
layers = [[20, 5],[10,0],[20,5]]
dx=1
m0=9.109e-31
hbar=1.05e-34
e=1.6e-19



def run_schrodinger():
    """
    Runs the schrodinger solver.
    It generates the structure/grid based on the layers variable
    """
    struct = Structure(layers,dx)
    energies = calc_energy(struct,1)
    #energies=[0.5,1,1.5,2,2.5,3]
    #psis = calc_psi(struct,energies)
    x=struct.x

    plt.figure()    
    plt.hold(True)
    plt.plot(x,struct.V,'k')
    #for i in range(len(energies)):
    #    plt.plot(x,np.ones(len(x))*energies[i],'--')
    #	plt.plot(x,psis[i],'-')
    plt.xlabel("distance")
    plt.ylabel("energy")
	
	
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
        
        print "length x = {}".format(len(self.x))
        print "length V = {}".format(len(self.V))
	
def round2int(x):
    return int(x+0.5)
	
	
def calc_energy(struct, numstates):
    """
    Given a structure, calculate the bound state energies.
    It does this by finding the energies where the value of the 
    wavefunction at x=inf (or practically, at the edge of the 
    simulation domain) == 0. If it is not a bound state then the
    wavefunction will diverge at infinity.
    """
    #this is very simplistic at the moment
    #the aestimo code implements several zero finding routines
    #manually. Surely there is a scipy routine which does this better
    
    e0=0
    delta_e = 0.001
    e_max = 1
    E = [0.0]*numstates
    
    for i in range(0,numstates):	
        #first we find the frist energy where psi(inf) is zero
        #this corresponds with a bound state as psi is finite
        psi_inf=calc_psi(struct,[e0])
        psi_inf=psi_inf[0,-1]
        psis=[psi_inf]
        es=[e0]
        while e0<e_max:
            e0+= delta_e
            psi_inf=calc_psi(struct,[e0])
            psi_inf=psi_inf[0,-1]
            es.append(e0)
            psis.append(psi_inf)
            print "e={}, y2 = {}".format(e0,psi_inf)
            if psis[-1]*psis[-2]<0:
                break
        else:
            #error
            print "could't find a bound state"
        plt.figure()
        plt.plot(es,psis)
        plt.show()
        plt.xlabel("energy")
        plt.ylabel("psi_inf")
        #linearly interpolate the last 2 points to get a better estimate
        e0 -= abs(psis[-1])/(abs(psis[-2])+abs(psis[-1]))*delta_e
        E[i] = e0
        #starting point for the next bound state
        e0+=delta_e
    
    return E
		
def calc_psi_inf(struct,E):
    """
    calculate the value of the wavefunction at the end of the domain
    This is used to find the bound state energies as this will
    tend to zero for a bound state and inf for non-bound states
    """
    V=struct.V
    meff=struct.meff
    psi0=0.0
    psi1=1.0
    #I'm not exactly sure what this is doing...
    #I think it is just doing a dodgy solution of schrodinger and only taking the last point
    print "  psi2 = "    
    for i in range(1,struct.numpoints,1):
        psi2=((2*(dx/hbar)**2*(V[i]-E) + 2/meff[i])*psi1 - psi0/meff[i])*meff[i]
        psi0=psi1
        psi1=psi2
        print "{} ".format(psi2)
    return psi2
		
def calc_psi(struct,E):
    """
    Given a structure and a list of bound state energies, calculate
    the wavefunctions for each one
    """
    V=struct.V
    meff=struct.meff
    dx=struct.dx
    #meff=np.ones(struct.numpoints)
    #hbar=1
    
    psi=np.zeros((len(E),struct.numpoints))
    for i in range(0,len(E)):
        e0=E[i]*e
        psi[i,0]=0
        psi[i,1]=1
        for j in range(1,struct.numpoints-1):
            newval=2*(meff[j+1]/hbar**2*(V[j+1]-e0)*dx**2 + 1)*psi[i,j] - psi[i,j-1]
            print newval
            psi[i,j+1]=newval
    
    return psi
    
    
def test():
    struct = Structure(layers,dx)
    psi = calc_psi2(struct,0.01)
    return psi
    
    
def calc_psi2(struct, E):
    V=struct.V
    meff=struct.meff[0]
    x=range(struct.numpoints)
    print x
    
    def f(y, x, params):
        print "y = "
        print y
        print "x = "
        print x
        psi, phi = y
        E, m, V = params
        derivs = [phi, 2*m/hbar**2*(V[x]-E)*psi]
        return derivs
    
    #initial conditions
    y0=[0,1e-4]
    params = [E, meff, V]
    psoln=odeint(f, y0, x, args=(params,))
    return psoln
    
	
if __name__ == "__main__":
    run_schrodinger()
    
    
    
    
    
    
    
    
    
    