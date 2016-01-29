"""
numerical schrodinger-poisson solver
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#layers format [thickness nm, energy eV]
layers = [[10, 2],[5,0],[10,2]]
dx=1
m0=9.109e-31

def run_schrodinger():
	"""
	Runs the schrodinger solver.
	It generates the structure/grid based on the layers variable
	
	"""
	struct = Structure(layers,dx)
	energies = calc_energy(struct,1)
	psis = calc_psi(struct,energies)
	x=struct.x
	
	plt.figure()
	plt.hold(True)
	plt.plot(x,struct.V,'k')
	for i in range(len(energies)):
		plt.plot(x,energies[i],'--')
		plt.plot(x,psis[i],'-')
	plt.xlabel("distance")
	plt.ylable("energy")
	
	
class Structure(object):
	def __init__(self,layers,dx):
		self.layers = layers
		self.dx = dx
		
		index = 0
		for layer in self.layers:
			layer_thickness = layer[0]
			layer_V = layer[1]
			layer_m = 1
			layer_numpoints = round2int(layer_thickness/self.dx)
			self.meff[index:index+layer_numpoints] = layer_m
			self.V[index:index+layer_numpoints] = layer_V
			index+=layer_numpoints
		self.numpoints = len(self.V)
		self.thickness = self.numpoints*self.dx
		self.x = np.linspace(0,self.thickness,self.dx)
	
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
	delta_e = 0.01
	e_max = 5
	E = [0.0]*numstates
	
	for i in range(0,numstates):	
		#first we find the frist energy where psi(inf) is zero
		#this corresponds with a bound state as psi is finite
		y2=calc_psi_inf(struct,e0)
		while e0<e_max:
			y1=y2
			e0+= delta_e
			y2=calc_psi_inf(struct,e0)
			if y1*y2<0:
				break
		finally:
			#error
			print "could't find a bound state"
		#linearly interpolate the last 2 points to get a better estimate
		e0 -= abs(y2)/(abs(y1)+abs(y2))*delta_e
		E[i] = e0
		#starting point for the next bound state
		e0+=delta_e
	return E
		
def psi_at_inf(struct,E):
	"""
	calculate the value of the wavefunction at the end of the domain
	This is used to find the bound state energies as this will
	tend to zero for a bound state and inf for non-bound states
	"""
	V=struct.V
	psi0=0.0
	psi1=1.0
	#I'm not exactly sure what this is doing...
	#I think it is just doing a dodgy solution of schrodinger and only taking the last point
	for i in range(1,struct.numpoints,1):
		psi2=((2*(dx/hbar)**2*(V[i]-E) + 2/meff[i])*psi1 - psi0/m[i])*m[i]
		psi0=psi1
		psi1=psi2
	return psi2
		
def calc_psi(struct,E):
	"""
	Given a structure and a list of bound state energies, calculate
	the wavefunctions for each one
	"""
	dx=struct.dx
	V=struct.V
	meff=struct.meff
	x=range(0,struct.numpoints)
	
	plt.figure()
	plt.hold(True)
	for i in range(0,len(E)):
		e0=E[i]
		
		def g(y,x):
			y0=y[0]
			y1=y[1]
			y2=2*meff[x]*(V[x]-e0)*y0/hbar**2
			return y1,y2
			
		#initial conditions for psi and psi' at x=0
		init = 0,0
		ans = odeint(g,init,x)
		psi[i] = ans[:,0]
		plt.plot(x,psi[:,0])
	plt.show()
	plt.xlabel("distance")
	plt.ylabel("psi")
	
	return psi
	
	
	