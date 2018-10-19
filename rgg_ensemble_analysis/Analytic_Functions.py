"""

Contains functions which evaluate the analytic approximation of the 
algebraic connectivity of graphs.


"""

#Import standard python modules:
import numpy as np
import itertools
import math
import scipy
from scipy import special



def Theory_Algebraic(N,Kappa,d) :
	
	"""
	Theoretical approximation for the algebraic connectivity
	of periodic RGGs

	Parameters
	-----------

	N : int

	Network size

	Kappa :  float

	Expected mean degree of the ensemble

	d: int

	Dimension of the embedding space

	Returns
	---------

	Mu2 : float

	Approximate value of the algebraic connectivity. 

	"""

	# Calculate the radius from the epxcted mean degree:
	r = (1.0 / ((np.pi) ** 0.5)) * ((((Kappa) / N) * scipy.special.gamma((d + 2.0) / 2.0)) ** (1.0 / d))

	#Compute the algebraic connectivity:
	Mu2 = Kappa- N*(r**(d/2.0))*scipy.special.jv( (d/2.0) , 2*math.pi*r )

	return Mu2

def Lattice_Theory_Algebraic(N,Kappa,d) : 

	"""
	
	Approximate expression for the algebraic connectivity of a periodic 
	random geometric graph.
	
	This formula is presented in chapter 4 of:
	Nyberg, Amy. The Laplacian Spectra of Random Geometric Graphs. Diss. 2014.
	(Also see Van Mieghem, Piet. Graph spectra for complex networks. Cambridge University Press, 2010.)
	
	
	Parameters
	-----------
	
	N : int
	
	Network size
	
	Kappa :  float
	
	Mean degree of the ensemble
	
	d: int
	
	Dimension of the embedding space
	
	Returns
	---------
	
	Mu2 : float
	
	Approximate value of the algebraic connectivity. 
	
	"""
	
	Mu2 = (1.0/6.0)*( (math.pi/(N**(1.0/d)) )**2 )*( ( Kappa + 1.0)**( (d + 2.0)/float(d) ) )
	return Mu2
	


def ER_Theory(N,Kappa) :
	"""
	Apprixmate expression for the mean algerbaic connectivity of
	ER graphs. 
	
	This formula is presented in:
	
	Jamakovic, A., and Piet Van Mieghem. "On the robustness of complex networks by using the algebraic connectivity." NETWORKING 2008 Ad Hoc and Sensor Networks, Wireless Networks, Next Generation Internet (2008): 183-194.
	
	
	
	Parameters
	-----------
	
	N : int
	
	Network size
	
	Kappa :  float
	
	Mean degree of the ensemble
	
	
	Returns
	----------
	
	Mu2 : float
	
	Expected algebraic connectivity of the ensemble.
	
	"""
	Mu2 = Kappa - ( 2*Kappa*(1.0 - (Kappa/N))*math.log(N)  )**0.5 + (( (Kappa*(1.0 - (Kappa/N)))/math.log(N) )**0.5)*(  math.log( (2*math.pi*math.log((N**2)/(2*math.pi)))  ) - 0.5772)
	return Mu2	


	
