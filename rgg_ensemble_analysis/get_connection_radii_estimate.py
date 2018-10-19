"""

Function to estimate the required connection radii
for RGGs. This can be used to estimate the connection
radius required to obtain RGGs with a given expected degree
for RGGs with arbitary position distributions. 

"""

#Import packages
import numpy as np
from scipy.spatial.distance import pdist
import periodic_distances as Cypdist
import scipy
import socket
import math
from scipy import special

def Get_Required_Kappa(Kappa_Desired , Network_Size, Position_Samples , d , Boundaries = 'S' , positions = None , pattern="uniform") :

	"""
	
	Function to determine the mean degree parameter required to obtain Kappa_Desired
	for an RGG with the specified conditions. 
	
	Given a connection radius, r, the mean degree parameter k_tilde is the volume of the ball in
	d dimensions multiplied by the number of nodes. This is the expected mean degree for
	a homogenous RGG in a region with volume 1. 
	
	This function supports generation of example point distributions for uniform in [0,1]^d with 
	periodic or solid boundary conditions. For Non-uniform point distributions we must have the
	postions inputted as an argument.
	
	Parameters
	--------------
	
	Kappa_Desired : float
	
	Choosen mean degree
	
	Network_Size : int
	
	Number of nodes in the graph of interest.
	
	Position_Samples : int
	
	Reduced number of samples we take to estimate the distribution of distances. Given M sample positions we will obtain
	M(M-1)/2 samples from the distribution of distances.
	
	d : int 
	
	dimension of the embedding space. 
	
	boundaries : str (optional)
	
	choice of boundary condition. For RGGs in [0,1]^d
	we can choose between 'S' (solid) and 'P' (periodic
	or toroidal). 
	
	positions : list
	
	list of of lists containing positions. Used if we
	have previously sampled the RGG positions. 
	
	pattern : str
	
	Point pattern to sample node positions from:
	"uniform" and "gaussian" are possible choices.
	
	Returns 
	----------
	
	d_threshold : float
	
	Estimated conneciton radius required
	to generate a network with the desired expected degree.
	
	
	"""

	#Make sure the required number of position samples is less than the network size:
	#An error will also occur when K >> N as the list coord will be greater than it should be.
	
	if Position_Samples > Network_Size :
		print("Error: position samples greater than network size")
	
	#If we don't get given a sample of positions then generate some:
	if positions is None and pattern == "uniform" :
		positions = np.random.uniform(0, 1.0, (Position_Samples, d))
	elif positions is None and pattern == "gaussian":
		positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=Position_Samples)

	
	#Compute the array of distances (periodic distance computation sped up using Cython).
	if Boundaries == 'S' : 
		Distances = pdist(positions)
	else :
		Distances = Cypdist.Periodic_Distances( positions , 1.0 )
		Distances = list(Distances)

	#Sort the distances:
	Sorted_Distances = np.sort(Distances)

	#Find the value required to get the edges above the required threshold:
	Edges_Required = (Kappa_Desired*Network_Size)/2.0

	#How far through the list to we have to be? F represents the 'fraction' through the list that we have got. 
	F = (2*Edges_Required)/(Network_Size*(Network_Size-1))
	List_Coord = int(0.5*F*(Position_Samples*(Position_Samples-1)))

	#Get the 'Edge coun'th' element (but check that it isn't greater than the list length)
	if List_Coord > len(Sorted_Distances) :
		d_threshold = Sorted_Distances[-1]
	else :
		d_threshold = Sorted_Distances[List_Coord]
	
	return d_threshold
	
