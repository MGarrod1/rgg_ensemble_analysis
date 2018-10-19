"""

Periodic distances.

Function to compute the distances between
points in [0,1]^d with periodic (toriodal)
boundaries. Cythonized to speed things up.

"""

from libcpp.vector cimport vector as cpp_vector

def Periodic_Distances( positions , Dom_Size = 1.0 ) :

	"""
	Computes the distance betweeen two points in the 
	domain [0,Dom_Size)^d with periodic boundaries. 
	
	Arguments
	-----------
	
	positions : list of lists
	
	list of positions of nodes.
	
	Dom size : float 
	
	The size of the square domain which the positions are embedded in.
	
	
	Returns
	----------
	
	Distances : list
	
	Array of pairwise distances between points in the domain.
	
	"""
	
	cdef int i
	cdef int j
	cdef int q
	cdef int N


	Distances = [ ] 
	q = len( positions[0][:] ) 
	N = len(positions) 
	
	for i in range(N) : 
		for j in range(i+1,N) : 

			dij = 0 
			#Loop over number of dimensions
			for k in range(q):
				# Compute the absolute distance
				dist = abs( positions[i, k] - positions[j, k] )
				if dist>0.5*Dom_Size : 
					dist = Dom_Size - dist
				# Add to the total distance
				dij = dij + dist**2
			dij = dij**0.5 
			Distances.append(dij) 
			
			
	return Distances
