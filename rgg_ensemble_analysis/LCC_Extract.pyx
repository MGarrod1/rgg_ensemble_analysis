"""

Code to extract the largest connected component 
for a network specified in terms of a sparse adjacency
matrix. 

Adjacency map class adapted from code by Till Hoffmann.

See: http://tillahoffmann.github.io/2016/04/18/Cpp-containers-in-cython.html

"""


#Import required packages:
import time
import numpy as np
from scipy.spatial.distance import pdist
from libc.math cimport abs
from scipy import sparse
from scipy.sparse.linalg import eigsh
import scipy
import numpy as np
from scipy import special
from libcpp cimport bool
from numpy.math cimport INFINITY
import time

cdef extern from "<algorithm>" namespace "std":
	iter find [iter, T](iter first, iter last, const T& val)

#Initialise a string (this is needed for compilation to work) 
STUFF = "HI"

# Import the map and vector templates from the STL
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector
from libcpp.utility cimport pair as cpp_pair

ctypedef cpp_vector[int] cpp_neighbourhood
ctypedef cpp_map[int, cpp_neighbourhood] cpp_adjacency_map
ctypedef cpp_pair[int, cpp_neighbourhood] cpp_item

# Import a few operators because they aren't supported by cython syntax
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc

#Define adjacency map class for detecting the giant component:
cdef class AdjacencyMap:
	
	cdef:
		cpp_adjacency_map container
		
	def __init__(self, int[:, :] adjacency_list):
		cdef:
			int i, ego, alter
			cpp_neighbourhood neighbourhood
			
		# Iterate over all entries of the adjacency list
		for i in range(adjacency_list.shape[0]):
			ego = adjacency_list[i, 0]
			alter = adjacency_list[i, 1]
			
			# Check if the ego is already in the map 
			# (see http://stackoverflow.com/a/101980/1150961 for details)
			lb = self.container.lower_bound(ego)
			
			# Check if the key already exists
			if lb != self.container.end() and ego == deref(lb).first:
				# Add the node to the pair
				deref(lb).second.push_back(alter)
			else:
				# Insert a new key value pair
				neighbourhood = cpp_neighbourhood()
				neighbourhood.push_back(alter)
				self.container.insert(lb, cpp_item(ego, neighbourhood))
				
	def get(self, int ego):
		"""
		Get the neighbours of `ego` or `None` if `ego` isn't in the map.
		"""
		# Search the dictionary
		iterator = self.container.find(ego)
		# Return none if we didn't find anything
		if iterator == self.container.end():
			return None
		# Create a list of values from the vector
		values = []
		# Iterate over the neighbourhood and add values to the list
		neighbourhood = deref(iterator).second
		neighbourhood_iterator = neighbourhood.begin()
		while neighbourhood_iterator != neighbourhood.end():
			values.append(deref(neighbourhood_iterator))
			preinc(neighbourhood_iterator)
			
		return values
	
	def _get_many(self, int ego, int repeats):
		"""
		Simple function to illustrate overhead.
		"""
		cdef int i
		# Try to find the ego a large number of times
		for i in range(repeats):
			iterator = self.container.find(ego)

def Expander(int x,Edge_List) :

	"""
	
	Function to idenfity all of the nodes in the same
	component as node x. 
	
	
	Parameters
	------------
	
	x : int
	
	Label for the node
	
	Edge_List :  List
	
	List of tuples specifying the edges in the network
	Takes the form: [ ( 0 1 )  ( 0 2 ) ( 1 0 ) ( 2 0 ) ] 
	
	
	Returns 
	--------
	
	Component_List : List
	
	List of the labels of nodes within the same connected component as node i
	
	
	"""
	
	stl_adjacency_map = AdjacencyMap(Edge_List)
	
	
	cdef cpp_vector[int] Component_List  
	cdef int change = 1 
	cdef int j 
	cdef int i 
	
	Component_List.push_back(x)

	cdef cpp_vector[int] Vec  
	
	#Keep Finding Niehgbours of nodes until
	#the component size remains constant.
	while change >= 1   :
		len_1 = len(Component_List)
		
		Current_array_store = Component_List
		
		#Cycle through node labels in the current component:
		for j in range( len(Component_List) ) :

			current_node = Component_List[j]
			

			#If the node has niehgbours:  
			iterator = stl_adjacency_map.container.find(current_node)  
			
			
			#Check that the iterator isn't emtpy: 
			if iterator != stl_adjacency_map.container.end() : 

				#Cycle through the neighbour list checking which ones are already in the neighbour list:
				for i in deref(iterator).second :

					#If they are not already then add them:
					if find(Component_List.begin(),Component_List.end(),i) == Component_List.end() : 
						Component_List.push_back(i)   
				
		#Compute the size of the new component with extra neighbours:
		len_2 = len(Component_List)
		change = len_2 - len_1
	
	return Component_List
	
	

def Components(A) : 

	"""
	Function to obtain the different connected
	components within a graph. 
	
	Parameters
	-----------
	
	A : Sparse adjacency matrix
	
	Returns
	----------
	
	comp :  List
	
	List of Lists containing the indices of the nodes contained
	in the different components in the network
	
	"""
	
	
	N = A.shape[0]
	
	#Begin by turning sparse matrix into tuples:
	#Each edge will be represented as a tuple
	Edge_List = np.transpose(np.nonzero(A)).astype(np.int32) 
	
	#Boolean to which is true if nodes are in the list of nodes
	#already visited.
	T_F = False
	
	#Empty array containing a list of the nodes which are already specified within a component:
	Done_List = [  ]
	
	#Empty array containing the list of components:
	Comps = [ ] 
	
	cdef int i 
	
	#go through nodes and do the expander
	for i in range(N) :
		
		
		T_F = i in Done_List
		
		if T_F == False :
			#Given a node label and list of edges we 
			#can identify the component which that node belongs to
			Current_Comp = Expander(i,Edge_List) 
				  
			Done_List = Done_List + Current_Comp 
			Comps.append( Current_Comp )
		
	return Comps


