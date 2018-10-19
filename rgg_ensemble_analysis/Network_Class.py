"""

Defines a class for networks which can 
be used to compute network properties
and return network matrices. 

"""

#Import required packages:
import numpy as np
from scipy import sparse
import networkx as nx
import Sparse_Eigenvalues as SEV
import matplotlib
import LCC_Extract as LCC
matplotlib.use('Agg')


def networkx_to_netclass(G) :
	
	"""
	
	Function to convert a networkx 
	graph object to a network class.
	
	"""
	
	A = nx.to_scipy_sparse_matrix(G)
	
	return Network(A)


class Network :
	
		
	def __init__( self , A ) :

		"""

		Class to store information about network and compute relevant properties.

		We store the adjacency in sparse format in order to reduce memory requirements
		and only convert to dense matrix format when required to run certain functions
		(e.g converting to networkx graph format)


		Parameters
		-------------

		A: Scipy Sparse CSR matrix.

		Adjacency matrix of an undirected graph.

		"""
		
		#Define adjacency matrix and size of network:
		self.Adjacency = A
		self.Network_Size = A.shape[0]

	def Mean_Degree(self) :
		Degree_array = np.array( self.Adjacency.sum(axis=0) )
		self.mean_degree = np.mean( Degree_array[0] ) 
		return self.mean_degree


	def Min_Degree(self) :
		Degree_array = np.array( self.Adjacency.sum(axis=0) ) 
		self.Min_Degree = min( Degree_array[0]  )
		return self.Min_Degree
		
	def Max_Degree(self) : 
		Degree_array = np.array( self.Adjacency.sum(axis=0) ) 
		self.Max_Degree = max( Degree_array[0]  ) 
		return self.Max_Degree
		
	def Edges(self) :
		self.Edges = len( self.Adjacency.data)/2
		return self.Edges
		
	def Degrees(self) :
	
		"""
		
		Returns
		----------
		
		degrees :  list
		
		A list containing the degrees of the network.
		
		"""
		degrees = np.array( self.Adjacency.sum(axis=0) )[0]
		return degrees

	def Laplacian(self) :
	
		"""
		
		Constructs the laplacian matrix of the network
		
		We use the unnormalized Laplacian matrix which is defined as:
		
		L = D - A ,
		
		where D is the diagonal matrix of degrees and A is the adjacency matrix of the network 
		
		
		Returns
		----------
		
		L : scipy sparse matrix
		
		Laplacian matrix for the network
		
		"""
	
	
		degrees = np.array( self.Adjacency.sum(axis=0) )
		Is = np.arange(0,self.Network_Size,1)
		D = sparse.coo_matrix((degrees[0], (Is, Is)) , shape = (self.Network_Size,self.Network_Size) , dtype=float).tocsr()
		L = - self.Adjacency + D
		return L
		

	def Mean_Path_Length(self) :
		
		"""
		
		Uses the function from networkx to compute the average shortest path
		length of the network.
		
		Average shortest path length is computed for the largest connected component as
		it is undefiend for unconnected graphs.
		
		Returns
		---------
		
		Path_Len : float 
		
		Average shortest path length between nodes in the network.
		
		"""

		#Check if we have already extracted the LCC:
		if 'A_LCC' not in vars(self).keys() :
			self.A_LCC = self.Largest_Connected_Component()

		#Use networkx to get the average shortest path length:
		G=nx.from_numpy_matrix(self.A_LCC.todense())
		Path_Len = nx.average_shortest_path_length(G)
		
		self.Mean_Path_Length = Path_Len
		return Path_Len

	def Diameter(self):
		# Check if we have already extracted the LCC:
		if 'A_LCC' not in vars(self).keys():
			self.A_LCC = self.Largest_Connected_Component()

		# Use networkx to get the average shortest path length:
		G = nx.from_numpy_matrix(self.A_LCC.todense())
		Diameter = nx.diameter(G)

		self.Diameter = Diameter
		return Diameter

		
	def Components(self) :
	
		"""

		Return a list of the connected components in the network
		
		Returns
		-------------

		Components : list of lists

		List of lists containing labels of node nodes in
		each of the respective components in the graph.


		"""
		Components = LCC.Components(self.Adjacency)

		return Components
		
	def Largest_Connected_Component(self, positions = None ) :
	
		"""

		Function to get the adjacency matrix of the largest connected
		component of the graph.
		
		Parameters
		---------------
		
		positions : list (optional)

		If this argument is provided then we also return the positions of the nodes
		in the largest connected component.


		Returns
		---------------
		
		A_LCC : scipy sparse matrix

		Adjacency matrix of the LCC for the network.

		LCC_Positions : list

		positions of the nodes in the LCC
		
		"""

		print("Extracting LCC")

		#Get components in the graph and return the largest one:
		index, LCC_Vec = max(enumerate(self.Components() ), key = lambda tup: len(tup[1]))

		#Convert to a matrix:
		A_LCC = self.Adjacency[LCC_Vec, :].tocsc()[:, LCC_Vec]
		A_LCC = A_LCC.tocsr()

		#Store the adjacency of the LCC and the number of nodes in it
		#as class variables:
		self.A_LCC = A_LCC
		self.N_LCC = A_LCC.shape[0]
		self.F_LCC = float(self.N_LCC)/float(self.Network_Size)
		
		
		
		#Extract the positions if we require them:
		if positions is not None : 
			LCC_Positions = [ ]
			for i in LCC_Vec :
				LCC_Positions.append(positions[i])
			
			return A_LCC , LCC_Positions
		
		else :
			return A_LCC

	def Fiedler_Largest_Connected_Component(self, positions=None):

		"""

		Extract the LCC using the fiedler vector to find the
		null space of the laplacian.

		Parameters
		--------------

		return_positions : bool (optional)

		Specify whether we also return the positions of the nodes in the LCC.

		Returns
		--------------

		Adjacency matrix corresponding to the largest connected component of the network.

		"""

		Mu2, Vector = SEV.Smallest_Non_Zero( self.Laplacian(), return_eigenvector=True)

		LCC_Comp_Vec = [  ]
		
		#Get values close to non zero:
		for i in range(len(Vector)) :
			if Vector[i] < 0.5 :
				LCC_Comp_Vec.append(i)

		Fiedler_of_LCC = Vector[LCC_Comp_Vec]
		A_LCC = self.Adjacency[LCC_Comp_Vec, :].tocsc()[:, LCC_Comp_Vec]
		A_LCC = A_LCC.tocsr()
		N_LCC = len(LCC_Comp_Vec)

		# Extract the positions if we require them:
		if positions is not None:
			LCC_Positions = []
			for i in LCC_Comp_Vec:
				LCC_Positions.append(positions[i])

			return A_LCC, LCC_Positions , Fiedler_of_LCC

		else:
			return A_LCC
	
	
	def LCC_Network(self, positions = None ) :
		"""
		Call the network class object for the 
		largest connected component of the network.

		Also return the corresponding positions if required.


		Uses the fiedler vector based method to extract the LCC
		from the network.

		Must return fiedler vector of the LCC rather than the original fiedler
		vector.

		"""
		if positions is not None :
			#LCC_Adjacency, LCC_Positions , Fiedler = self.Fiedler_Largest_Connected_Component(positions=positions)

			LCC_Adjacency, LCC_Positions = self.Largest_Connected_Component(positions=positions)

			#return Network(LCC_Adjacency), LCC_Positions , Fiedler


			return Network(LCC_Adjacency), LCC_Positions
		else :
			LCC_Adjacency = self.Largest_Connected_Component()
			return Network(LCC_Adjacency)


	def LCC_Min_Degree(self):
		
		"""
		Return the minimum degree of the
		Laplacian matrix.
		
		"""
		
		
		
		# Check if we have already extracted the LCC:
		if 'A_LCC' not in vars(self).keys():
			self.A_LCC = self.Largest_Connected_Component()
		Degree_array = np.array(self.A_LCC.sum(axis=0))
		self.LCC_Min_Degree = min(Degree_array[0])
		return self.LCC_Min_Degree

	def LCC_Laplacian(self) :
		
		
		"""
		
		Return the laplacian of the LCC for the network. 
		
		"""

		# Check if we have already extracted the LCC:
		if 'A_LCC' not in vars(self).keys():
			self.A_LCC = self.Largest_Connected_Component()

		degrees = np.array( self.A_LCC.sum(axis=0) )
		Is = np.arange(0,  self.A_LCC.shape[0],  1)
		D = sparse.coo_matrix((degrees[0], (Is, Is)) , shape = (self.A_LCC.shape[0],self.A_LCC.shape[0]) , dtype=float).tocsr()
		
		return  -self.A_LCC + D
		
	def Largest_Comp_Fraction(self) :
		"""Fraction of the Nodes in the Largest connected component"""
		return float(self.Largest_Connected_Component().shape[0])/float(self.Adjacency.shape[0])

		
	def Algebraic_Connectivity(self) :
	
		"""
		
		Compute the algerbaic connectivity of the LCC of
		the network of interest. Uses a sparse linear algebra method.

		The algebraic connectivity is the second smallest eigenvalue
		of the Laplacian matrix of the network.

				
		"""

		#Get the graph laplacian for the LCC:
		LCC_Laplace = self.LCC_Laplacian()

		#Using the original mathod:
		self.algebraic_connectivity = SEV.Smallest_Non_Zero(LCC_Laplace)

		return self.algebraic_connectivity
		
	def Fiedler_and_Mu2(self) :
	
		"""
		Return the algebraic connectivity and fiedler vecor for the network.
		
		Assumes LCC as it is not useful to compute Fiedler partition if the graph is not fully
		connected.
		"""

		Laplacian = self.LCC_Laplacian()
		#Laplacian = self.Laplacian()
		
		Mu2, Fiedler = SEV.Smallest_Non_Zero( Laplacian , return_eigenvector = True )
		self.Fiedler_Vector = Fiedler
		
		
		return  Mu2 , Fiedler
	
	
	def Fraction_In_Fiedler_Partition(self) :
		"""
		Return the fraction of nodes in the smaller partition associated with
		taking positive and negative components of the Fiedler vector. 
		
		Can take a pre-computed Fiedler vector as an optional argument.
		"""
		
		#Compute the Fiedler vector if we have not previously:
		if 'Fiedler_Vector' not in vars(self).keys():
			Mu2 , Fiedler = self.Fiedler_and_Mu2()
		else : 
			Fiedler = self.Fiedler_Vector
		
		Num_In_Partition = np.sum( [ 0.5*(np.sign(i) + 1.0) for i in Fiedler  ])
		Fraction_in_smaller_partition = min (  Num_In_Partition , len(Fiedler) - Num_In_Partition )/len(Fiedler)
		self.Fraction_In_Fiedler_Partition = Fraction_in_smaller_partition
		return Fraction_in_smaller_partition
	
	def Largest_Laplacian_Eigenvalue(self) : 
		return SEV.Largest_EV( self.Laplacian() )
	
	def Largest_Adjacency_Eigenvalue(self) :
		return SEV.Largest_EV( self.Adjacency ) 
		
	def Adjacency_Spectrum(self) :
		return np.linalg.eigvals( self.Adjacency.todense() )
		
	def Laplacian_Spectrum(self) :
		return np.linalg.eigvals( self.Laplacian().todense() )
		
	def Normed_Laplace_Spectrum(self) :
		return np.linalg.eigvals(  self.Normalized_Laplacian().todense() )


   
