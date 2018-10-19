"""

Code to generate an RGG using a KD-Tree
based approach to speed up the nearest 
nieghbour computation.

Allows for different choices of boundary
conditions and node distributions.

Currently contains functions for generating
positions and storing them as data.

Also contains class for using constructing
and computing the properties of ensembles of
RGGs. 

"""

#Import packages:
import matplotlib.pyplot as plt
import numpy as np
import time
import networkx as nx
from scipy import special
import scipy
import scipy.spatial as sp
import itertools
from scipy import sparse
import math
import socket
import csv
import glob
import sqlite3
import os
import sys
import datetime


#Local packages:
import Network_Class as NetClass
import sql_utils as sql



def Make_RGG_Given_Positions(positions,r,boundary='S') :

	"""
	
	Generate an RGG with an algorithm using a KD tree (this algorithm is much faster than
	the naive O(N^2) approach since using a KD tree speeds upnearest nieghbour queries)
	
	We focus on RGGs in the unit square/hypercube ([0,1]^d).
	
	Parameters
	---------------
	
	positions : list
	
	Nxd List of lists containing the positions of the nodes.
	
	r : float 
	
	connection raidus for the RGG
	
	boundary : str (optional)
	
	String specifying whether to make the RGG with periodic
	or solid boundaries. 'P' - periodic boundaries. 'S' - solid boundaries.
	
	Note that choosing solid boundaries does not restrict the points to be in 
	a specific domain. 
	
	
	Returns
	---------------
	
	G : networkx graph
	
	"""

	N = len(positions)

	#Choose the boundaries:
	if boundary == 'P':
		tree = sp.cKDTree(positions, copy_data=True, boxsize=1.)
	elif boundary =='S' :
		tree = sp.cKDTree(positions)
	else :
		print("unknown boundary condition: error. Available options are 'P' or 'S'")


	n_pairs = tree.query_pairs(r)

	# Make a new networkx graph and add nodes:
	G = nx.Graph()
	for j in range(N) :
		G.add_node(j)

	# Add the edges in the list of pairs:
	for i in n_pairs:
		G.add_edge(i[0], i[1])

	return G



def Radius_From_Degree(N,Kappa,d) : 

	"""
	
	Compute the expected radius for a homogeneous RGG in
	d dimensions.
	
	Details about this formula can be found in: Dall, Jesper, and Michael Christensen. "Random geometric graphs." Physical review E 66.1 (2002): 016121.
	
	Parameters
	-------------
	
	N : int
	
	number of nodes in the graph. 
	
	Kappa : 
	
	Expected mean degree of the graph.
	
	d : int
	
	dimension.
	
	Returns
	-------------
	
	r : float
	
	connection radius.
	
	"""
	r = (1.0/((math.pi**0.5)) )*(( ((Kappa)/N)*scipy.special.gamma( (d +2.0)/2.0  )   )**(1.0/d ) )
	return r 


def degree_from_radius(N,r,d) : 
	return (  (N*np.pi**(d/2.0))  /( scipy.special.gamma( (d +2.0)/2.0  )) )*(r**d)

class Position_Data_Set :

	"""

	Class to keep track of the stored RGG positions
	for specific parameter sets.

	Includes methods to generate new positions.
	
	"""

	def __init__(self,N,d, directory="data/positions/",pattern="uniform") :
		
		"""
		
		Parameters
		---------------
		
		N : int 
		
		Number of positions
		
		d : int
		
		Dimension of the space
		
		directory : str (optional)
		
		Path to file where the positions are stored.
		
		pattern : str (optional)
		
		pattern or distribution of point process to sample from.
		
		Possible options:
		- "uniform" : uniform distribution in [0,1]^d
		- "gaussian" : Gaussian distribution with unit covariance. 
		
		"""
		
		
		# Make the directory if it does not exist:
		if not os.path.exists(directory):
			os.makedirs(directory)

		#Set class variables:
		self.N = N
		self.d = d
		self.directory = directory
		self.pattern = pattern
		
	
	def Generate_Positions(self):
		
		"""
		Generate positions in the specified pattern.

		Possible patterns:
		- "uniform" : uniform distribution in [0,1]^d
		- "gaussian" : Gaussian distribution with unit covariance.

		These are saved to file in a .csv file.

		Filename convention: pattern_N_d_number


		"""
		

		#Sample the positions:
		if self.pattern == "uniform":
			positions = np.random.uniform(0.0, 1.0, (self.N, self.d))
		elif self.pattern == "gaussian":
			positions = np.random.multivariate_normal(np.zeros(self.d), np.eye(self.d), size=self.N)
		else : 
			print("Error: unrecognized pattern. Available options:\n-`gaussian'\n-'uniform'")
			sys.exit(1)

		#Get the path to save generated positions to:
		path = self.Make_Position_Save_Name()

		#Write to file:
		with open(path, 'wb') as csvfile:
			writer = csv.writer(csvfile)

			# Write positions to all the rows:
			for i in positions:
				writer.writerow(i)

	def Make_Position_Save_Name(self,index=None):

		"""
		Construct a string for the filename of the .csv file
		to store the positions in given parameters.

		Filename convention: [pattern]_[N]_[d]_[number]

		This function will check the existing file names and make another with the
		next index.

		Need to make sure directoy is set correctly w.r.t to the folder we run the code in.


		Parameters
		----------------

		N : int

		Number of positions

		d : int

		Dimension

		pattern : str (optional)

		Specify the type of point pattern. Default is "uniform"


		Returns
		----------------

		fname : str

		Filename to save positions.


		"""

		#Find the number of existing files if the index has not been specified:
		if index is None:

			num_of_files = len( self.Enumerate_Position_File_Names() )
			index = num_of_files + 1

		#Make the file name of the .csv file:
		name = self.pattern + "_N{}_d{}_{:04d}".format(self.N, self.d, index)
		fname = self.directory + name + ".csv"
		
		return fname

	def Enumerate_Position_File_Names(self):
		
		"""

		Return the file names for the currently existing positions.

		Returns
		-------------
		
		Pos_File_Names : List of str
		
		List of strings containing names of .csv files containing
		positions in the directory associated with this class. 

		
		"""

		#Use glob to get the position names:
		Pos_File_Names = glob.glob(self.directory + self.pattern + "_N{}_d{}_*.csv".format(self.N, self.d))

		#Store this as a variable we can access in the class:
		self.Pos_File_Names = Pos_File_Names
		
		#Return the sorted list:
		return sorted(Pos_File_Names)

	def Load_Positions(self,index = 0 ):
		
		"""

		Load positions from specified file path.

		Parameters
		-------------

		index : int

		Index of positions. Integar values goes from 0 to the number of positions
		wit given parameters in the file. 

		Returns
		------------

		positions : list

		Set of positions.


		"""
		
		#This will note work if no positions exist so have to check:
		file_path = self.Enumerate_Position_File_Names()[index]
		print("loading {}".format(file_path) ) 
		
		#Empty list to store positions:
		positions = []

		with open(file_path, 'rb') as csvfile:
			reader = csv.reader(csvfile)

			# Read in positions while converting to float:
			for row in reader:
				positions.append([float(k) for k in row])

		return positions


class RGG(NetClass.Network) :

	"""
	
	Random Geometric Graph Class

	Inherits from network class and 
	also includes functions for saving
	data to an SQL database. 
	
	"""


	def __init__(self,N,r,d,boundary='S',Positions=None,pattern="uniform") :

		"""

		Generates and RGG.

		Some output needed to describe the files that we read in and out.

		Parameters
		----------------

		N : int

		Number of nodes in the network
		
		r : float
		
		Connection radius to use in the generation of the graph
		
		d : int
		
		Dimension of the embedding space
		
		boundary : str (optional)
		
		String specifying whether to make the RGG with periodic
		or solid boundaries. 'P' - periodic boundaries. 'S' - solid boundaries.
	
		Note that choosing solid boundaries does not restrict the points to be in 
		a specific domain. 
		
		Positions : list (optional)
		
		Can generate an RGG with a pre-specified set of positions. 


		"""

		#Load in the positions. if Positions=None then generate randomly.
		if Positions is	None :
			if pattern == "uniform":
				positions = np.random.uniform(0.0, 1.0, (N, d))
			elif pattern == "gaussian":
				positions = np.random.multivariate_normal(np.zeros(d), np.eye(d), size=N)

		#If positions is not none then load in the choosen positions:
		else :
			positions=Positions

		#Make the RGG:
		graph = Make_RGG_Given_Positions(positions,r,boundary=boundary)
		
		#Set some class variables:
		self.graph = graph
		self.dimension = int(d)
		self.boundary = boundary
		self.radius = r
		
		#Default name for table of data (used if we save to database)
		self.graph_data_table_name = "graph_properties"

		#Make the class:
		A = nx.to_scipy_sparse_matrix(graph)
		RGG = NetClass.Network.__init__(self,A)


	def get_single_entry_variables(self,Allowed_Types = {'float64', 'int64', 'float', 'str', 'int'}) :
		
		"""
		
		Returns the variables which can be stored in 
		a single cell within a database. 
		
		I.e. floats, ints and strings. 
		(Not Matrices or other large objects)
		
		Parameters
		--------------
		
		Allowed_Types :  set
		
		Set of strings specifying the variable types
		that we would like to return.
		
		Returns
		------------
		
		variable_dict : dict
		
		Dictionary with variable names as keys and their
		respective values as entries. 
		
		"""

		
		variable_dict = { }
		for i in vars(self).keys() :
			#first check that the stored varialbe is a float int (ie. not something dumb like a matrix or a class)
			if type(vars(self)[i]).__name__ in Allowed_Types :
				
				variable_dict[i] = vars(self)[i]


		return variable_dict

	def Save_to_SQL(self, database_name, extra_metadata=None,connect_to_server=True):

		"""
		Save the computed properties to an SQL database.

		Must find a way to fill empty columns in with
		zeros or missing data.

		For generation params such as dimension will need to save somewhere else.

		Parameters
		------------------
		
		database_name : str
		
		Name of the database to save to.
		(not used if connect_to_server == True )
		
		extra_meta_data : dict
		
		Dictionary containing data about graph creation such 
		as the dimension which are of interest but not recorded
		as variables assoicated with the network class.
		(Can also record the time of creation etc..). 
		
		connect_to_server : bool (optional)
		
		Set to true if the database of interest is hosted on
		and external server. 

		"""


		#Make a dictionary which contains all variables of the class
		#which are either float int or string. 
		variable_dict = self.get_single_entry_variables()
		
		#Add any extra metadata to the dictionary to be saved:
		for k in extra_metadata.keys() : 
			variable_dict[k] = extra_metadata[k]



		#Noisy saving of variables:
		print("saving variables:")
		print(variable_dict)
		
		#Must make the SQL table if it doesn't exist:
		sql.Make_SQL_Table(database_name,variable_dict,self.graph_data_table_name,connect_to_server=connect_to_server)

		#Save the data:
		sql.Save_Dict_to_SQL(database_name,variable_dict,self.graph_data_table_name,connect_to_server=connect_to_server)


def look_for_max_index(database_name,parameter_dict) :

	"""
	Look for indices assoicated with
	stored positions in a database containing
	data about RGGs.
	
	"""
	
	#Make a graph ensemble:
	current_graph_data = Graph_Ensemble(database_name,parameter_dict)

	#Call the function ot get the list of position indices:
	current_position_indices = current_graph_data.get_prop_array('position_index')

	#Find the largest current index associated with graphs in database.
	current_max_index = len(current_position_indices)
	print("current_max= {}".format(current_max_index) ) 
	print("current indices = {}".format(current_position_indices) ) 
	
	return current_max_index



def Sample_RGGs(N, r, d, boundary,degree_scaling_parameter,Samples, methods_to_call, database_name,graph_data_table_name,pattern,Save_Positions=True,connect_to_server=True) :
	"""

	Sample RGGs multiple RGGs from the same ensemble
	and save their properties to an SQL database.

	To do:
	
	-Note: we can use the command dir(rgg_graph) to get 
	the methods for a particular class.
	
	
	Parameters
	-----------------
	
	N : int
	
	Size of the graph
	
	r : float
	
	Connection raidus.
	
	d : int
	
	Dimension of the embedding space.
	
	boundary : str
	
	Choice of boundary condition. 
	Current options 'P' (periodic) or 'S' (solid). 
	
	degree_scaling_parameter : float
	
	Scaling parameter for degree. This is included so that it
	can be saved in the table of data about the graph. 
	
	Samples: int
	
	Number of RGGs to sample
	
	methods_to_call : list of str
	
	List of strings where each string is a method within
	network_class or RGG class which can be called to compute
	network properties.
	
	database_name : str
	
	path to the database where the data will be saved.
	
	pattern : str
	
	pattern of the point process. Allowed values are:
	"uniform" and "gaussian". 
	
	
	"""
	
	#Make a position data class if we are saving the positions:
	if Save_Positions == True : 
	
		#Check how many positions we already have:
		Pos_Data_Class = Position_Data_Set(N, d,pattern=pattern)
		current_number_stored = len(Pos_Data_Class.Enumerate_Position_File_Names() ) 
	
		#Check which position indices are already in the database for these params:
		parameter_dict = { 'N' : N , 'd' : d , 'boundary' : boundary , 'r' : r } 
		
		#This only works if graphs have already been saved in the database: 
		try : 
			current_max_index = look_for_max_index(database_name,parameter_dict)
		except : 
			current_max_index = 0 
	
		#Sample more positions:
		for k in range(Samples):
			Pos_Data_Class.Generate_Positions()

	
	#Loop for sampling the graphs:
	for p in range(Samples):

		start_time = time.time()
		metadata = { }
		metadata['sim_start_time'] = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
		print("sample={}/{}".format(p,Samples))
		
		if Save_Positions == True :
		
			#Set the position index to generate as the current max:
			#(minus one to account for pytonic indiexing)
			pos_index = p+current_max_index
			
			# Load in positions:
			positions = Pos_Data_Class.Load_Positions(index=pos_index)
		
			# Make an RGG class object: (passing in params which define the ensemble)
			rgg_graph = RGG(N, r, d, Positions=positions, boundary=boundary)
			
			metadata['position_index'] = pos_index
			metadata['degree_scaling_parameter'] = degree_scaling_parameter

		else : 
			rgg_graph = RGG(N, r, d, boundary=boundary,pattern=pattern)
			metadata['degree_scaling_parameter'] = degree_scaling_parameter

		rgg_generation_time = time.time() - start_time
		metadata['rgg_generation_time'] = rgg_generation_time

		#Loop to compute the chosen properties of the graph:
		for method in methods_to_call:
			function = getattr(RGG, method)
			val = function(rgg_graph)
			
		#time taken to make the RGG and call all the methods of interest:
		total_sample_time = time.time()-start_time
		metadata["total_sim_time"] = total_sample_time
		metadata["sim_end_time"] =  datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
		
		# Save stuff to the database: (make sure we save the position index).
		rgg_graph.graph_data_table_name = graph_data_table_name
		
		#Save to database:
		rgg_graph.Save_to_SQL(database_name,extra_metadata= metadata,connect_to_server=connect_to_server)


		print("time taken= {} s".format(total_sample_time))



class Graph_Ensemble : 
	
	
	"""
	
	An RGG ensemble is a collection of graphs with the same
	set of parameters. This class reads in data about these 
	ensembles which has been precomputed and saved into an SQL
	database.
	
	This will be used in the analysis after the data has been generated.  
	
	"""

	def __init__(self,database_name,graph_data_table,parameter_dict,verbose=True ,connect_to_server=False) : 
		
		"""
		
		Parameters
		---------------
		
		database_name : str
		
		Path to the SQL datbase where the graph data
		is stored. 
		Not used if connect_to_server == True
		
		graph_data_table : str
		
		Name of the table within the database where graph
		data is stored.
		
		parameter_dict : dict
		
		Dictionary containing the parameter set
		for this specific ensemble. 
		
		Keys are: 'boundary' , 'N' , 'd' and 'r'
		which are the parameters for RGGs in [0,1]^d. 
		
		verbose : bool (optional)
		
		Set this to true if we want to print info about
		the nubmer of ensemble members and the properties
		of the graph.
		
		connect_to_server : bool
		
		Set to true if we are reading from a database connected to
		some exteneral server. 
		
		
		"""
		 
		self.parameter_dict = parameter_dict
		self.N = parameter_dict['N']
		self.boundary = parameter_dict['boundary']
		self.d = parameter_dict['d']
		self.radius = parameter_dict['r']
		self.degree_scaling_parameter = parameter_dict['degree_scaling_parameter']
		self.database_name = database_name
		
		#tuple used when reading in data:
		self.values = (self.boundary,self.N,self.d,self.radius)
		
		#Default name for table of data:
		self.graph_data_table_name = graph_data_table

		#Parameter to control whether SQL database is hosted on a server via MySQL
		self.connected_to_server = connect_to_server
		
		#Find the # of ensemble members by reading from the Network_Size column: 
		graph_numbers = self.get_prop_array("Network_Size")
		num_of_members = len(graph_numbers)
		self.num_of_members = num_of_members


		
		#Print info about the ensemble if required:
		if verbose ==True : 
		
			print("Found {} ensemble members in {} with (bc,N,d,r)=({},{},{},{})".format(num_of_members,database_name,self.boundary,self.N,self.d,self.radius) )
			column_names = sql.find_column_names(self.database_name,self.graph_data_table_name,connect_to_server=connect_to_server)
			
			#List of which table headers (ie. different graph properties) that we have?
			#print("Found the following graph properties:")
			#for k in column_names:
			#	print("-{}".format(k))

	def get_prop_array(self,required_prop) :
		
		"""
		
		Get the array of values associated with
		a certain property.
		
		Parameters
		--------------
		
		required_prop : str
		
		String specifying the graph property required. This will be one of the
		Class variables associated with the network class. e.g:
		'Algebraic_Connectivity' , 'LCC_Min_Deg'
		
		"""
		
		#Exit if we don't find the property of interest in the database:
		column_names = sql.find_column_names(self.database_name,self.graph_data_table_name,connect_to_server=self.connected_to_server)

		if required_prop not in column_names : 
			print("\nError: property not found in database.\nFailed to find {} in:\n{}".format(required_prop,column_names) )
			sys.exit(1)
		

		#open the connection to the SQL database:
		#We use differetn syntax for MySQL and sqlite3 databases
		if self.connected_to_server==True :
		
			connection, cursor = sql.connect_to_maths_cluster_db()
			cursor.execute("SELECT {} FROM {} WHERE boundary=%s AND Network_Size=%s AND dimension=%s AND radius=%s".format(required_prop, self.graph_data_table_name), self.values)


		else : 
			#Open the connection to the databse:
			connection = sqlite3.connect(self.database_name)
			cursor = connection.cursor()
		
			#Pull out the array with the specified properties: (sqlite syntax)
			cursor.execute("SELECT {} FROM {} WHERE boundary=? AND Network_Size=? AND dimension=? AND radius=?".format(required_prop,self.graph_data_table_name),self.values)


		result = cursor.fetchall()
		property_array = [i[0] for i in result]
		
		#close the connection:
		connection.close()
		
		
		return property_array
		
	def get_stat(self,required_prop,statistic_func) : 
		
		
		"""
		
		Gets some statistic of interest for one of the
		graph properties.
		
		Parameters
		--------------
		
		required_prop : str
		
		String specifying the graph property of interest. This will be one of the
		Class variables associated with the network class. e.g:
		'Algebraic_Connectivity' , 'LCC_Min_Deg'
		
		statistic_func : function
		
		Some function which operates on lists of values
		in order to produce.
		
		e.g np.mean , np.var. 
		
		
		"""

		return statistic_func( self.get_prop_array(required_prop) ) 
		
	
	def bivar_stat(self,property_1,property_2,bi_var_stat_func) : 
	
		"""
		
		Return a statistic of two properties. For example
		the correlation between the variables.
		
		Parameters
		------------
		
		property_1 : str
		
		First property of interest
		
		property_2 :str
		
		Second property of interest
		
		bi_var_stat_func : function
		
		Some function which takes in two arrays
		and computes a statistic.
		
		E.g np.corrcoef , scipy.stats.pearsonr
		
		"""
		
		return bi_var_stat_func( self.get_prop_array(property_1) , self.get_prop_array(property_2) ) 



class multi_ensemble : 

	"""
	This class is used when analysing a range
	of RGG ensembles with different parameters.
	Allows slices through parameter space.
	
	"""
	
	def __init__(self,database_name,graph_data_table,parameter_dict_list,connect_to_server=False) :
	
		"""
		
		Parameters
		--------------
		
		database_name : str
		
		String specifying file path to the database where the graph parameters are stored
		
		graph_data_table : str
		
		Name of the table in the database to load graph data from.
		
		parameter_dict_list : list of dict
		
		List of dictionaries containing the parameters of interest.
		
		"""
	
		self.Ensembles = [ Graph_Ensemble(database_name,graph_data_table,k,connect_to_server=connect_to_server) for k in parameter_dict_list ]
		
	
	def pull_out_prop_stat_as_function( self, required_prop , statistic_func , sweep_parameter , second_property = None ,  fixed_param_dict = None ) :
	
		"""
		
		Sweep through the parameter space pulling out a chosen
		statistic of a given graph metric as a function of the 
		chosen parameter.
		
		Includes the option to include a second property for which we 
		can compute the bivariate function of.
		
		Parameters
		--------------
		
		required_prop : str
		
		string specifying the graph metric to compute
		
		statistic_func : function
		
		Some function which operates on lists of values
		in order to produce.
		Can be a function of two lists if the `second_property'
		is supplied. 
		
		e.g np.mean , np.var. 
		
		
		sweep_parameter : str
		
		Allowed values: N ,d ,r , boundary

		
		second_property : str (optional)
		
		Specify a second property to compute. 
		
		fixed_param_dict : dict (optional)

		Dictionary containing other parameters which we want to
		keep fixed. In the format:

		e.g { 'boundary' : 'P' } or { 'boundary' : 'S' , 'degree_scaling_param' : 1.5 } 
		
		
		
		Returns
		------------
		
		parameter_array : list
		
		List of values of the sweep parameter found in the 
		ensemble
		
		statistic_array : list
		
		List of values of the chosen statistic of the distribution
		of required_prop for each of the ensembles.
		
		"""

		#Restrict to the set of ensembles where the param is fixed:
		if fixed_param_dict is not None : 
		
			Required_Ensembles = [ ]
			for ens in self.Ensembles :
				
				#Check all elements in the array are true:
				truth_array = [ ens.parameter_dict[fixed_param_dict.keys()[i]] == fixed_param_dict.values()[i] for i in range(len(fixed_param_dict)) ] 
				
				if all(truth_array) == True :  
					#print("adding {}".format(ens.parameter_dict) ) 
					Required_Ensembles.append(ens)
					
					
		else: 
			Required_Ensembles = self.Ensembles
		
		#Check whether to compute a second property:
		if second_property is not None : 
			statistic_array = [ ens.bivar_stat(required_prop,second_property,statistic_func) for ens in Required_Ensembles ]
		else : 
			statistic_array = [ ens.get_stat(required_prop,statistic_func) for ens in Required_Ensembles ]
			
		parameter_array = [ ens.parameter_dict[sweep_parameter] for ens in Required_Ensembles ]
		
		print("Number of samples data:")
		print("First_dict = {}".format(ens.parameter_dict))
		print( [ ens.num_of_members for ens in Required_Ensembles ] ) 



		return parameter_array , statistic_array
		


class graph_plot :

	"""
	Class for plotting graphs. Allows the attrbiutes to be
	edited if required. 
	
	Focus on spatial graphs. 
	
	Makes use of networkx's plotting function.
	
	"""
	

	
	
	def __init__( self , graph , positions = None ) :
	
		"""
	
		Parameters 
		---------------
	
		graph : networkx graph
	
		the graph for plotting
	
		positions : list 
	
		set of positions. Currently only supports plotting in 2D. 
	
		"""
		
	
		self.graph = graph
		if positions is not None and positions != 'spectral':
			
			#If we have more than 2 position coords then take only the first two:
			if len(positions[0]) != 2 : 
				print("More than 2d found. Taking first two dimensions")
				self.positions = [  [ i[0] , i[1] ]  for i in positions ] 
			else : 
				self.positions = positions
		
		#If positions are not specified use the spring layout:
		elif positions == None : 
			self.positions=nx.spring_layout(graph)
		
		elif positions == 'spectral' :
			self.positions=nx.spectral_layout(graph)
		
		
		
		
		#Set default parameters:
		self.node_size = 35
		self.node_line_width = 0.1
		self.color_scheme = 'jet'
		self.edge_width = 1.0
		self.node_color = 'b'
		self.node_shape = 'o'
		self.node_transparency = 1.0
		self.edge_color = 'k'
		self.edge_style = 'dashed' #solid|dashed|dotted,dashdot
		self.edge_transparency = 1.0
		
	
	def make_plot(self) : 
		nx.draw_networkx_nodes(self.graph, self.positions ,  node_size = self.node_size, node_color = self.node_color,node_shape = self.node_shape,alpha = self.node_transparency, cmap = self.color_scheme, linewidths = self.node_line_width)
		plt.axis('off')
		nx.draw_networkx_edges(self.graph, self.positions, width = self.edge_width,edge_color = self.edge_color,style = self.edge_style,alpha =self.edge_transparency)
	
	
	def Add_Edges(self,graph) :	
		plt.axis('off')
		nx.draw_networkx_edges(graph, self.positions, width = self.edge_width,edge_color = self.edge_color,style = self.edge_style,alpha =self.edge_transparency)
		
	def save_plot(self,file_path,form="pdf") :
		plt.figure(figsize=(10,10))
		self.make_plot()
		plt.savefig(file_path + ".{}".format(form) , bbox_inches = 'tight',format=form)
		
	def show_plot(self) :
		plt.figure(figsize=(10,10))
		self.make_plot()
		plt.show()
		
	def change_node_size(self,size) :
		self.node_size = size
