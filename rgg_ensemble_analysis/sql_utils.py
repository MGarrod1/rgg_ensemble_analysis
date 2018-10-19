"""

Functions to be used for reading
and writing data to SQL databases

"""
import sqlite3
import sys
import pdb
import pickle

"""If importing MySQLdb doesn't work then we
must isntead import the connector"""
try : 
	import MySQLdb
except : 
	print("Importing mysql as MySQLdb")
	import mysql.connector as MySQLdb


def connect_to_maths_cluster_db() : 

	#Details for accessing the database on the maths cluster:
	file_name = "maths_clust_database_details"
	
	#Load details in from a pickled dictionary
	with open(file_name + '.pickle', 'rb') as handle:
		details = pickle.load(handle)
	
	host = details["host"]
	user = details["user"]
	passwd = details["passwd"]
	db = details["db"]

	connection = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db)
	cursor = connection.cursor()
	
	return connection , cursor



def find_sql_type_string(variable) :
	
	"""
	
	Given a python variable return the
	sql type string appropriate for saving it

	Current variable lengths used for SQL e.g VARCHAR(10)
	are arbitary for now. 
	
	e.g if we want to save an int nto the SQL database

	Parameters
	----------------

	variable : some variable in python

	Returns
	----------------

	string specifying the variable type
	to save this as in an SQL database.
	
	"""


	if type(variable).__name__ == 'str' :
		return "VARCHAR(30)"
	elif type(variable).__name__ == 'float64' or type(variable).__name__== 'float' :
		#return "FLOAT(5,4)"
		return "DOUBLE"

	elif type(variable).__name__ == 'int' or type(variable).__name__ == 'int64' :
		return "INT(20)"
	else:
		#Exit if we detect a variable which is not one of supported types.
		print("Error: Found a {}\nCurrently only supports int, str, float data types".format(type(variable)) ) 
		sys.exit(1)

def convert_to_supported_type(variable) :
	"""
	Converts variables into types which are
	supported by MySQL.
	
	Types such as float64 are not supported
	by SQL. 
	
	Inspired by:
	https://stackoverflow.com/questions/17053435/mysql-connector-python-insert-python-variable-to-mysql-table
	
	"""
	
	if type(variable).__name__ == 'float64' :
		return float(variable)
		
	elif type(variable).__name__ == 'int64' :
		return int(variable)
		
	else :
		return variable
	
	

def Make_SQL_Table(database_name,variable_dict,table_name,connect_to_server=True):
	
	"""
	
	Make an SQL database containing the graph
	properties which we are interested in.

	Parameters
	-----------------

	database_name : str

	String specifying the name and location of the SQL database.

	variable_names : dict

	Dictionary containing variable types as keys. 
	
	table_name : str
	
	Name of the table to create within the database. 
	
	connect_to_server : Bool (optional)
	
	True is the database is hoted on some server. 

	"""
	
	if connect_to_server==True : 
		
		connection, cursor = connect_to_maths_cluster_db()
	
	else : 
		#Open the connection to the databse:
		connection = sqlite3.connect(database_name)
		cursor = connection.cursor()
	
	#Alernative command if we don't have an integar primary key:
	sql_command = """CREATE TABLE IF NOT EXISTS {} (""".format(table_name)
	

	#Loop through the variables and their types in order to build the table.
	for k in variable_dict.keys() :
		sql_command += k + " " + find_sql_type_string(variable_dict[k]) + ","
	#remove final comma:
	sql_command = sql_command[:-1]
	
	sql_command = sql_command + ");"

	cursor.execute(sql_command)

	connection.commit()

	connection.close()


def find_column_names(database_name,table_name,connect_to_server=True) : 
	
		"""
		
		Return the names of the columns of a specific table
		in the database of interest.
		
		Parameters
		--------------
		
		database_name : str
		
		file path to the databse of interest
		
		tabel_name : str
		
		name of the table within the database. 
		
		Returns
		-------------
		
		column_names : list of str
		
		Containing the names of the different graph properties
		stored in the SQL database. 
		
		"""
		
		
		if connect_to_server==True : 
		
			connection, cursor = connect_to_maths_cluster_db()
	
		else : 
			#Open the connection to the databse:
			connection = sqlite3.connect(database_name)
			cursor = connection.cursor()

		#Extecute a command so that the cursor is active:
		stuff = cursor.execute('select * from {}'.format(table_name) )

		#Pull out column names:
		column_names = [ i[0] for i in cursor.description ]

		connection.close()
		
		return column_names

	
	
def Save_Dict_to_SQL(database_name,variable_dict,table_name,connect_to_server=True) :

	"""
	
	Save a row of values to an SQL database.
	
	Parameters
	-----------------
	
	database_name : str
	
	file path to the database
	
	variable_dict :  dict
	
	Dictionary containing keys corresponding to the relevant column headers
	and values correesponding to those to be saved in the table. 
	
	table_name : str
	
	Name of table in the database to save the dictionary to. 
	
	"""
	
	
	#Open the database:
	if connect_to_server==True : 
		connection, cursor = connect_to_maths_cluster_db()
	
	else : 
		#Open the connection to the databse:
		connection = sqlite3.connect(database_name)
		cursor = connection.cursor()
	
	params = list( variable_dict.values() )
	
	#Convert params to supported variable types:
	params = tuple([ convert_to_supported_type(i) for i in params ])
	
	
	#MySQL and sqlite use differetn syntax:
	
	if connect_to_server == True : 
		null_string = '%s'
		#Make string of ?'s of the right length for insertion into the table:
		for i in range(len(params)-1):
			null_string += ',%s'

		cursor.execute("insert into {} values({})".format(table_name,null_string) , params)
	
	else : 
		null_string = '?'
		for i in range(len(params)-1):
			null_string += ',?'
	
		cursor.execute("insert into {} values({})".format(table_name,null_string) , params)

	# Commit the changes to the database:
	connection.commit()

	# Is it necessary to close the connection every single time:
	connection.close()
	
	

	
def Pull_Value_From_Row(database_name,table_name,column_name,row_num,connect_to_server=True) : 

	"""
	
	Pull a specific property from a given row number.
	
	Parameters
	-----------------
	
	database_name : str
	
	file path to the database
	
	table_name : str
	
	name of the table witin the database
	
	column_name : str
	
	Name of the column of interest
	
	row_num : int
	
	row number. Note that pythonic convention starts
	list at entry 0 while SQL starts tables at row 1.
	
	Returns
	-------------
	
	value : variable contained in the specified row of
	a given column. 
	
	Usage example: (works if the .db file exists). 
	value = Pull_Value_From_Row("database.db","Results_Table","N",23)
	
	"""
	
	#open the connection to the SQL database:
	if connect_to_server==True : 
		connection, cursor = connect_to_maths_cluster_db()
	
	else : 
		#Open the connection to the databse:
		connection = sqlite3.connect(database_name)
		cursor = connection.cursor()
	
	
	
	#Pull out the specific value from the chosen row:
	cursor.execute("SELECT {} FROM {}".format(column_name,table_name))
	rows = cursor.fetchall()
	#pdb.set_trace()
	#value =  cursor.fetchall()[0][0] 
	value = rows[row_num-1][0]
	
	#If the value is unicode then we convert to a Python string:
	#This only matters in Python 2.X. In Python 3 unicode has been
	#renamed to string. 
	"""
	See: https://stackoverflow.com/questions/19877306/nameerror-global-name-unicode-is-not-defined-in-python-3
	"""
	
	if sys.version_info[0] < 3:
		if type(value) == unicode : 
			value = str(value) 
	
	#close the connection:
	connection.close()
	
	return value
	
def get_num_of_rows(database_name,table_name,connect_to_server=False) :
	
	"""
	
	Return the number of rows in a specific table in an SQL
	database.
	
	Parameters
	---------------
	
	database_name : str
	
	file path to the database
	
	table_name : str
	
	name of the table witin the database
	
	Returns
	------------
	
	num_rows : int
	
	Number of rows in the sql table.
	
	
	"""
	#open the connection to the SQL database:
	if connect_to_server==True : 
		
		connection, cursor = connect_to_maths_cluster_db()
	
	else : 
		#Open the connection to the databse:
		connection = sqlite3.connect(database_name)
		cursor = connection.cursor()

	cursor.execute("SELECT Count(*) FROM {}".format(table_name) ) 
	num_rows =  cursor.fetchall()[0][0]
	
	#close the connection:
	connection.close()
	
	return num_rows 
	
	
def Get_Param_Dicts_to_Sample(param_database,param_table_name,connect_to_server=False) :

	"""
	
	Function to read in the input parameters for simulations 
	from a table. 

	Parameters
	----------------

	parameter_database : str

	path to the database containing a table of parameters.


	Returns
	--------------

	parameter_dict_list : list of dict

	List of dictionaries containing the parameters to sample.
	
	Keys are the variable names and values are the corresponding
	values of the variable to sample. 

	e.g { 'N' : 1000 , 'r' : 0.1 , 'd' : 2 , 'boundary' : 'P' }
	
	"""
	
	#Empty arrays to store the output:
	parameter_dict_list = [ ] 

	num_of_rows = get_num_of_rows(param_database,param_table_name,connect_to_server=connect_to_server)
	
	for row_num in range(num_of_rows) : 
	
		#Empty dictionary:
		Input_Param_Dict = {}
		
		#Loop through the different input parameter names:
		for c_name in {'N','r','d','boundary','degree_scaling_parameter'} : 

			Input_Param_Dict[c_name] = Pull_Value_From_Row(param_database,param_table_name,c_name,row_num+1,connect_to_server=connect_to_server)
		
		parameter_dict_list.append(Input_Param_Dict) 
	
	return parameter_dict_list
