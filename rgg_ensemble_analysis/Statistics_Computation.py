"""

Functions to compute various statistics of the data
and their associated errors.

"""


#Compute statistics about the variance of an estimator:

import numpy as np 
from scipy import stats
import math

#The following functions as involved in estimating the standard 

def Moment(Sample,k) :

	"""
	This function computes the kth moment of the
	Sample.
	
	"""
	Mean = np.mean(Sample)
	
	Moment = 0.0
	for i in range(0,len(Sample) ) :
		Moment = Moment + (Mean - Sample[i] )**k 
		
	return Moment/(len(Sample))
	
def D4(Sample) :
	"""Compute the fourth central moment of a sample
	See: https://en.wikipedia.org/wiki/Central_moment
	for defintion"""
	M = float( len(Sample) )
	D4 = ((M-1)/(M**3))*( (M**2 - 3*M + 3)*Moment(Sample,4) + 3*(2*M -3)*Moment(Sample,2)**2 )
	
	return D4
	
def VarSE(Sample) :


	"""
	Returns the standard error on the variance of the sample given.
	
	Formula taken from: 
	Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels. "Previous estimates of mitochondrial DNA mutation level variance 	did not account for sampling error: comparing the mtDNA genetic bottleneck in mice and humans." The American Journal of Human 	Genetics 86.4 (2010): 540-550.
	
	Parameters
	-----------
	
	Sample : list
	
	list of values
	
	Returns
	-----------
	
	SE : float
	
	Standard error on the variance of a sample
	
	
	"""


	M = float( len(Sample) )
	
	SE = ( (1.0/M)*(D4(Sample) - ( (M-3)/(M-1) )*np.var(Sample)**2  )     )**0.5

	return SE
	
def CV(Sample) : 

	"""
	Computes the coefficient of variation of the values.
	
	Arguments
	------------
	
	Sample : list
	
	
	"""
	return np.std(Sample)/np.mean(Sample)
	
def CV_Error(Sample) :
	"""
	Compute the standard error on the coefficient of variation.
	"""
	try :
		part1 = (VarSE(Sample)/(  2.0*(np.var(Sample)**0.5)*np.mean(Sample) ) )**2
		part2 = ( (((np.var(Sample))**0.5 )*stats.sem(Sample)/(np.mean(Sample)**2)   )  )**2
		Coeff_ERR = ( part1   +  part2    )**0.5
	except :
		print("Error encountered")
		print("Sample length = " + str(len(Sample)))
		print("Setting error to zero by defulat")
		Coeff_ERR = 0.0


	return Coeff_ERR
	

	
#Bootsrapping Tools:

#Given a sample we can estimate the standard error in the variance
def Bootstrap_SE(Sample,Redraws) :
	Data_Points = len(Sample)
	Var = np.var(Sample)

	Vars = [ ]
	
	for j in range(0,Redraws) : 
		Sample2 = [  ]
		for i in range(0,(len(Sample))) :
			#Draw a random integar less than the sample size
			Entry = int(math.floor( np.random.uniform(0, (len(Sample)), 1) ) )
			Sample2.append(Sample[Entry])      
		New_Var = np.var(Sample2)
		Vars.append(New_Var)
	
	  
	
	Var_SE = stats.sem(Vars)
		
	return (Redraws**0.5)*Var_SE
	
def Bootstrap_Mean_SE(Sample,Redraws) :
	Data_Points = len(Sample)
	
	Mean = np.mean(Sample)

	Means = [ ]
	
	for j in range(0,Redraws) : 
		Sample2 = [  ]
		
		for i in range(0,(len(Sample))) :
			#Draw a random integar less than the sample size
			Entry = int(math.floor( np.random.uniform(0, (len(Sample)), 1) ) )
			Sample2.append(Sample[Entry])  
				
		New_Mean = np.mean(Sample2)
		Means.append(New_Mean)
	
	  
	
	Mean_SE = stats.sem(Means)
	Boot_Mean = np.mean(Means)
		
	return (Redraws**0.5)*Mean_SE
	


def Get_Spearman_Correlation(Array_1, Array_2, P_Val_Threshold=(10 ** (-4))):
	
	"""
	Returns the correlation between two arrays.

	If p val > 10^-4 we simply report the correlation as zero
	
	Parameters
	--------------
	
	Array_1 : list
	
	
	Array_2 : list
	
	P_Val_Threshold :float
	
	Set a threshold for the p-value. If the p value
	if greater than this threshold then we can set
	the correlation to zero (ie. we are not confident
	that the correlation exists). 

	"""
	
	P_Val_Threshold = 1.1
	Correlation, pval = stats.spearmanr(Array_1, Array_2)

	if pval < P_Val_Threshold:
		return Correlation

	else:
		return 0.0


def Bootstrap_Correlation_Confid_Int(Sample_1, Sample_2, Redraws=50):
	"""

	Use the bootstrapping method to estimate the confidence interval on the
	Spearman Correlation

	Returns the 95% confidence interval by default

	(Does p-value threshold want to be an additional argument?)

	Parameters
	--------------

	Sample_1 : list

	First sample

	Sample_2 : list

	Second sample

	Redraws : int

	Number of times to same with replacement from the joint distribution
	of sample_1 and sample_2


	Returns
	-----------

	95% confidence interval on the Spearman correlation.

	"""

	Data_Points = len(Sample_1)

	Original_Correlation = Get_Spearman_Correlation(Sample_1, Sample_2, P_Val_Threshold=20.0)

	Correlations = []
	Differences = []

	for j in range(0, Redraws):
		Sample2 = []

		Redraw_Sample_1 = []
		Redraw_Sample_2 = []

		# Redraw the samples:
		for i in range(0, (len(Sample_1))):
			#Redraw pairs of values:
			Entry = int(math.floor(np.random.uniform(0, (len(Sample_1)), 1)))
			Redraw_Sample_1.append(Sample_1[Entry])
			Redraw_Sample_2.append(Sample_2[Entry])

		Redrawn_Correlation = Get_Spearman_Correlation(Redraw_Sample_1, Redraw_Sample_2, P_Val_Threshold=20.0)
		Correlations.append(Redrawn_Correlation)
		Differences.append(Redrawn_Correlation - Original_Correlation)

	# sort the list of differences:
	Sorted_Differences = np.sort(Differences)
	
	return Sorted_Differences[int(math.ceil(0.95 * len(Differences)))]


