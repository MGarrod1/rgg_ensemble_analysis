"""

Example code for reproducing some of the results
from the paper.

For a particular RGG ensemble. Picks a certain point in parameter
space. We then draw some samples from the algebraic connectivity
distriution and 
1) display the mean and CV.
2) Plot a histogram (like that from figure 4...). 


Allows the user to check different datapoints
from the plots in the paper.


"""

#Import statements:
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import math

import rgg_ensemble_analysis

import rgg_ensemble_analysis.Analytic_Functions as Analytic
import rgg_ensemble_analysis.RGG as RGG
import rgg_ensemble_analysis.Statistics_Computation as StatComp
import rgg_ensemble_analysis.get_connection_radii_estimate as connection


#Choose parameters for the ensemble:
Samples = 50
N = 1000
d = 10
boundary = 'P'

#Expected mean degree parameter
#The code below estimates the radius value
#required in order to obtain this specific degree value.
#Kappa = 50.0
C = 1.5
#Kappa = 13.8155
Kappa = C*math.log(N)

#Create a dictionary for the input.
Input_Param_Dict = {}
Input_Param_Dict['N'] = N
Input_Param_Dict['d'] = d
Input_Param_Dict['boundary'] = boundary



#Should in fact derive radius from the scaling param:
Input_Param_Dict['degree_scaling_parameter'] = Kappa

#Number of samples for estimating the connection radius:
radii_samps = 10

#Set the radius:
#Perform several radius estimates and take the mean:
radii_estimates = [ connection.Get_Required_Kappa(Kappa , N, 1000 , d , Boundaries = boundary , positions = None , pattern="uniform") for k in range(radii_samps) ]
r = np.mean(radii_estimates)
Input_Param_Dict['r'] = r


database_name = "example_database.db"
graph_data_table_name = "graph_properties"

#Choose the properties to sample:
pattern = "uniform"
methods_to_call = [ "Algebraic_Connectivity" ,'LCC_Min_Degree' , 'Mean_Degree' ] 

### Sample RGGs ###

RGG.Sample_RGGs(N, r, d, boundary,Kappa,Samples, methods_to_call, database_name,graph_data_table_name,pattern,Save_Positions=False,connect_to_server=False)

### DATA READING IN ###
ensemble = RGG.Graph_Ensemble(database_name,graph_data_table_name,Input_Param_Dict)

#Check the number of nodes in the ensemble:
#print("N={}".format(ensemble.N) )

mu2_array = ensemble.get_prop_array('algebraic_connectivity')
kmin_array = ensemble.get_prop_array('LCC_Min_Degree')


Mean_Mu2 = np.mean(mu2_array)
CV_Mu2 = StatComp.CV(mu2_array)
CV_ER_Mu2 = StatComp.CV_Error(mu2_array)


#Compute the coefficient of variation and its error:
print("Number of samples = {}".format(len(mu2_array)) ) 

print("Mean Mu2 = {}".format(Mean_Mu2) ) 
if boundary =='P' : 
	print("Theoretical Mean Mu2 (periodic) = {}".format( Analytic.Theory_Algebraic(N,Kappa,d) ) ) 

print("CV = {} +/- {}".format(CV_Mu2 ,CV_ER_Mu2  ) ) 

#Also report the spearman correlation between mu2 and kmin.
#Note that we may need many samples to obtain a confident
#estimate of this quantity.
#correlation_mu2_kmin = StatComp.Get_Spearman_Correlation(mu2_array, kmin_array ) 
correlation_mu2_kmin, pval = stats.spearmanr(mu2_array, kmin_array)
print("correlation = {} (P val = {})\n\n".format(correlation_mu2_kmin,pval) ) 


###Plot the histogram ####
file_path = "hist_example"

d = ensemble.d


#Code below bins mu2 values according to the corresponding
#values of the minimim degree (kmin). 

#Zip mu2 and kmin together:
combined_array = zip(mu2_array,kmin_array)

#get the set of different kmin values:
distinct_kmin_vals = set([ int(k) for k in set(kmin_array) ] )

#Make a dictionary where keys are the values:
kmin_mu2_dict = { }
for k in distinct_kmin_vals :
	kmin_mu2_dict[k] = [ ]

for i in range(len(mu2_array)) :
	kmin_mu2_dict[kmin_array[i]].append(mu2_array[i])
	
#Make histogram bins:
num_of_bins = 20
BIN_POINTS = np.linspace(0.9*min(mu2_array), 1.1*max(mu2_array), num = num_of_bins)

plt.clf()
plt.figure(1)
ax = plt.gca()

#colors picked from: https://www.w3schools.com/colors/colors_picker.asp?color=80ced6
colors = ['#ff9933' , '#66ff66' , '#0099ff'  , '#ff66ff', 'r' , 'b' , 'g' , 'm' ,  'y' , 'k' , 'y']

colors = [ 'k' , 'b' , 'r' , 'g' , 'm']
hatching_patterns = [ '/', '\\', '|', '-' ,'+', 'x', 'o', 'O', '.', '*' ]
line_styles = [ '-', '--', '-.', ':']

col_index = 0
all_counts = [ ]
for k in kmin_mu2_dict.keys() :

	#n, bins, hist = plt.hist(kmin_mu2_dict[k], bins=30, alpha=0.6,label="mindeg = {}".format(k))
	counts , bins = np.histogram(kmin_mu2_dict[k] , bins= BIN_POINTS )

	#bins is len(counts)+1 so need to take half points:
	mid_points = [ (bins[i+1]+bins[i])/2.0 for i in range(len(bins)-1)  ]

	#Must do my own normalisation over the total counts so that the normed version is P(mu2 for all degrees.
	#normalize by dividing by the total # of counts:
	counts = [ i/float(len(mu2_array)) for i in counts ]

	#Use this if we can actually get transparncey working?
	#ax.fill_between(mid_points, counts, alpha=.4, label="$\kappa_{min}$ = " + str(k),color=colors[col_index],hatch=hatching_patterns[col_index])
	#ax.plot(mid_points, counts, c=colors[col_index])

	#Alternative if we only want different styles:
	ax.plot(mid_points, counts, c = colors[col_index] , linestyle= line_styles[col_index] , linewidth=2,label = "$\kappa_{min}$ = " + str(k))

	col_index +=1


	all_counts = np.concatenate((all_counts,counts))

"""
old version of the loop. 


for k in kmin_mu2_dict.keys() :
	#n, bins, hist = plt.hist(kmin_mu2_dict[k], bins=30, alpha=0.6,label="mindeg = {}".format(k))
	counts , bins = np.histogram(kmin_mu2_dict[k] , bins= BIN_POINTS )

	#bins is len(counts)+1 so need to take the mid points:
	mid_points = [ (bins[i+1]+bins[i])/2.0 for i in range(len(bins)-1)  ]

	#Normalize histograms by the total counts:
	counts = [ i/float(len(mu2_array)) for i in counts ]

	ax.fill_between(mid_points, counts, alpha=.4, zorder=5,label="$\kappa_{min}$ = " + str(k),color=colors[col_index],hatch=hatching_patterns[col_index])
	ax.plot(mid_points, counts, c=colors[col_index] )

	col_index +=1
	all_counts = np.concatenate((all_counts,counts))
"""

plt.title("$E(\mu_2)$ = {} , $CV(\mu_2)$ = {} +\- {}".format( round(Mean_Mu2,3), round(CV_Mu2,3),  round(CV_ER_Mu2,3) ) ) 
plt.legend(loc=2)
plt.xlabel("$\mu_2$", fontsize=20)
plt.ylabel("$P(\mu_2)$", fontsize=20)
plt.ylim(0.0, 1.1*max(all_counts))
plt.savefig(file_path,bbox_inches="tight",format="pdf")
