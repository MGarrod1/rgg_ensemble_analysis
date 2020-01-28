
This package provides tools for generating and analysing ensembles of Random Geometric Graphs (RGGs). 

The code is used to perform the numerical simulations in:

Garrod, Matthew, and Nick S. Jones. "Large algebraic connectivity fluctuations in spatial network ensembles imply a predictive advantage from node location information." Physical Review E 98.5 (2018): 052316. [https://doi.org/10.1103/PhysRevE.98.052316](https://doi.org/10.1103/PhysRevE.98.052316)

# Installation:

(sudo) python setup.py install


## Construction of RGGs
RGGs are constructed using a KD-TREE based method which allows efficient computation of the set of points within connection radii for a particular node. Once an RGG is generated it is possible to generate a `network class' object. The network class has methods for computing properties of the RGG such as the mean degree and average shortest path length. Many of the functions are based on those from networkx, however, the package also includes computation of eigenvalues using scipy's sparse linear algebra library.

## Analysis of ensembles

Once we have sampled multiple RGGs we can analyse the statistical properties by using the `ensemble class.' The ensemble class and accompanying functions allow computation of descriptive statistics for particular ensembles such as the mean and standard deviation of various networkx metrics.

