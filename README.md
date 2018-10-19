
This package provides tools for generating and analysing ensembles of Random Geometric Graphs (RGGs). 

# Installation:

(sudo) python setup.py install


## Construction of RGGs
RGGs are constructed using a KD-TREE based method which allows efficient computation of the set of points within connection radii for a particular node. Once an RGG is generated it is possible to generate a `network class' object. The network class has methods for computing properties of the RGG such as the mean degree and average shortest path length. Many of the functions are based on those from networkx, however, the package also includes computation of eigenvalues using scipy's sparse linear algebra library.

## Analysis of ensembles

Once we have sampled multiple RGGs we can analyse the statistical properties by using the `ensemble class.' The ensemble class and accompanying functions allow computation of descriptive statistics for particular ensembles such as the mean and standard deviation of various networkx metrics.

