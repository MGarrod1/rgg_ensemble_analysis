"""
This code includes functions for computing 
the eigenvalues of large sparse matrices.

The main eigenvalue of interest is the second
smallest eigenvalue of the graph laplacian
(known as the algerbaic connectivity).

"""



#Import packages:
from scipy.sparse.linalg import LinearOperator, eigsh, minres
import numpy as np

def Largest_EV(M,return_eigenvector = False) :

	"""
	
	Computes the largest eigenvalue of a sparse matrix M
	
	Parameters
	------------
	
	M : scipy sparse matrix
	
	Returns
	----------
	
	EV : double 
	
	The largest eigenvalue of the matrix M
	
	
	"""
	M = M.asfptype()

	if return_eigenvector == True:
		evals_large , evecs_large = eigsh(M, 1, which='LM', tol=1E-6, return_eigenvectors=True)
		return evals_large[0], np.transpose(evecs_large)[-1]
	else:
		evals_large = eigsh(M, 1, which='LM', tol=1E-6, return_eigenvectors=False)
		return evals_large[0]
	


def Smallest_Non_Zero(M, return_eigenvector = False , v0 = None , Threshold = 1E-8 , zero_eig_limit = 100) :

		"""
		
		Method to find the smallest non zero eigenvalue of a positive semidefinite matrix.

		-This method can be used when we want to estimate the mu2 value
		without having to assume that we have initially extracted the LCC.
		
		-Recommend using this even for matrices which we expect to have a single
		zero eigenvalue as the eigsh method can sometimes return multiple values
		close to zero due to numerical error in the eigensolver.

		Terminates if the first 100 eigenvalues are found to be below
		the threshold.

		From M we can construct a linear operator which is used in the eigenvalue computation:
		(See https://github.com/scipy/scipy/issues/4170 or https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.sparse.linalg.LinearOperator.html for details )


		Parameters
		---------------

		M : scipy sparse matrix

		matrix to find eigenvalues of (This code has been desigend for this to be a
		graph laplacian matrix).

		return_eigenvector : bool (optional)

		option to return the eigenvector associated with the smallest non-zero eigenvalue.


		v0 : vector/array (optional)
		
		Vector to use as a starting point for iteration.

		Threshold : float (optional)

		Threshold for eigenvalues to be assumed to be zero.
		Set to 1E-8 by default. Has to be non-zero otherwise
		we will incorrectly return zero eigenvalues as non-zero
		due to numerical error.

		zero_eig_limit : int

		Terminate if we find more than 'zero_eig_limit' eigenvalues close to zero.

		"""
		
		Found_Non_Zero = False
		k = 0
		
		#Construct the linear operator to use in the solver:
		OPinv = LinearOperator(matvec=lambda v: minres(M, v, tol=1e-6)[0], shape=M.shape, dtype=M.dtype)

		

		"""Compute eigenvalues and eigenvectors until we find one
		which is greater than a certain threshold."""
		while Found_Non_Zero == False and k < zero_eig_limit :

			#Compute eigenvalues and eigenvectors: (looks for at least two eigenvalues)
			evals_small , evecs_small = eigsh(M, k+2, sigma= 0, which='LM', tol = 1E-6, return_eigenvectors = True , OPinv = OPinv , v0 = v0)

			#Check if any eigenvalues are above the threshold:
			if max(evals_small) > Threshold :
				Found_Non_Zero = True
				Mu2 = max(evals_small)
			k = k + 1

		if return_eigenvector == True :
			return Mu2 , np.transpose(evecs_small)[-1]
		else: 
			return Mu2
