
# input: datamatrix as loaded by numpy.loadtxt('shapes.txt')
# output:  1) the eigenvalues in a vector (numpy array) in descending order
#          2) the unit eigenvectors in a matrix (numpy array) with each column being an eigenvector (in the same order as its associated eigenvalue)
#
# note: make sure the order of the eigenvalues is descending, and the eigenvectors have the same order as their associated eigenvalues
def pca(data):
    import numpy as np
    Sigma = np.cov(data)
    evals, evecs = np.linalg.eig(Sigma)

    return evals, evecs