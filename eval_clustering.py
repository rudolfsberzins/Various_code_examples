import numpy as np
from kmeans import kmeans


# input: 1) Datamatrix as loaded by numpy.loadtxt()
#        2) Random datamatrix of same size as input 1
#        3) numpy array of length 10 with initial center indices
# output: vector (numpy array) conisting of the gap statistics for k=1..10
# note this function should be called 10 times and averaged for the report.
def eval_clustering(data, randomData, initialCenters):
    step = len(initialCenters)
    results = []
    for i in range(step):
        reg_c, reg_labs = kmeans(data, initialCenters[: i + 1])
        rand_c, rand_labs = kmeans(randomData, initialCenters[: i + 1])
        e = objective_function_kmeans(data, reg_c, reg_labs)
        e_rand = objective_function_kmeans(data, rand_c, rand_labs)
        results.append(gap_statistics(e_rand, e))

    results = np.array(results)
    return results


# input:  	1) datamatrix as loaded by numpy.loadtxt()
#			2) numpy array of the cluster centers.
#			3) vector (numpy array) consisting of the assigned cluster for each observation.
# output:	the k-means objective function value as specified in the assignment for the given k
# note k is infered based on number of elements in centers
def objective_function_kmeans(data, centers, clusterLabels):
    k = len(centers)
    sum = 0
    for i in range(k):
        indexes = [index for index, group in enumerate(clusterLabels) if group == i]
        for j in indexes:
            sum += (np.linalg.norm(data[j] - centers[i])) ** 2

    return sum


# input:  	1) vector (numpy array) of objective function values for k=1..10 for a random dataset
#			2) vector (numpy array) of objective function values for k=1..10 for a given dataset
# output:	vector (numpy array) of the computed gap statistics for each k=1..10
# note should calculate step 4 in 4.a in assignment description
def gap_statistics(Erand, E):
    a = np.log(Erand) - np.log(E)
    return a
