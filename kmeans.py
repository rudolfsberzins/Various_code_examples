import scipy.spatial
import numpy as np


# input:  	1) datamatrix as loaded by numpy.loadtxt('Irisdata.txt').
#			2) vector (numpy array) of indices for k seed observation.
# output: 	1) numpy array of the resulting cluster centers. 
#			2) vector (numpy array) consisting of the assigned cluster for each observation.
# note the indices are an index to an observation in the datamatrix
def kmeans(data, seedIndices):
    centroids = []

    for i in seedIndices:
        centroids.append(data[i])

    def cluster_centroids(data, clusters, k=None):
        if k is None:
            k = np.max(clusters) + 1
        result = np.empty(shape=(k,) + data.shape[1:])
        for i in range(k):
            np.mean(data[clusters == i], axis=0, out=result[i])
        return result

    k = len(centroids)

    for _ in range(100):
        # Squared distances between each point and each centroid.
        sqdists = scipy.spatial.distance.cdist(centroids, data, 'sqeuclidean')
        # Index of the closest centroid to each data point.
        clusters = np.argmin(sqdists, axis=0)

        new_centroids = cluster_centroids(data, clusters, k)
        if np.array_equal(new_centroids, centroids):
            break

        centroids = new_centroids

    return centroids, clusters
