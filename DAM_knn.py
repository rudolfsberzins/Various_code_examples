# input: train: 1) train data (without labels) in the form of a N by d numpy array, where N is the number of training data points and d is the number of dimensions
#               2) test: test data (without labels) in the form of a M by d numpy array, where M is the number of test data points and d is the number of dimensions
#               3) trainlabels: labels for training data in the form of a N by 1 numpy vector, where N is the number of training data points
#               4) k: paramater k
# output:1) distance matrix (numpy array) between test and training samples 
#        2) vector (numpy array) consisting of the predicted classes for the test data
#
# note: the labels should **not** be part of the train/test data matrices!
def knn(train, test, trainlabels, k):
    from collections import Counter
    import numpy as np
    def dot(v, w):
        """v_1 * w_1 + ... + v_n * w_n"""
        return sum(v_i * w_i
            for v_i, w_i in zip(v, w))

    def sum_of_squares(v):
        """v_1 * v_1 + ... + v_n * v_n"""
        return dot(v, v)

    def vector_subtract(v, w):
        """subtracts corresponding elements"""
        return [v_i - w_i
                for v_i, w_i in zip(v, w)]

    def squared_distance(v, w):
        """(v_1 - w_1) ** 2 + ... + (v_n - w_n) ** 2"""
        return sum_of_squares(vector_subtract(v, w))

    def distance(v, w):
        return np.sqrt(squared_distance(v, w))

    dis_mat = []

    for i in range(len(test)):
        temp_array = []
        for j in range(len(train)):
            temp_array.append(distance(train[j], test[i]))
        dis_mat.append(temp_array)


    dis_mat = np.array(dis_mat)

    labels_after_test = []

    for i in range(len(dis_mat)):
        # get sorted indexes
        inds = dis_mat[i].argsort()[0:k]

        # get labels
        label = []
        for j in inds:
            label.append(trainlabels[j])
        c = Counter(label)
        mc = c.most_common(1)[0][0]

        labels_after_test.append(mc)

    return dis_mat, np.array(labels_after_test)


