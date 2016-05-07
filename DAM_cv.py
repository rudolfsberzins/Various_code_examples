from __future__ import division
# input:  1) training data in the form of a N by d numpy array, where N is the number of training data points and d is the number of dimensions
#         2) training labels in the form of a N by 1 numpy vector, where N is the number of training data points
#         3) a random permutation of entries as a numpy array, e.g. np.random.permutation(len(trainlabels))
# output: 1) the optimal k
#         2) an error matrix (numpy array) of size (5,13) where column i consists of the accuracy for the 5 folds for k=i
#
# The random-permuted vector rand_perm should be used for generating 5 folds, where the first fold consists of the first N/5 elements from rand_perm, rounded up to the nearest integer; the second fold consists of the next N/5 elements, etc, and the fifth fold consists of the remaining elements
# note: to create the folds consider: KFold(len(trainlabels), n_folds=5) from sklearn.cross_validation (http://scikit-learn.org/stable/modules/generated/sklearn.cross_validation.KFold.html)
# note: once you have the folds use the rand_perm vector to get the random indices in the training data and labels
def cv(train, trainlabels, rand_perm):
    import numpy as np
    import math
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

    trainnum, dim = train.shape
    fold_size = math.ceil(trainnum/5)
    errors = np.zeros((5,13))
    k_to_test = range(1,27,2)

    for i in range(5):
        if i < 4:
            val_fold = rand_perm[range(int(i*fold_size), int((i+1)*fold_size))]
        else:
            val_fold = rand_perm[range(int(i*fold_size), int(trainnum))]

        val_data = train[val_fold]
        val_l = trainlabels[val_fold]
        rest_folds = np.setdiff1d(rand_perm, val_fold)
        rest_data = train[rest_folds]
        rest_l = trainlabels[rest_folds]

        a = []

        for k in k_to_test:
            dist, pred = knn(rest_data, val_data, rest_l, k)
            b = len(pred) - np.sum(np.equal(pred,val_l))
            a.append(b)

        errors[i,:] = a


    mean_errors = np.mean(errors,0)
    op_k_ind = np.argmin(mean_errors)
    optimal_k = k_to_test[op_k_ind]

    return optimal_k, errors


