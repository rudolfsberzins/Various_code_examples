
# input: 1) train: train data in the form of a N by d numpy array, where N is the number of training data points and d is the number of dimensions
#        2) test: test data in the form of a M by d numpy array, where M is the number of test data points and d is the number of dimensions    
# output: 1) centered and normalized train data as a numpy array
#         2) centered and normalized test data as a numpy array 
def cent_and_norm(train, test):
    import numpy as np
    def scale(data_matrix):
        """returns the means and standard deviations of each column"""
        num_rows, num_cols = np.shape(data_matrix)
        means = [np.mean(data_matrix[:,j])
                 for j in range(num_cols)]
        stdevs = [np.std(data_matrix[:,j])
                  for j in range(num_cols)]
        return means, stdevs

    def rescale(data_matrix):
        """rescales the input data so that each column
        has mean 0 and standard deviation 1
        leaves alone columns with no deviation"""
        means, stdevs = scale(data_matrix)

        def rescaled(i, j):
            if stdevs[j] > 0:
                return (data_matrix[i][j] - means[j]) / stdevs[j]
            else:
                return data_matrix[i][j]

        num_rows, num_cols = np.shape(data_matrix)
        return make_matrix(num_rows, num_cols, rescaled)

    def make_matrix(num_rows, num_cols, entry_fn):
        """returns a num_rows x num_cols matrix
        whose (i,j)th entry is entry_fn(i, j)"""
        return [[entry_fn(i, j) # given i, create a list
                 for j in range(num_cols)] # [entry_fn(i, 0), ... ]
                for i in range(num_rows)] # create one list for each i

    norm_train = np.array(rescale(train))

    def rescale_test(test):

        means, stdevs = scale(train)
        row, col = np.shape(test)
        def rescaled_test(i, j):
            if stdevs[j] > 0:
                return (test[i][j] - means[j]) / stdevs[j]
            else:
                return test[i][j]
        return make_matrix(row, col, rescaled_test)

    norm_test = np.array(rescale_test(test))

    return norm_train, norm_test

