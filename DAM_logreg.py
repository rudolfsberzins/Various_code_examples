from __future__ import division
# input:  1) train data (without labels) in the form of a N by d numpy array, where N is the number of training data points and d is the number of dimensions
#         2) trainlabels: labels for training data in the form of a N by 1 numpy vector, where N is the number of training data points
#         3) test: test data (without labels) in the form of a M by d numpy array, where M is the number of test data points and d is the number of dimensions
#
# output: 1) vector (numpy array) consisting of the predicted classes for the test data
#         2) the beta parameter of the model
# note: the labels should **not** be part of the train/test data matrices!
def logreg(train_data, train_labels, test_data):
    import sklearn.linear_model as lm
    import numpy as np
    rows, _ = train_data.shape
    ov = np.ones((rows,1))
    train_data_ones = np.concatenate((ov, train_data), axis=1)

    rows, _ = test_data.shape
    ov = np.ones((rows,1))
    test_data_ones = np.concatenate((ov, test_data), axis=1)

    model = lm.LogisticRegression()
    model_train = model.fit(train_data_ones, train_labels)

    predictions = model_train.predict(test_data_ones)
    w = model.coef_[0]
    return predictions, w
