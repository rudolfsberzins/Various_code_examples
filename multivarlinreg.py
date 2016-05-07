# input: 1) x: the independent variables (data matrix), as a N x M dimensional matrix as a numpy array
#        2) y: the dependent variable, as a N dimensional vector as a numpy array
#
# output: 1) the regression coefficients as a (M+1) dimensional vector as a numpy array
#
# note: the regression coefficients should include the w_0 (the free parameter), thus having dimension (M+1).
# note: The tested datamatrix is **NOT** extended with a column of 1's - if you prefer to do this, then you can do it inside the function by extending x.       
def multivarlinreg(x, y):
    import numpy as np
    N,_ = x.shape
    onevec = np.ones((N,1))
    X = np.concatenate((onevec, x), axis = 1)

    w = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(X), X)), np.transpose(X)), y)

    return w