# input: 1) x: the independent variable, as a N dimensional vector as a numpy array
#        2) y: the dependent variable, as a N dimensional vector as a numpy array
#        3) alpha: the alpha parameter
#        4) beta: the beta parameter
#
# output: 1) the root mean square error (rmse) 

def rmse(x, y, alpha, beta):
    import numpy as np
    t = np.add(alpha, beta*x)
    rmse = np.sqrt(np.mean((y-t)**2))
    return rmse
