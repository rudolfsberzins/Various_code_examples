# input: 1) x: the independent variable, as a N dimensional vector as a numpy array
#        2) y: the dependent variable, as a N dimensional vector as a numpy array
#
# output: 1) the alpha parameter
#         2) the beta parameter
def univarlinreg(x,y):
    import numpy as np
    N = len(x)
    beta_top = np.mean(x)*np.mean(y) - np.dot(x,y)/N
    beta_bottom = np.mean(x)**2 - np.dot(x,x)/N
    beta = beta_top/beta_bottom

    alpha = np.mean(y)-beta*np.mean(x)

    return alpha, beta