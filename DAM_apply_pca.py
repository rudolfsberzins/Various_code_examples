# input: Datamatrix as loaded by numpy.loadtxt('Irisdata.txt')
# output: Datamatrix of the projected data onto the two first principal components.

def apply_pca(data):
    from sklearn.preprocessing import StandardScaler
    X_std = StandardScaler().fit_transform(data)

    from sklearn.decomposition import PCA as sklearnPCA
    sklearn_pca = sklearnPCA(n_components=2)
    Y_sklearn = sklearn_pca.fit_transform(X_std)

    return Y_sklearn
