import numpy as np


def nn_search(X=None, x=None, y=None, z=None, 
                  n=1,
                  return_distances=False):
    """
    Vectorized nearest neighbor search accepting either
    an (N,D) matrix (first kwarg) or three (N,) vectors
    (representing x, y, and z coordinates) to identify
    the nearest neighbors to a selection of points. Uses
    Euclidean distance.

    This is naiively vectorized, so will operate quickly
    but will likely blow up for very large N.

    Core algorithm taken from the
    vectorized_nearest_neighbor.py routine in 
    astroML (Statistics, Data Mining, and Machine Learning
    in Astronomy; Ch. 2.5, pg 55).

    Parameters
    -----------
    X   :   (N,D) numpy array, optional
         (N,D) numpy array to perform the NN search.
         This must be provided unless at least two of the
         below x,y,z are provided.

    x,y,z : (N,) numpy arrays, optional
        If X is not provided, can provided 3 1D arrays
        which will be used to construct X. Must provide
        at least two of these if used

    n : integer, optional
        Order of nearest neighbor to search for. By default
        this finds the nearest neighbor (n=1), but can 
        easily return second-nearest (n=2), and so on.
        Default : 1.

    return_distances : bool, optional 
        By default, returns just the indexes corresponding
        to the nearest neighbor for each value. If True,
        also returns the distances to these values.
    """

    if X is None:
        none_count = int(x is None) + int(y is None) + int(z is None)
        if none_count > 1:
            print("If X is not provided, must provide at least 2 1D arrays")
            raise ValueError

        elif none_count == 0:
            l = [x,y,z]
        elif none_count == 1:
            if (x is None):
                l = [y,z]
            elif (y is None):
                l = [x,z]
            else:
                l = [x,y]
         

        X = np.array(l).T

    # now do the sort
    XXT = np.dot(X, X.T)
    Xii = XXT.diagonal()
    D   = Xii - 2 * XXT + Xii[:, np.newaxis]

    print(np.shape(D))
    #
    #
    #
    nn_sort = np.argsort(D, axis=1)[:,n]

    if return_distances:
        distances = D[nn_sort]
        return nn_sort, D[np.arange(len(nn_sort)), nn_sort]
    else:
        return nn_sort
