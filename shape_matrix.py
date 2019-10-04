import numpy as np 

def shape_matrix(z):
    N1 = 1/2 * (1-z)
    N2 = 1/2 * (1+z)
    N = np.array([[N1],[N2]])
    return N

#print(shape_matrix(0))