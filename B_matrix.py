import numpy as np 
from shape_matrix import *
from dN_dz_matrix import *

def B_matrix(z,R,dN_dZ,J): 
    N = shape_matrix(z)
    dN_dZ = dN_dz_matrix()
    dN_dR=(1/J)*dN_dZ
    B11 = np.asscalar(dN_dR[0])
    B12 = np.asscalar(dN_dR[1])
    B21 = np.asscalar(N[0]/(np.matmul(np.transpose(N),R)))
    B22 = np.asscalar(N[1]/(np.matmul(np.transpose(N),R)))
    B31 = np.asscalar(N[0]/(np.matmul(np.transpose(N),R)))
    B32 = np.asscalar(N[1]/(np.matmul(np.transpose(N),R)))
    B = [[B11,B12],[B21,B22],[B31,B32]]
    return B