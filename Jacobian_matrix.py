import numpy as np
from dN_dz_matrix import *

def Jacobian_matrix(x_ele,dN_dZ):
    # Jacobian_matrix  
    dN_dZ = dN_dz_matrix()
    J = np.matmul(np.transpose(dN_dZ),x_ele)
    return J


