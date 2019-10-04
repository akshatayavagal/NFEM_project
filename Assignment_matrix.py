import numpy as np
from Parameters import *

def Assignment_matrix(nelem):
    
    Assignment_matrixs = []
    for i in range(nelem):
        a = np.zeros((2,nelem+1))
        a[0][i]=1
        a[1][i+1]=1
        Assignment_matrixs.append(a)
    return Assignment_matrixs

