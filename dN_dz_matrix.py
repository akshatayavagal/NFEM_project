import numpy as np 
def dN_dz_matrix():
    dN_dZ1 = -1/2
    dN_dZ2 =  1/2
    dN_dZ = np.array([[dN_dZ1],[dN_dZ2]])
    return dN_dZ

print(dN_dz_matrix())