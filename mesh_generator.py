#  file_name = mesh_generator.py
#  Generate list of position of nodes according to a geometric series
#     for assignement in "Nonlinear Finite Element Methods" 
#     in summer term 2019
#     lecturer in charge: Dr. Geralf H�tter
import numpy as np 
import matplotlib.pyplot as plt
from Parameters import*

def mesh_generator(nelem):
    # Input Parameters
    #ro=75e-6                      #outer radius
    #ri=15e-6                      #radius of inclusion
    #nelem=10                      #number of elements
    meshrefinementfactor=5        #ratio of element sizes at outer and inner radius
    # ratio between element sizes of subsequent elements for a geometric series
    q = meshrefinementfactor**(1/(nelem-1))
    # size of first element
    lelem=(ro-ri) * (1-q) / (1-meshrefinementfactor*q)
    rnode = ri
    rnodes = np.array([ri])

    # loop over all elements
    for i in range(nelem):
        rnode += lelem
        rnodes = np.append(rnodes,rnode)
        lelem = lelem * q 

    # visualize location of nodes 
    #plt.plot(rnodes,np.zeros([nelem+1]),'*')
    #plt.show()
    return rnodes
