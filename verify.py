'''verification'''
import numpy as np 
from Parameters import *
from mesh_generator import *
from Main_program import *
from material_routine import *
U_s = []
sigma_rr=[]
for j in range(nelem+1):
    r=rnodes[j]
    U = ((ri**3)*volumetric_strain)/(3*r*r)
    sigma_verif=(-2*youngs_modulus_E*volumetric_strain*ri**3)/(3*(1+poission_ratio_v)*r**3)
    U_s.append(U)
    sigma_rr.append(sigma_verif)


fig, ax=plt.subplots(figsize=(10, 10))
plt.subplot(2,1,1)
plt.plot(U_s,rnodes,'-o',markevery=[0])
plt.plot(U_global,rnodes,'-r',markevery=[0])
plt.title('Displacement_diffrent nodes_elastic_area')
plt.xlabel('Displacement[m]')
plt.ylabel('nodes_length[m]')

plt.subplot(2,1,2)
plt.plot(,rnodes,'-o',markevery=[0])
plt.plot(sigma_rr,rnodes,'-r',markevery=[0])
plt.title('Displacement_diffrent nodes_elastic_area')
plt.xlabel('Displacement[m]')
plt.ylabel('nodes_length[m]')

plt.savefig("elastic.png")
plt.close()
