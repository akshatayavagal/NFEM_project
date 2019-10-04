import numpy as np
from shape_matrix import *
from dN_dz_matrix import *
from Jacobian_matrix import *
from B_matrix import *
from Assignment_matrix import *
from material_routine import *
from Parameters import*


def Element_routine(u_ele,t,x_ele,elemental_plastic_strain):
    ### gauss point intialization 
    z = 0
    #### call shape matrix function
    N = shape_matrix(z)
    #### call dn_dz matrix function
    dN_dZ = dN_dz_matrix()
    #### call jacobian matrix function
    J = Jacobian_matrix(x_ele,dN_dZ)
    #### call B matrix function
    B = B_matrix(z,x_ele,dN_dZ,J)
    B = np.array(B)
    strain_current=B.dot(u_ele)

    ##### call material routine
    C_t,elemental_plastic_strain,sigma_next_step = material_routine(lemda,mue,strain_current,elemental_plastic_strain,sigma_yield)
    
    ##### assign the weight function to gauss point
    if z == 0:
        w_alpha = 2
    
    ##### get the element stiffness matrix 
    N_r=np.matmul(np.transpose(N),x_ele)
    K_t_ele = w_alpha*np.transpose(B).dot(C_t).dot(B)*N_r**2 * J

    ##### get the element internal force matrix
    F_int_ele = w_alpha*np.transpose(B).dot(sigma_next_step)*N_r**2*J

    return K_t_ele,F_int_ele,elemental_plastic_strain,sigma_next_step 

