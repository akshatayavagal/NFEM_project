import numpy as np
import math as math
import sys
from Parameters import *
from mesh_generator import *
from shape_matrix import *
from dN_dz_matrix import *
from Jacobian_matrix import *
from B_matrix import *
from Assignment_matrix import *
from element_routine import *


def main_program(nelem):
    # take initial values from the parameters
    no_steps = math.ceil(t_total/delta_t)
    # initialization of global displacement
    U = [np.zeros([nelem+1,1])]
    # Plastic strain initialization
    plastic_strain = [np.zeros([nelem+1,3,1])]
    # make list for sigma history for every time step
    sigma_history = []
    Assignment_matrixs = Assignment_matrix(nelem)
    # get rnodes list which contain all the nodes positions
    rnodes = mesh_generator(nelem)



    # main loop for every iteration of time steps:
    for m in range(no_steps+1):
        ## if_else loop to take full loading when t_total is not equally divided by delta_t
        if m == no_steps:
            t = 1
        else:
            t = t_ini + m*delta_t

        with open('output_1.txt','a') as f:
            f.write("New itteration with value of t = %s\n" %str((t)))
        ## if time step is equal to t_ini this loop skips first itteration for time step
        if t == 0:
            continue
        ## this one takes the value of gloabal displacement from last itteration.
        U_global = U[-1].copy()
        ## boundary condition at first node of the first element.
        U_global[0] = (1/3)*t*volumetric_strain*ri

        ## intialization of the Newton raphson scheme:
        k=1
        while True:
            ### take plastic strain value from last time step itteration.
            current_plastic_strain = plastic_strain[-1]
            ### initialization of gloabal stiffness matrix and hlobal G_matrix
            K_t = np.zeros([nelem+1,nelem+1])
            G = np.zeros([nelem+1,1])
            plstic_strain_history = []
            sigma = []

            ### loop for itteration over every element:
            for i in range(nelem):
                with open('output_1.txt','a') as f:
                    f.write("Element no. = %s\n" %str((i)))
                ### take values of plastic strain from list according to element
                elemental_plastic_strain=(current_plastic_strain[i])
                ### find current displcement over approximation
                u_ele = np.matmul(Assignment_matrixs[i],U_global)
                ### takes position of both nodes of relevent element from meshgenerator.py
                x_ele = [[rnodes[i]],[rnodes[i+1]]]
                ### call Element Routine
                K_t_ele,F_int_ele,elemental_plastic_strain,sigma_next_step = Element_routine(u_ele,t,x_ele,elemental_plastic_strain)
                ### append the values of elements stresses and nodes displacements
                sigma.append(sigma_next_step)
                plstic_strain_history.append(elemental_plastic_strain)
                ### find global stiffness matrix from element stiffness matrix
                A_kt_A = np.transpose(Assignment_matrixs[i]).dot(K_t_ele).dot(Assignment_matrixs[i])
                K_t = K_t + A_kt_A
                ### find global F_int from the elemental f_int
                A_F_int = np.matmul(np.transpose(Assignment_matrixs[i]),F_int_ele)
                ### update G_matrix
                G = G + A_F_int


            #### initialize the delta displacement
            delta_U = np.zeros([nelem+1,1])
            #### reduced the global stiffness matrix and global G_matrix
            K_t_reduce = K_t[1:,1:]
            G_reduce = G[1:]

            #### find delta displcement from the below equation.
            delta_U[1:] = np.matmul(np.linalg.inv(K_t_reduce),-G_reduce)
            #### Update the U_global by adding the delta displcaement
            U_global = U_global + delta_U
            with open('output_1.txt','a') as f:
                f.write("value of k is : %s\n" %str((k)))
            with open('output_1.txt','a') as f:
                f.write(" value of U_global \n: %s\n" %str((U_global)))

            #### condtion for newton raphson scheme
            if  np.amax(abs(G_reduce)) < 0.005*np.amax(abs(A_F_int)) or np.amax(abs(delta_U)) < 0.005*np.amax(abs(U_global)) or k>5:
                U.append(U_global)
                #### if value of k is 6 this program doesnt run further
                if k==6:
                    sys.exit()

                break
            #### append U with U_global to keep history of displcements
            U.append(U_global)
            k=k+1

        #### append the values of plastic strain and stresses
        plastic_strain.append(plstic_strain_history)
        sigma_history.append(sigma)

    return plastic_strain,sigma_history,U

#### function to get the Analytical values of the stresses and displacement.
def analytical_solution(nelem):
    rnodes = mesh_generator(nelem)
    U_s = []
    sigma_rr=[]
    gauss_local = []
    for j in range(nelem+1):
        r=rnodes[j]

        if j < nelem:
            gauss_l = (rnodes[j]+rnodes[j+1])/2
            gauss_local.append(gauss_l)
            sigma_verif=(-2*youngs_modulus_E*volumetric_strain*gauss_local[0]**3)/(3*(1+poission_ratio_v)*gauss_l**3)
            sigma_rr.append(sigma_verif)
        U_ana = ((ri**3)*volumetric_strain)/(3*r*r)

        U_s.append(U_ana)

    return U_s,sigma_rr,gauss_local

U_s,sigma_rr,gauss_local = main_program(nelem)