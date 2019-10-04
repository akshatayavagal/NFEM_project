import numpy as np
import math as math
from Parameters import *
import matplotlib.pyplot as plt
#from mesh_generator import *
#  file_name = mesh_generator.py
#  Generate list of position of nodes according to a geometric series
#     for assignement in "Nonlinear Finite Element Methods" 
#     in summer term 2019
#     lecturer in charge: Dr. Geralf H�tter

# Input Parameters
#ro=75e-6                      #outer radius
#ri=15e-6                      #radius of inclusion
#nelem=10                      #number of elements
 


# visualize location of nodes 
#plt.plot(rnodes,np.zeros([nelem+1]),'*')
#plt.show()
#######################################################################################
#from shape_matrix import *
def shape_matrix(z):
    N1 = 1/2 * (1-z)
    N2 = 1/2 * (1+z)
    N = np.array([[N1],[N2]])
    return N
#####################################################################################
#from dN_dz_matrix import *
def dN_dz_matrix():
    dN_dZ1 = -1/2
    dN_dZ2 =  1/2
    dN_dZ = np.array([[dN_dZ1],[dN_dZ2]])
    return dN_dZ
########################################################################################
#from Jacobian_matrix import *
def Jacobian_matrix(x_ele,dN_dZ):
    # Jacobian_matrix  
    dN_dZ = dN_dz_matrix()
    J = np.matmul(np.transpose(dN_dZ),x_ele)
    return J
#########################################################################################
#from B_matrix import *
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
############################################################################################
#from Assignment_matrix import *



##################################################################################################################
def Element_routine(u_ele,t,x_ele,elemental_plastic_strain): 
    z = 0
    N = shape_matrix(z)
    dN_dZ = dN_dz_matrix()
    J = Jacobian_matrix(x_ele,dN_dZ)
    B = B_matrix(z,x_ele,dN_dZ,J)
    B = np.array(B)
    strain_current=B.dot(u_ele)
    C_t,elemental_plastic_strain,sigma_next_step = material_routine(lemda,mue,strain_current,elemental_plastic_strain,sigma_yield)
    
    if z == 0:
        w_alpha = 2
    
    ##### K_t_ele 
    N_r=np.matmul(np.transpose(N),x_ele)
    K_t_ele = w_alpha*np.transpose(B).dot(C_t).dot(B)*N_r**2 * J

    ##### F_int_ele
    #F_int_ele = np.matmul(K_t_ele,u_ele)
    F_int_ele = w_alpha*np.transpose(B).dot(sigma_next_step)*N_r**2*J

    return K_t_ele,F_int_ele,elemental_plastic_strain,sigma_next_step 

###############################################################################################################

def material_routine(lemda,mue,strain_current,elemental_plastic_strain,sigma_yield):

    def material_tangent_elastic_matrix(lemda,mue):
        Ct11 = Ct22 = Ct33 = lemda + 2*mue
        Ct12 = Ct13 = Ct23 = Ct21 = Ct31 = Ct32 = lemda
        C_elastic= np.array([[Ct11,Ct12,Ct13],[Ct21,Ct22,Ct23],[Ct31,Ct32,Ct33]])
        return C_elastic

    C_elastic = material_tangent_elastic_matrix(lemda,mue)
    sigma_trial=C_elastic.dot(np.subtract(strain_current,elemental_plastic_strain))
    sigma_trial_sprl=np.sum(sigma_trial)*(1/3)
    sigma_trial_devi=np.array(sigma_trial)-sigma_trial_sprl
    sigma_equi=np.sqrt(np.transpose(sigma_trial_devi).dot(sigma_trial_devi)*(1.5))
    sigma_equi_trial=np.asscalar(sigma_equi)
    phi_cal = sigma_equi_trial - sigma_yield
    #print(phi_cal)
    if phi_cal<=0:
        
        delta_lemda=0
        Ct=C_elastic
        sigma_next_step=sigma_trial
        with open('output_1.txt','a') as f:
            f.write("sigma_next_step =\n %s\n" %str(sigma_next_step))
        element_plastic_strain=elemental_plastic_strain+delta_lemda*np.sign(sigma_next_step)
        with open('output_1.txt','a') as f:
            f.write("mode : elastic \n  p_strain =\n %s\n" %str((element_plastic_strain)))

    else: 
        # with open('output_1.txt','a') as f:
        #     f.write("mode : plastic\n")

        delta_lemda = (sigma_equi_trial-sigma_yield) / (3*mue)
        C_1=(3*lemda+2*mue)/3
        C_2=2*mue*(sigma_equi_trial-3*mue*delta_lemda)/sigma_equi_trial
        C_3=np.array(3*mue*(sigma_trial_devi).dot(np.transpose(sigma_trial_devi))/(sigma_equi_trial**2))
        Ct=C_1-C_3
        for i in range(3):
            for j in range(3):
                if i==j:
                    Ct[i][j]=Ct[i][j]+C_2*(2/3)
                else:
                    Ct[i][j]=Ct[i][j]+C_2*(-1/3)
        sigma_next_step=sigma_trial_sprl+(((sigma_equi_trial-(3*mue*delta_lemda))*sigma_trial_devi)/sigma_equi_trial)
        
        sigma_next_sprl=np.sum(sigma_next_step)*(1/3)
        sigma_next_devi=np.array(sigma_next_step)-sigma_next_sprl
        sigma_next_equi=np.sqrt(np.transpose(sigma_next_devi).dot(sigma_next_devi)*(1.5))
        sigma_equi_trial=np.asscalar(sigma_next_equi)

        
        #print(sigma_next_step)
        with open('output_1.txt','a') as f:
            f.write("sigma_next_step =\n %s\n" %str(sigma_next_step))
        elemental_plastic_strain=elemental_plastic_strain-delta_lemda*(1.5)*sigma_next_devi/sigma_equi_trial
        #element_plastic_strain=strain_current-(np.linalg.inv(C_elastic).dot(sigma_next_step))
        with open('output_1.txt','a') as f:
            f.write("mode : plastic \n  p_strain =\n %s\n" %str((elemental_plastic_strain)))
    
    return Ct,elemental_plastic_strain,sigma_next_step

#######################################################################################################################

def main_program(i):
    no_steps = math.ceil(t_total/delta_t)
    U = [np.zeros([nelem+1,1])]
    # Plastic strain initialization
    plastic_strain = [np.zeros([nelem+1,3,1])]
    sigma_history = []
    ### main loop for 
    for m in range(no_steps+1):
        if m == no_steps:
            t = 1
        else:
            t = t_ini + m*delta_t

        with open('output_1.txt','a') as f:
            f.write("New itteration with value of t = %s\n" %str((t)))
                
        if t == 0:
            continue
        U_global = U[-1]
        U_global[0] = (1/3)*t*volumetric_strain*ri
        k=1
        while True:
            current_plastic_strain = plastic_strain[-1]
            #print(current_plastic_strain)
            K_t = np.zeros([nelem+1,nelem+1])
            G = np.zeros([nelem+1,1])
            plstic_strain_history = []
            sigma = []
            for i in range(nelem):
                with open('output_1.txt','a') as f:
                    f.write("Element no. = %s\n" %str((i)))
                elemental_plastic_strain=(current_plastic_strain[i])
                #print(elemental_plastic_strain)
                u_ele = np.matmul(Assignment_matrixs[i],U_global)
                x_ele = [[rnodes[i]],[rnodes[i+1]]]
                K_t_ele,F_int_ele,elemental_plastic_strain,sigma_next_step = Element_routine(u_ele,t,x_ele,elemental_plastic_strain)
                sigma.append(sigma_next_step)
                plstic_strain_history.append(elemental_plastic_strain)
                # with open('output_1.txt','a') as f:
                #     f.write("element_plastic_strain : %s\n" %str((element_plastic_strain)))
                #####

                A_Kt = np.matmul(np.transpose(Assignment_matrixs[i]),K_t_ele)
                A_kt_A = np.matmul(A_Kt,Assignment_matrixs[i])
                K_t = K_t + A_kt_A


                #####
                A_F_int = np.matmul(np.transpose(Assignment_matrixs[i]),F_int_ele)
                G = G + A_F_int

            
        
            delta_U = np.zeros([nelem+1,1])
            delta_U[0] = 0 
            #delta_U[0] = (1/3)*t*volumetric_strain*rnodes[0]
            K_t_reduce = K_t[1:,1:]
            G_reduce = G[1:]

            #####
            delta_U[1:] = np.matmul(np.linalg.inv(K_t_reduce),-G_reduce)
            #delta_U = delta_U.insert(0,delta_U_0)
            U_global = U_global + delta_U
            with open('output_1.txt','a') as f:
                f.write("value of k is : %s\n" %str((k)))
            with open('output_1.txt','a') as f:
                f.write(" value of U_global \n: %s\n" %str((U_global)))
            
    #value of k is: %s\n  
            
            # print('value of k is :',k)
            if  np.linalg.norm(G_reduce,np.inf) < 0.005*np.linalg.norm(A_F_int,np.inf) or np.linalg.norm(delta_U,np.inf) < 0.005*np.linalg.norm(U_global,np.inf) or k>5:
                U.append(U_global)
                break
            U.append(U_global)
            k=k+1
        plastic_strain.append(plstic_strain_history)
        sigma_history.append(sigma)
        U_s,sigma_rr,gauss_local = analytical_solution(nelem)
    return plastic_strain,sigma_history,U,U_s,sigma_rr,gauss_local



#####################################################################################################################

# print(U[-1])
# print(len(sigma_history))
# print(sigma_history[-1])

# sigma_rr = []
# for i in sigma_history[-1]:
#     fi = np.asscalar(i[0])
#     sigma_rr.append(fi)
# print(sigma_rr)
#############################################################################

def analytical_solution(nelem):
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

#############################################################################################
fig, ax=plt.subplots(ncols=2,figsize=(10, 10))

for i in np.array([10]):
    nelem = i
    
    ####################################################################
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
    Assignment_matrixs = []
    for i in range(nelem):
        a = np.zeros((2,nelem+1))
        a[0][i]=1
        a[1][i+1]=1
        Assignment_matrixs.append(a)
    
    
    
    plastic_strain,sigma_history,U,U_s,sigma_rr,gauss_local = main_program(i)
    sigma_rrr = []
    for i in sigma_history[-1]:
        fi = np.asscalar(i[0])
        sigma_rrr.append(fi)
    #print(len(sigma_rrr))
    #print(len(gauss_local))
    ax[0].plot(rnodes,U[-1],'*--', label = 'Elements:{}'.format(nelem))
    ax[0].set(xlabel = 'r [mm]',ylabel='U_r [mm]',title = 'Displacement_study')
    ax[0].legend()
    ax[1].plot(gauss_local,sigma_rrr,'*--',label = 'Elements:{}'.format(nelem))
    ax[1].set(xlabel = 'r_c [mm]',ylabel='sigma_rr [mm]',title = 'stress_study')
    ax[1].legend()
    if nelem == 20:
        ax[0].plot(rnodes,U_s,'o-',label = 'Analytical_Elements:{}'.format(nelem))
        ax[0].legend()
        ax[1].plot(gauss_local,sigma_rr,'o-',label = 'Analytical_Elements:{}'.format(nelem))
        ax[1].legend()
    

    # plt.subplot(2,1,2)
    # sigma_rrr = []
    # for i in sigma_history[-1]:
    #     fi = np.asscalar(i[0])
    #     sigma_rrr.append(fi)
    # plt.plot(sigma_rr,np.arange(nelem),'-o',markevery=[0])
    # plt.plot(sigma_rrr,np.arange(nelem),'-r',markevery=[0])
plt.savefig("elastic_convergence")   
plt.show()
    






##########################################################################
# fig, ax=plt.subplots(figsize=(10, 10))
# plt.subplot(2,1,1)
# plt.plot(U_s,rnodes,'-o',markevery=[0])
# plt.plot(U[-1],rnodes,'-r',markevery=[0])
# plt.title('Displacement_diffrent nodes_elastic_area')
# plt.xlabel('Displacement[m]')
# plt.ylabel('nodes_length[m]')
# plt.savefig("elastic_3.png")
# plt.show()