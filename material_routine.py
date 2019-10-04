import numpy as np
from Parameters import*



def material_routine(lemda,mue,strain_current,elemental_plastic_strain,sigma_yield):
    ##### material tangent for elastic matrix
    def material_tangent_elastic_matrix(lemda,mue):
        Ct11 = Ct22 = Ct33 = lemda + 2*mue
        Ct12 = Ct13 = Ct23 = Ct21 = Ct31 = Ct32 = lemda
        C_elastic= np.array([[Ct11,Ct12,Ct13],[Ct21,Ct22,Ct23],[Ct31,Ct32,Ct33]])
        return C_elastic

    C_elastic = material_tangent_elastic_matrix(lemda,mue)
    #### calculation of stress trial
    sigma_trial=C_elastic.dot(np.subtract(strain_current,elemental_plastic_strain))
    sigma_trial_sprl=np.sum(sigma_trial)*(1/3)
    sigma_trial_devi=np.array(sigma_trial)-sigma_trial_sprl
    sigma_equi=np.sqrt(np.transpose(sigma_trial_devi).dot(sigma_trial_devi)*(1.5))
    #### calculation of stress equivivalent
    sigma_equi_trial=np.asscalar(sigma_equi)
    #### calculation of elastic or plastic region of material
    phi_cal = sigma_equi_trial - sigma_yield
    if phi_cal<=0:
        
        delta_lemda=0
        ####calculation of Ct for elastic 
        Ct=C_elastic
        sigma_next_step=sigma_trial
        with open('output_1.txt','a') as f:
            f.write("sigma_next_step =\n %s\n" %str(sigma_next_step))
        element_plastic_strain=elemental_plastic_strain+delta_lemda*np.sign(sigma_next_step)
        # with open('output_1.txt','a') as f:
        #     f.write("mode : elastic \n  p_strain =\n %s\n" %str((element_plastic_strain)))

    else:

        delta_lemda = (sigma_equi_trial-sigma_yield) / (3*mue)
        ####calculation Ct for plastic region
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
        #### stress calculation for plastic region 

        sigma_next_step = sigma_trial_sprl+(((sigma_equi_trial-(3*mue*delta_lemda))/sigma_equi_trial)*sigma_trial_devi)
        sigma_next_devi = ((sigma_equi_trial-3*mue*delta_lemda)/sigma_equi_trial)*sigma_trial_devi
        sigma_next_equi = (sigma_equi_trial - 3*mue*delta_lemda)

        with open('output_1.txt','a') as f:
            f.write("sigma_next_step =\n %s\n" %str(sigma_next_step))
        ####plastic strain calculation 
        elemental_plastic_strain=elemental_plastic_strain+delta_lemda*(1.5)*sigma_next_devi/sigma_next_equi
        # with open('output_1.txt','a') as f:
        #     f.write("mode : plastic \n  p_strain =\n %s\n" %str((elemental_plastic_strain)))
    
    return Ct,elemental_plastic_strain,sigma_next_step