                        #####       variant = 8        ######

#####     Material Parameters         #####

youngs_modulus_E = 70000                   # mpa
poission_ratio_v = 0.30                   # unit_less
sigma_yield = 140                          # mpa
volumetric_strain = 0.02                   # unit_less

#####       Equations                 ##### 
            
lemda = poission_ratio_v*youngs_modulus_E/(1-2*poission_ratio_v)/(1+poission_ratio_v)
mue = youngs_modulus_E/2/(1+poission_ratio_v)

##### parameters require for meshing  #####

nelem=10                                   # unit_less
ri = 15                                # micro_m
ro = 75                                 # micro_m

#####          times steps            #####      

delta_t = 0.1
t_ini = 0.0
t_total = 1.0