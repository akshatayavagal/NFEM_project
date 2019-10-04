import numpy as np 
t_ini = 0
t_total = 1
delta_t = 0.17
t_steps = []
for i in range(int(t_total/delta_t)+1):
    x = i*delta_t
    t_total = t_total-x
    t_steps.append(x)
    print(t_steps)




#no_steps = np.linspace(0,t_total+n_s,n_ss+2)
