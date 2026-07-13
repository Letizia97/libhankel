import libhankel
import numpy as np

nu = 0
params_broad_peak = np.array([10e5, 1000, 0.01, 2, 2])
x_arr = np.array([
    30.,   49.,   69.,   88.,   108.,  127.,  147.,  167.,  186.,  206.,  225., 
    225.,  245.,  265.,  284.,  304.,  323.,  343.,  362.,  382.,  402.,  421.,  
    441.,  460.,  480.,  
])

strategy_p_dict = {"n_eval": 250, "eps_rel": 1e-9}

result = libhankel.hankel_transform(
    nu,
    "broad_peak",
    x_arr,
    params_broad_peak,
    "QWE_Chave",
    strategy_p_dict,
)

print("Result of calling the hankel_transform with broad_peak form factor", result)
