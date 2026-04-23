import numpy as np
import libhankel

nu = 0
params_spheres = np.array([10.0,  1.0])
x_arr = np.array([1.0,  2.0,  3.0,  4.0,  5.0,  6.0])

strategy_p_dict = {
    'n_eval': 250, 
    'eps_rel': 1e-9
}

result = libhankel.hankel_transform(
    nu,
    libhankel.form_factor_sphere_callback,
    x_arr,
    params_spheres,
    "QWE_Chave",
    strategy_p_dict
)

print(
    "Result of calling the hankel_transform with sphere form factor", 
    result
)