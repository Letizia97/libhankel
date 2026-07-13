import libhankel
import numpy as np

nu = 0
params_broad_peak = np.array([10e5, 1000, 0.01, 2, 2])
x_arr = np.array(
    [
        30.0,
        49.0,
        69.0,
        88.0,
        108.0,
        127.0,
        147.0,
        167.0,
        186.0,
        206.0,
        225.0,
        225.0,
        245.0,
        265.0,
        284.0,
        304.0,
        323.0,
        343.0,
        362.0,
        382.0,
        402.0,
        421.0,
        441.0,
        460.0,
        480.0,
    ]
)

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
