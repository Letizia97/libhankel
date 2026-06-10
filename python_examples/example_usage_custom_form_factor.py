import numpy as np
import math
import libhankel

nu = 0
params_gdab = np.array([10.0,  0.5, 1e-4])
x_arr = np.array([
    15.,          18.54166667,  22.08333333,  25.625,       29.16666667,
    32.70833333,  36.25,        39.79166667,  43.33333333,  46.875,
    50.41666667,  53.95833333,  57.5,         61.04166667,  64.58333333,
    68.125,       71.66666667,  75.20833333,  78.75,        82.29166667,
    85.83333333,  89.375,       92.91666667,  96.45833333,  100.,   
])

strategy_p_dict = {
    'n_eval': 250, 
    'eps_rel': 1e-9
}


def sf_poch(a, x):
    return math.gamma(a + x) / math.gamma(a)

def dab(q, params):
    """
    Example Debye-Anderson-Brumberger (DAB) model.
    """
    XI = params[0]
    H = params[1]
    ETA = params[2]

    factor1 = (2 * XI) ** 3
    factor2 = sf_poch(H, 1.5)

    numer = (factor1 * factor2 * ETA) ** 2 * (np.pi ** 3)
    denom = (1 + (q * XI) ** 2) ** (1.5 + H)

    return numer / denom



result = libhankel.hankel_transform(
    nu,
    "gdab",
    x_arr,
    params_gdab,
    "QWE_Chave",
    strategy_p_dict
)

print(
    "Result of calling the hankel_transform "
    "with g_dab form factor (BUILT IN)\n", 
    result
)


result = libhankel.hankel_transform(
    nu,
    dab,
    x_arr,
    params_gdab,
    "QWE_Chave",
    strategy_p_dict
)
print(
    "Result of calling the hankel_transform with g_dab "
    "form factor (USER-DEFINED python function)\n", 
    result
)