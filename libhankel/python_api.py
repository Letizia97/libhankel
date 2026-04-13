from ._clib import libhankel

from ctypes import (
    c_double
)


def form_factor_g_dab(
    single_q: float,
    param_array: list[float]
) -> float:   
    """
    Computes the gdab form factor on a single q.

    Parameters
    ----------
    single_q : float
        A single q where to compute the form factor.
    param_array : list[float]
        List of parameters needed by this form factor. They must be in the following order:
        -   XI  correlation length (e.g., 10.0)
        -   H   Hurst exponent (e.g., 0.5)
        -   ETA scattering length density contrast (e.g., 1e-4)        
        
    Returns
    -------
    float
        The resulting form factor value.
    """
    return libhankel.form_factor_g_dab(
        c_double(single_q), 
        (c_double * 50)(*param_array)
    )








# def hankel_transform(
#     nu,
#     form_factor_f,
#     x_arr,
#     f_params, 
#     output, 
#     len_x, 
#     strategy_name,
#     strategy_params 
# ):
#     """
#     Inputs
#         nu:    either 0 or 1
    
#     """

#     res = libhankel.hankel_transform(
#         c_int(nu),

#     )
#     return res

