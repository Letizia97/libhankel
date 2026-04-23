
# Standard library imports
import ctypes
import warnings
from ctypes import c_double, c_int, c_size_t, POINTER, CFUNCTYPE, CDLL

# Third-party imports
import numpy as np

# Local application/library imports
from ._clib import libhankel, StrategyParams


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

def form_factor_broad_peak(
    single_q: float,
    param_array: list[float]
) -> float:   
    """
    Computes the broad peak form factor on a single q.

    Parameters
    ----------
    single_q : float
        A single q where to compute the form factor.
    param_array : list[float]
        List of parameters needed by this form factor. They must be in the following order:
        -   I0  
        -   XI  correlation length (e.g., 10.0)   
        -   Q0
        -   M
        -   P
     
    Returns
    -------
    float
        The resulting form factor value.
    """
    return libhankel.form_factor_broad_peak(
        c_double(single_q), 
        (c_double * 50)(*param_array)
    )

def form_factor_sphere(
    single_q: float,
    param_array: list[float]
) -> float:   
    """
    Computes the sphere form factor on a single q.

    Parameters
    ----------
    single_q : float
        A single q where to compute the form factor.
    param_array : list[float]
        List of parameters needed by this form factor. They must be in the following order:
        -   R   radius (> 0)
        -   ETA scattering length density contrast (e.g., 1e-4)  
      
    Returns
    -------
    float
        The resulting form factor value.
    """
    return libhankel.form_factor_sphere(
        c_double(single_q), 
        (c_double * 50)(*param_array)
    )


def hankel_transform(
    nu,
    form_factor_f,
    x_arr,
    f_params, 
    strategy_name,
    strategy_params_dict, 
):
    """
    Computes the Hankel transform. 

    Parameters
    ----------
        nu : int              
            Order of the Bessel function. Either 0 or 1.

        form_factor_f : Callable[[float, list[float]], float]   
            Function to Hankel transform i.e. the form factor
            (can use e.g. form_factor_g_dab).

        x_arr : ndarray            
            Numpy array of x values.

        f_params : ndarray of shape (n,) , with n < 50    
            Array containing the form factor parameters (ordered as per docs). 

        strategy_name : str
            The strategy name - consult docs for options. 

        strategy_params_dict : dict[str, float]
            The parameters needed for the selected strategy to work 
            (please consult the docs to know which params are required by which strategy)

    """

    # Prepare array data x
    N = x_arr.size
    ptr_x_arr = ctypes.cast(
        x_arr.ctypes.data,
        ctypes.POINTER(c_double)
    )
    
    # Prepare form factor 
    CArrayType50 = c_double * 50
    CALLBACK = ctypes.CFUNCTYPE(
        c_double,
        c_double,
        ctypes.POINTER(CArrayType50)
    )
    c_form_factor_f = CALLBACK(form_factor_f)

    # Prepare f_params
    params = np.asarray(f_params, dtype=np.float64)
    params = np.ascontiguousarray(params)

    c_array_f_params = CArrayType50()
    for i in range(len(params)):
        c_array_f_params[i] = params[i]

    # Prepare output array
    output_arr = np.zeros_like(x_arr)
    ptr_output_arr = ctypes.cast(
        output_arr.ctypes.data,
        ctypes.POINTER(c_double)
    )

    # Prepare strategy_name string
    ptr_strategy_name = ctypes.c_char_p(strategy_name.encode("utf-8"))

    # Exclude params not needed
    field_names = {name for name, _ in StrategyParams._fields_}
    ignored = set(strategy_params_dict) - field_names
    if ignored:
        warnings.warn(
            f"Ignoring unknown StrategyParams fields: {sorted(ignored)}",
            UserWarning,
            stacklevel=2,
        )

    # Prepare parameters for the strategy
    strategy_params_struct = StrategyParams(**{
        k: v for k, v in strategy_params_dict.items()
        if k in field_names
    })

    # Compute transform
    status = libhankel.hankel_transform(
        c_int(nu),
        c_form_factor_f,
        ptr_x_arr,
        c_array_f_params,
        ptr_output_arr,
        c_size_t(N),
        ptr_strategy_name,
        strategy_params_struct
    )

    if status != 0:
        raise Exception(f"Compute failed with status {status}")

    res = np.ctypeslib.as_array(ptr_output_arr, shape=(N,))
    return res

