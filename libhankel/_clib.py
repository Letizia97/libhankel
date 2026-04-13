from ctypes import (
    CDLL, 
    c_int, 
    c_double, 
    POINTER, 
    CFUNCTYPE, 
    c_char_p,
    Structure
)

libhankel = CDLL("liblibhankel.so")


###############################################
# For hankel_transform function
###############################################
Array50 = c_double * 50
PtrToArray50 = POINTER(Array50)

FUNC_TYPE = CFUNCTYPE(
    c_double,      # return type
    c_double,      # first arg
    PtrToArray50   # second arg: double (*)[50]
)

class strategy_params(Structure):
    _fields_ = [
        ("n_eval", c_double),
        ("eps_rel", c_double),
        ("f_max", c_double)
    ]

# libhankel.hankel_transform.argtypes = [
#     c_int,
#     FUNC_TYPE, 
#     POINTER(c_double),
#     PtrToArray50,
#     POINTER(c_double),
#     c_int,
#     c_char_p,
#     strategy_params
# ]

# libhankel.hankel_transform.restype = c_int


###############################################
# For form_factor_g_dab function
###############################################
libhankel.form_factor_g_dab.argtypes = [
    c_double,       # first arg
    PtrToArray50    # second arg: double (*)[50]
]
libhankel.form_factor_g_dab.restype = c_double


###############################################
###############################################
