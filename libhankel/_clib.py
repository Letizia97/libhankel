from ctypes import (
    CDLL, 
    c_int, 
    c_double, 
    POINTER, 
    CFUNCTYPE, 
    c_char_p,
    Structure,
    c_size_t
)

libhankel = CDLL("build/src/liblibhankel.so")

Array50 = c_double * 50
PtrToArray50 = POINTER(Array50)
FORM_FACTOR_FUNC_TYPE = CFUNCTYPE(
    c_double,      # return type
    c_double,      # first arg
    PtrToArray50   # second arg: double (*)[50]
)

class StrategyParams(Structure):
    _fields_ = [
        ("n_eval", c_int),
        ("eps_rel", c_double),
        ("f_max", c_double),
    ]


libhankel.hankel_transform.argtypes = [
    c_int,
    FORM_FACTOR_FUNC_TYPE, 
    POINTER(c_double),
    PtrToArray50,
    POINTER(c_double),
    c_size_t,
    c_char_p,
    StrategyParams
]
libhankel.hankel_transform.restype = c_int


libhankel.form_factor_g_dab.argtypes = [
    c_double,       # first arg
    PtrToArray50    # second arg: double (*)[50]
]
libhankel.form_factor_g_dab.restype = c_double

libhankel.form_factor_broad_peak.argtypes = [
    c_double,       # first arg
    PtrToArray50    # second arg: double (*)[50]
]
libhankel.form_factor_broad_peak.restype = c_double

libhankel.form_factor_sphere.argtypes = [
    c_double,       # first arg
    PtrToArray50    # second arg: double (*)[50]
]
libhankel.form_factor_sphere.restype = c_double


form_factor_g_dab_callback = FORM_FACTOR_FUNC_TYPE(
    ("form_factor_g_dab", libhankel)
)
form_factor_broad_peak_callback = FORM_FACTOR_FUNC_TYPE(
    ("form_factor_broad_peak", libhankel)
)
form_factor_sphere_callback = FORM_FACTOR_FUNC_TYPE(
    ("form_factor_sphere", libhankel)
)



