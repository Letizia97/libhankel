import libhankel
import numpy as np

nu = 0

# Tabulated form factor data: sphere-like form factor
# q: scattering vector (1/Angstrom)
# f: form factor intensity
q_data = np.array([0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0])
f_data = np.array([1.0, 0.98, 0.92, 0.72, 0.48, 0.28, 0.15, 0.05, 0.01, 0.002])

# Transform evaluation points
x_arr = np.array([
    1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    21.0, 22.0, 23.0, 24.0, 25.0
])

# Strategy parameters
strategy_p_dict = {"n_eval": 250, "eps_rel": 1e-9}

# Call hankel_transform with tabulated form factor as a dict
# exponent is optional: pass a value for power_law tail, or omit to auto-fit
result = libhankel.hankel_transform(
    nu,
    {
        "q": q_data,
        "f": f_data,
        "interp_type": "linear",
        "tail": "zero"
    },
    x_arr,
    [],  # no parameters needed
    "QWE_Chave",
    strategy_p_dict
)

print("Result of Hankel transform with tabulated form factor:")
print(result)
