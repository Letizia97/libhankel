import libhankel
import numpy as np

nu = 0

# Tabulated form factor data: sphere-like form factor
# q: scattering vector (1/Angstrom)
# f: form factor intensity
q_data = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
f_data = np.array([100.0, 99.5, 97.5, 90.0, 70.0, 20.0, 5.0, 1.0, 0.1, 0.01])

# Transform evaluation points
x_arr = np.array([
    0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0
])

# Strategy parameters
strategy_p_dict = {"n_eval": 100, "eps_rel": 1e-6}

# Call hankel_transform with tabulated form factor as a dict
# exponent is optional: pass a value for power_law tail, or omit to auto-fit
result = libhankel.hankel_transform(
    nu,
    {
        "q": q_data,
        "f": f_data,
        "interp_type": "cubic",
        "tail": "power_law"
    },
    x_arr,
    [],  # no parameters needed
    "QWE_Chave",
    strategy_p_dict
)

print("Result of Hankel transform with tabulated form factor:")
print(result)
