import libhankel
import numpy as np

nu = 0
params_gdab = np.array([10.0, 0.5, 1e-4])
x_arr = np.array([
    15.,          18.54166667,  22.08333333,  25.625,       29.16666667,
    32.70833333,  36.25,        39.79166667,  43.33333333,  46.875,
    50.41666667,  53.95833333,  57.5,         61.04166667,  64.58333333,
    68.125,       71.66666667,  75.20833333,  78.75,        82.29166667,
    85.83333333,  89.375,       92.91666667,  96.45833333,  100.,
])

# FFTLog takes no strategy parameters; pass an empty dict
result = libhankel.hankel_transform(nu, "gdab", x_arr, params_gdab, "FFTLog", {})

print("Result of calling hankel_transform with FFTLog strategy:\n", result)
