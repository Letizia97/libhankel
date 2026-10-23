import libhankel
import numpy as np

params_gdab = np.array([10.0, 0.5, 1e-4])
epsrel = 1e-6

result = libhankel.hankel_G0("gdab", params_gdab, epsrel)

print("Result of calling hankel_G0 with the built-in g_dab form factor:", result)
