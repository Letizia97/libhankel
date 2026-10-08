import sys

from setuptools import Extension, setup

include_dirs = [
    "include",
    "external_libs",
    ".",
    "src/utils",
    "external_libs/utils",
]

# Boost.Math is used header-only. Its headers are on the compiler's default
# search path on Linux, but not on macOS, where Homebrew installs under
# /opt/homebrew (Apple Silicon) or /usr/local (Intel). For anywhere else, pass
# CPPFLAGS=-I/path/to/boost.
if sys.platform == "darwin":
    include_dirs += ["/opt/homebrew/include", "/usr/local/include"]
elif sys.platform == "win32":
    include_dirs += ["C:/vcpkg/installed/x64-windows/include"]

module = Extension(
    "libhankel",
    define_macros=[("LIBHANKEL_BUILDING_DLL", None)],
    sources=[
        "src/py_interface.c",
        "src/form_factors.c",
        "src/hankel_transform.c",
        "src/hankel_DE_quadrature.c",
        "src/hankel_DHT.c",
        "src/hankel_G0.c",
        "src/hankel_QWE.c",
        "external_libs/DE-quadrature/intde.c",
        "external_libs/qwe_Chave/qwe_Chave.c",
        "external_libs/qwe_Key/qwe_Key.c",
        "external_libs/utils/tanhsinh.c",
        "src/utils/sasfit_integrate.c",
        "src/utils/sf_functions.c",
        "src/utils/validate_x.c",
        "src/utils/boost_bessel.cpp",
        "src/utils/tabulated_ff.c",
        "src/utils/interp_linear.c",
        "src/utils/interp_cubic.cpp",
        "src/utils/interp_loglinear.c",
    ],
    include_dirs=include_dirs,
)

setup(
    name="libhankel",
    version="0.1",
    ext_modules=[module],
)
