import os
import shutil
import subprocess

from setuptools import Extension, setup


def boost_include_dirs():
    """Locate the header-only Boost.Math headers used by src/utils/boost_bessel.cpp.

    Returns the extra include directories needed to compile against Boost, or
    an empty list when the compiler can already find it.  Debian and Ubuntu
    install Boost into /usr/include, which is searched by default; Homebrew
    keeps it under its own prefix, which is not, so ask brew where it is.
    BOOST_INCLUDEDIR overrides both, for unusual installs.
    """
    from_env = os.environ.get("BOOST_INCLUDEDIR")
    if from_env:
        return [from_env]

    if os.path.isdir("/usr/include/boost"):
        return []

    brew = shutil.which("brew")
    if brew:
        result = subprocess.run(
            [brew, "--prefix", "boost"], capture_output=True, text=True, check=False
        )
        if result.returncode == 0:
            include_dir = os.path.join(result.stdout.strip(), "include")
            if os.path.isdir(include_dir):
                return [include_dir]

    return []


module = Extension(
    "libhankel",
    sources=[
        "src/py_interface.c",
        "src/form_factors.c",
        "src/hankel_transform.c",
        "src/hankel_DE_quadrature.c",
        "src/hankel_DHT.c",
        "src/hankel_QWE.c",
        "external_libs/DE-quadrature/intde.c",
        "external_libs/qwe_Chave/qwe_Chave.c",
        "external_libs/qwe_Key/qwe_Key.c",
        "external_libs/utils/tanhsinh.c",
        "src/utils/sasfit_integrate.c",
        "src/utils/sf_functions.c",
        "src/utils/validate_x.c",
        "src/utils/boost_bessel.cpp",
    ],
    include_dirs=[
        "include",
        "external_libs",
        ".",
        "src/utils",
        "external_libs/utils",
    ]
    + boost_include_dirs(),
)

setup(
    name="libhankel",
    version="0.1",
    ext_modules=[module],
)
