from setuptools import Extension, setup

module = Extension(
    "libhankel",
    sources=[
        "src/py_interface.c",
        "src/form_factors.c",
        "src/hankel_transform.c",
        "src/hankel_DE_quadrature.c",
        "src/hankel_DHT.c",
        "src/hankel_FBT_Ogata.c",
        "src/hankel_QWE.c",
        "external_libs/DE-quadrature/intde.c",
        "external_libs/qwe_Chave/qwe_Chave.c",
        "external_libs/qwe_Key/qwe_Key.c",
        "external_libs/utils/tanhsinh.c",
        "src/utils/sasfit_integrate.c",
        "src/utils/sf_functions.c",
        "src/utils/boost_bessel.cpp",
    ],
    include_dirs=[
        "include",
        "external_libs",
        ".",
        "src/utils",
        "external_libs/utils",
    ],
)

setup(
    name="libhankel",
    version="0.1",
    ext_modules=[module],
)
