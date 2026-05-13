from setuptools import setup, Extension

setup(
    name="libhankel",
    version="0.1.0",
    packages=["libhankel"],
    ext_modules=[
        Extension(
            "libhankel._core",
            sources=[
                "libhankel/_core.c",
                "src/form_factors.c",     
                "src/hankel_transform.c",    
                "src/hankel_DE_quadrature.c",
                "src/hankel_DHT.c",
                "src/hankel_FBT_Ogata.c",
                "src/hankel_QWE.c",
                "external_libs/DE-quadrature/intde.c",
                "external_libs/qwe_Chave/qwe_Chave.c",
                "external_libs/qwe_Key/qwe_Key.c",
                "external_libs/Ogata/FTB.cpp",
                "external_libs/utils/tanhsinh.c",
                "src/utils/sasfit_integrate.c",
                "src/wrappers/wrapper_FBT.cpp"
                
            ],
            include_dirs=[
                "include",
                "external_libs",
                ".",
                "src/utils",
                "external_libs/Ogata",
                "external_libs/utils"                 
            ],
            libraries=["gsl", "gslcblas"], 
        )
    ],
)