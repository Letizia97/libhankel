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
            ],
            include_dirs=[
                "include"              
            ],
            libraries=["gsl", "gslcblas"], 
        )
    ],
)