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
                "src/form_factors.c",        # ✅ now from top-level csrc
            ],
            include_dirs=[
                "src"               # ✅ tells compiler where add.h is
            ],
        )
    ],
)