from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
    name='pyAOMlite',
    ext_modules=cythonize([
        Extension(
            'src.mulliken', sources=[
                'src/mulliken.pyx',
                'src/mulliken_functions.c'
            ]),
        Extension(
            'src.anIres', sources=[
                'src/anIres.pyx',
                'src/anIres_functions.c'
            ])
    ]),
    compiler_directives={'language_level': "3"}
)
