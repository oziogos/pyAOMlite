from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np


setup(
    name='pyAOMlite',
    ext_modules=cythonize([
        Extension(
            'src.anIres', sources=[
                'src/anIres.pyx',
                'src/anIres_functions.c'
            ]),
        Extension(
            'src.aom_overlap', sources=[
                'src/aom_overlap.pyx',
                'src/mulliken_functions.c'
            ],
            include_dirs=[np.get_include()]
        )
    ],
        compiler_directives={'language_level': '3'}
    ),
)
