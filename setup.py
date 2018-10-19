from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name='rgg_ensemble_analysis',
    version='0.1',
    author='Matthew Garrod',
    packages=find_packages(),
    ext_modules=cythonize('rgg_ensemble_analysis/*.pyx',language="c++")
)
