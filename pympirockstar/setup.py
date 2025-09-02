from setuptools import setup, Extension
from Cython.Build import cythonize
import os

# Path to library directory
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
lib_dir = os.path.join(root_dir, 'src')

extensions = [
    Extension(
        name="pympirockstar._rockstar",
        sources=["_rockstar.pyx"],
        language="c++",
        include_dirs=[lib_dir],
        library_dirs=[lib_dir],
        libraries=["mpi_rockstar"],
    )
]

setup(
    name="pympirockstar",
    ext_modules=cythonize(extensions),
)
