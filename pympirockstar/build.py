import os
from setuptools import Extension
from Cython.Build import cythonize


def get_extensions():
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    lib_dir = os.path.join(root_dir, "src")
    lib_name = os.environ.get("ROCKSTAR_LIB", "mpi_rockstar")

    ext = Extension(
        name="pympirockstar._rockstar",
        sources=["_rockstar.pyx"],
        language="c++",
        include_dirs=[lib_dir],
        library_dirs=[lib_dir],
        libraries=[lib_name],
    )
    return cythonize([ext])
