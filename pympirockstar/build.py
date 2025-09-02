import os
from pathlib import Path
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize


class build_ext(_build_ext):
    """Custom build_ext to cythonize and link against Rockstar."""

    def finalize_options(self):
        super().finalize_options()
        # Transform any Cython sources into C++ code so we can link against
        # the Rockstar library.  cythonize will generate ``.cpp`` sources when
        # ``language="c++"`` is specified which avoids compilation errors from
        # C++ headers being included in ``.c`` files.
        self.extensions = cythonize(self.extensions, language="c++")

    def build_extensions(self):
        root_dir = Path(__file__).resolve().parents[1]
        lib_dir = root_dir / "src"
        lib_name = os.environ.get("ROCKSTAR_LIB", "mpi_rockstar")
        # Use the MPI C++ compiler if provided so the extension can link
        # against MPI-enabled Rockstar builds.
        mpicxx = os.environ.get("MPICC") or os.environ.get("MPICXX")
        if mpicxx and hasattr(self, "compiler"):
            self.compiler.set_executable("compiler_cxx", mpicxx)
            self.compiler.set_executable("linker_so", mpicxx)

        for ext in self.extensions:
            ext.include_dirs.append(str(lib_dir))
            ext.library_dirs.append(str(lib_dir))
            ext.libraries.append(lib_name)
            ext.language = "c++"
        super().build_extensions()
