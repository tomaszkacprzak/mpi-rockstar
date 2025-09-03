# pympirockstar

Python bindings for the MPI-Rockstar halo finder.

## Prerequisites

* A compiled MPI-Rockstar library in `../src`.
  Build the library from the repository root with:

  ```bash
  make libmpi_rockstar -C src               # standard build
  # or
  make libmpi_rockstar_hdf5 -C src          # with HDF5 support
  ```

## Installation

Install the wrapper with `pip`, using your MPI C/C++ compiler.  From the
repository root:

```bash
MPICC=mpicxx pip install -e ./pympirockstar
```

If you built the HDF5 variant, set `ROCKSTAR_LIB=mpi_rockstar_hdf5` before
running `pip` so the extension links against the correct library.

## Demo

After installation, run the example under MPI:

```bash
mpiexec -n 2 python -m pympirockstar.demo
```

The demo loads `examples/parallel_256.cfg`.  Adjust the path or configuration
as needed for your environment.
