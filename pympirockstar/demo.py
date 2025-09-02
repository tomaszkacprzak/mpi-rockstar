"""Example showing how to call MPI-Rockstar from Python."""

from mpi4py import MPI
from pympirockstar import run_config

if __name__ == "__main__":
    MPI.Init()
    try:
        # Replace the configuration path with one appropriate for your setup
        run_config("../examples/parallel_256.cfg")
    finally:
        MPI.Finalize()
