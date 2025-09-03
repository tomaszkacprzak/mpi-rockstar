"""Example showing how to call MPI-Rockstar from Python.

This script supports execution both as ``python -m pympirockstar.demo``
and directly via ``python demo.py`` from within the source tree.  The
latter requires a small ``sys.path`` tweak so the ``pympirockstar``
package can be located when it has not yet been installed.

The underlying MPI-Rockstar library performs its own MPI initialisation,
so this script merely imports ``mpi4py`` for its ``MPI`` namespace
without calling :func:`MPI.Init` or :func:`MPI.Finalize` explicitly.
"""

import os
import sys

try:  # pragma: no cover - import resolution differs by execution context
    from pympirockstar import run_config
except ImportError:  # Fallback when executed from the package directory
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    from pympirockstar import run_config

# MPI import without initialising
os.environ['MPI4PY_RC_INITIALIZE'] = '0'
from mpi4py import MPI 

if __name__ == "__main__":
    
    
    MPI.Init()
    # Replace the configuration path with one appropriate for your setup
    try:
        run_config("../examples/parallel_256.cfg")
    finally:
        MPI.Finalize()
