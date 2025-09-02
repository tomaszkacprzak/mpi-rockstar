# cython: language_level=3
from libc.stdlib cimport malloc, free
from libc.string cimport strdup

cdef extern from "mpi_rockstar.h":
    void rockstar_mpi_main(int argc, char **argv) except +

def run(args):
    """Run MPI-Rockstar with a list of command line arguments.

    Parameters
    ----------
    args : list of str
        Arguments as would be passed on the command line. The program name
        should be the first element.
    """
    cdef int argc = len(args)
    cdef char **cargs = <char **>malloc((argc + 1) * sizeof(char *))
    if cargs == NULL:
        raise MemoryError()

    cdef Py_ssize_t i
    for i in range(argc):
        arg = args[i].encode('utf-8')
        cargs[i] = strdup(arg)
    cargs[argc] = NULL

    try:
        rockstar_mpi_main(argc, cargs)
    finally:
        for i in range(argc):
            free(cargs[i])
        free(cargs)
