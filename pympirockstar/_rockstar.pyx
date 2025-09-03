# cython: language_level=3
# distutils: language=c++
from libc.stdlib cimport malloc, free
from libc.string cimport strdup

cdef extern from "error.h":
    cdef cppclass rockstar_error(exception):
        int code
        const char *file
        int line

cdef extern from "mpi_rockstar.h":
    void rockstar_mpi_main(int argc, char **argv) except +

cdef class RockstarError(Exception):
    cdef public int code
    cdef public str file
    cdef public int line

    def __init__(self, int c, const char *f, int l):
        self.code = c
        self.file = f.decode('utf-8')
        self.line = l
        Exception.__init__(self, f"Rockstar error {c} at {self.file}:{self.line}")

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
        try:
            rockstar_mpi_main(argc, cargs)
        except rockstar_error as e:
            raise RockstarError(e.code, e.file, e.line)
    finally:
        for i in range(argc):
            free(cargs[i])
        free(cargs)
