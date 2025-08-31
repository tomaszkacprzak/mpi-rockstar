/* Declare only the single entry point from the library to avoid exposing
 * internal headers or functions, but initialize MPI locally. */

#include <mpi.h>

extern "C" {
void rockstar_mpi_main(int argc, char **argv);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    rockstar_mpi_main(argc, argv);
    MPI_Finalize();
    return 0;
}
