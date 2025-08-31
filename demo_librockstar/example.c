#include <mpi.h>

/* Declare only the function needed from the library to avoid exposing
 * internal headers. */
extern "C" {
int rockstar_mpi_main(int argc, char **argv);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    rockstar_mpi_main(argc, argv);
    MPI_Finalize();
    return 0;
}
