#include <mpi.h>
#include "config.h"
#include "mpi_rockstar.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    do_config("parallel_256.cfg");
    mpi_main(0, NULL);
    MPI_Finalize();
    return 0;
}
