#include <mpi.h>

extern "C" {
#include "config.h"
#include "mpi_rockstar.h"
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    char cfg[] = "rockstar_lightcone.config";
    do_config(cfg);
    mpi_main(0, NULL);
    MPI_Finalize();
    return 0;
}
