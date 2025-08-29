#include <mpi.h>
#include <stdio.h>

extern "C" {
#include "config.h"
#include "mpi_rockstar.h"
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    char *cfg = argv[1];
    do_config(cfg);
    mpi_main(0, NULL);
    MPI_Finalize();
    return 0;
}
