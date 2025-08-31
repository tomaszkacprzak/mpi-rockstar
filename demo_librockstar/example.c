/* Declare only the single entry point from the library to avoid exposing
 * internal headers or functions. */
extern "C" {
void rockstar_mpi_main(int argc, char **argv);
}

int main(int argc, char **argv) {
    rockstar_mpi_main(argc, argv);
    return 0;
}
