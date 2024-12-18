#ifndef _MPI_MAIN_H_
#define _MPI_MAIN_H_

#ifdef __cplusplus
extern "C" {
#endif

void init_mpi(int argc, char *argv[]);
void mpi_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif /* _MPI_MAIN_H_ */
