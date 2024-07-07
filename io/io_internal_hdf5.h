#ifndef _IO_INTERNAL_HDF5_H_
#define _IO_INTERNAL_HDF5_H_
#ifdef ENABLE_HDF5
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>
#include "../particle.h"

void output_hdf5(int64_t id_offset, int64_t snap, int64_t chunk,
                 float *bounds, int64_t output_particles);
void load_hdf5_header(int64_t snap, int64_t chunk,
                      struct binary_output_header *bheader);
void load_hdf5_halos(int64_t snap, int64_t chunk,
                     struct binary_output_header *bheader,
                     struct halo **halos, int64_t **part_ids,
                     int64_t coalesced);

#endif /* ENABLE_HDF5 */
#endif /* _IO_INTERNAL_HDF5_H_ */
