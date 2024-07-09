#ifndef _IO_INTERNAL_HDF5_H_
#define _IO_INTERNAL_HDF5_H_
#ifdef ENABLE_HDF5
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>
#include "../particle.h"

void write_hdf5_dataset(hid_t HDF_FileID, char *dataid, hid_t type,
                        hsize_t rank, hsize_t *dims, void *data);
void write_hdf5_header(hid_t HDF_FileID, struct binary_output_header *bheader);
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
