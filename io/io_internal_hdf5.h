#ifndef _IO_INTERNAL_HDF5_H_
#define _IO_INTERNAL_HDF5_H_
#ifdef ENABLE_HDF5
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>
#include "../particle.h"

void set_buffer(void *buffer, const int *to_write,
                int64_t offset, int64_t stride, hsize_t type);
void load_buffer(void *buffer, struct halo *halos, int64_t to_read,
                 int64_t offset, int64_t stride, hsize_t type);
void write_hdf5_dataset(hid_t HDF_FileID, char *dataid, hid_t type,
                        hsize_t rank, hsize_t *dims, void *data);
void read_hdf5_dataset(hid_t HDF_FileID, char *dataid,
                       hid_t type, void *buffer);
void read_hdf5_halos(hid_t HDF_FileID, struct halo *halos,
                     struct binary_output_header *bh);
void write_hdf5_header(hid_t HDF_FileID, struct binary_output_header *bheader);
void read_hdf5_header(hid_t HDF_FileID, struct binary_output_header *bh, char *filename);
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
