#ifndef _IO_HDF5_H_
#define _IO_HDF5_H_
#ifdef ENABLE_HDF5
#include <hdf5.h> /* HDF5 required */
#include <inttypes.h>

hid_t check_H5Fopen(char *filename, unsigned flags);

hid_t check_H5Gopen(hid_t HDF_FileID, char *gid, char *filename);

hid_t check_H5Dopen(hid_t HDF_GroupID, char *dataid, char *gid, char *filename);
hid_t check_H5Dopen2(hid_t HDF_GroupID, char *dataid);
hid_t check_H5Dget_space(hid_t HDF_DatasetID);
void  check_H5Dread(hid_t HDF_DatasetID, hid_t type, void *buffer, char *dataid,
                    char *gid, char *filename);
void  check_H5Dread2(hid_t HDF_DatasetID, hid_t type, void *buffer, char *dataid);

hid_t check_H5Aopen_name(hid_t HDF_GroupID, char *dataid, char *gid,
                         char *filename);
hid_t check_H5Aopen(hid_t HDF_GroupID, char *dataid, char *filename);
hid_t check_H5Aget_space(hid_t HDF_AttrID);
void  check_H5Aread(hid_t HDF_AttrID, hid_t type, void *buffer, char *dataid,
                    char *gid, char *filename);
void  check_H5Aread2(hid_t HDF_AttrID, hid_t type, void *buffer);

void    check_H5Sselect_all(hid_t HDF_DataspaceID);
int64_t check_H5Sget_simple_extent_ndims(hid_t HDF_DataspaceID);
void check_H5Sget_simple_extent_dims(hid_t HDF_DataspaceID, hsize_t *dimsize);
void check_H5Sset_extent_simple(hid_t HDF_DataspaceID, hsize_t rank,
                                hsize_t *dims, hsize_t *maxdims);

hid_t check_H5Fcreate(char *filename, unsigned flags);
void  check_H5Fclose(hid_t HDF_FileID);

hid_t check_H5Screate_simple(hsize_t rank, hsize_t *dims, hsize_t *maxdims);
hid_t check_H5Screate(H5S_class_t type);
void  check_H5Sclose(hid_t HDF_DataspaceID);

hid_t check_H5Tcopy(hid_t type);
void  check_H5Tset_size(hid_t type, size_t size);
void  check_H5Tclose(hid_t HDF_TypeID);

hid_t check_H5Dcreate(hid_t HDF_GroupID, char *dataid, hid_t type,
                      hid_t HDF_DataspaceID);
void  check_H5Dwrite(hid_t HDF_DatasetID, hid_t type, void *buffer);
void  check_H5Dclose(hid_t HDF_DatasetID);

hid_t check_H5Acreate(hid_t HDF_DatasetID, char *dataid, hid_t type,
                      hid_t HDF_DataspaceID);
void  check_H5Awrite(hid_t HDF_AttrID, hid_t type, void *buffer);
void  check_H5Aclose(hid_t HDF_AttrID);

hid_t check_H5Gcreate(hid_t HDF_FileID, char *groupid);
void  check_H5Gclose(hid_t HDF_GroupID);

#endif /* ENABLE_HDF5 */
#endif /* _IO_HDF5_H_ */
