#ifdef ENABLE_HDF5
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "meta_io.h"
#include "io_hdf5.h"
#include "../config_vars.h"
#include "../rockstar.h"
#include "../groupies.h"
#include "../check_syscalls.h"
#include "../version.h"
#include "../halo.h"


void output_hdf5(int64_t id_offset, int64_t snap, int64_t chunk,
                 float *bounds, int64_t output_particles) {
    float                       max[3] = {0}, min[3] = {0};
    struct halo                 tmp;
    char                        filename[1024];
    int64_t                     i, j, id = 0;
    int64_t                     *ids, *p_start;
    struct binary_output_header bheader;
    hid_t                       file_id, attr;
    hid_t                       attribute_id, dataset_id, datatype_id, dataspace_id;
    herr_t                      status;
    hsize_t                     ndims;
    hsize_t                     dims1[1], dims2[2];
    float                       *buff;
    int64_t                     *bufi;


    memset(&bheader, 0, sizeof(struct binary_output_header));

    ids = (int64_t *) malloc(sizeof(int64_t) * num_halos);
    p_start = (int64_t *) malloc(sizeof(int64_t) * num_halos);


    // Output Halos
    if (num_halos) {
        memcpy(min, halos[0].pos, sizeof(float) * 3);
        memcpy(max, halos[0].pos, sizeof(float) * 3);
    }
    for (i = 0; i < num_halos; i++) {
        if (!_should_print(halos + i, bounds))
            continue;
        for (j = 0; j < 3; j++) {
            if (min[j] > halos[i].pos[j])
                min[j] = halos[i].pos[j];
            if (max[j] < halos[i].pos[j])
                max[j] = halos[i].pos[j];
        }
        ids[i] = id + id_offset;
        p_start[i] = bheader.num_particles;
        bheader.num_halos++;
        bheader.num_particles += tmp.num_p;
        id++;
    }

    get_output_filename(filename, 1024, snap, chunk, "hdf5");

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dims1[0] = num_halos;

    // id
    dataspace_id = H5Screate_simple(1, dims1, NULL); 
    datatype_id = H5Tcopy(H5T_NATIVE_LLONG);
    dataset_id = H5Dcreate2(file_id, "/id", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    free(ids);

    // p_start
    dataspace_id = H5Screate_simple(1, dims1, NULL); 
    datatype_id = H5Tcopy(H5T_NATIVE_LLONG);
    dataset_id = H5Dcreate2(file_id, "/p_start", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_start);
    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    free(p_start);

    // pos
    buff = (float *) malloc(sizeof(float) * 6 * num_halos);

    dims2[0] = num_halos;
    dims2[1] = 6;

    for (i = 0; i < num_halos; i++) {
        if (!_should_print(halos + i, bounds))
            continue;
        for (j = 0; j < 6; j++){
            buff[i * 6 + j] = halos[i].pos[j];
        }
    }

    dataspace_id = H5Screate_simple(2, dims2, NULL); 
    datatype_id = H5Tcopy(H5T_NATIVE_FLOAT);
    dataset_id = H5Dcreate2(file_id, "/pos", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    free(buff);

    // corevel
    buff = (float *) malloc(sizeof(float) * 3 * num_halos);

    dims2[0] = num_halos;
    dims2[1] = 3;

    for (i = 0; i < num_halos; i++) {
        if (!_should_print(halos + i, bounds))
            continue;
        for (j = 0; j < 3; j++){
            buff[i * 3 + j] = halos[i].corevel[j];
        }
    }

    dataspace_id = H5Screate_simple(2, dims2, NULL); 
    datatype_id = H5Tcopy(H5T_NATIVE_FLOAT);
    dataset_id = H5Dcreate2(file_id, "/corevel", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    free(buff);

    // bulkvel
    buff = (float *) malloc(sizeof(float) * 3 * num_halos);

    dims2[0] = num_halos;
    dims2[1] = 3;

    for (i = 0; i < num_halos; i++) {
        if (!_should_print(halos + i, bounds))
            continue;
        for (j = 0; j < 3; j++){
            buff[i * 3 + j] = halos[i].bulkvel[j];
        }
    }

    dataspace_id = H5Screate_simple(2, dims2, NULL); 
    datatype_id = H5Tcopy(H5T_NATIVE_FLOAT);
    dataset_id = H5Dcreate2(file_id, "/bulkvel", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    free(buff);

    // m
    buff = (float *) malloc(sizeof(float) * num_halos);

    for (i = 0; i < num_halos; i++) {
        if (!_should_print(halos + i, bounds))
            continue;
        buff[i] = halos[i].m;
    }

    dataspace_id = H5Screate_simple(1, dims1, NULL); 
    datatype_id = H5Tcopy(H5T_NATIVE_FLOAT);
    dataset_id = H5Dcreate2(file_id, "/m", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    free(buff);



    // Output header
    dataspace_id = H5Screate(H5S_NULL);
    datatype_id = H5Tcopy(H5T_NATIVE_INT);
    dataset_id = H5Dcreate2(file_id, "/Header", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    attribute_id = H5Screate(H5S_SCALAR);
    attr = H5Acreate2(dataset_id, "Omega_m", H5T_NATIVE_FLOAT, attribute_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, H5T_NATIVE_INT, &Om);

    status = H5Aclose(attr);
    status = H5Sclose(attribute_id);

    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);



    status = H5Fclose(file_id);
}

#endif /* ENABLE_HDF5 */
