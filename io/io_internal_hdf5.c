#ifdef ENABLE_HDF5
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <hdf5.h> /* HDF5 required */
#include "meta_io.h"
#include "io_internal_hdf5.h"
#include "io_internal.h"
#include "io_hdf5.h"
#include "../config_vars.h"
#include "../rockstar.h"
#include "../groupies.h"
#include "../check_syscalls.h"
#include "../version.h"
#include "../halo.h"

void write_hdf5_bheader(struct binary_output_header bheader) {
    /*
    dataspace_id = H5Screate(H5S_NULL);
    datatype_id = H5Tcopy(H5T_NATIVE_INT);
    dataset_id = H5Dcreate2(file_id, "/Header", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    attribute_id = H5Screate(H5S_SCALAR);
    attr = H5Acreate2(dataset_id, "Omega_m", H5T_NATIVE_DOUBLE, attribute_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &Om);

    status = H5Aclose(attr);
    status = H5Sclose(attribute_id);

    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    */
}

void output_hdf5(int64_t id_offset, int64_t snap, int64_t chunk,
                 float *bounds, int64_t output_particles) {
    float                       max[3] = {0}, min[3] = {0};
    struct halo                 tmp;
    char                        filename[1024];
    int64_t                     i, j, num_write, id = 0;
    int64_t                    *ids, *p_start;
    int                        *to_write; 
    struct binary_output_header bheader;
    hid_t                       file_id, attr;
    hid_t                       attribute_id, dataset_id, datatype_id, dataspace_id;
    herr_t                      status;
    hsize_t                     ndims;
    hsize_t                     dims1[1], dims2[2];
    float                       *buff;
    int64_t                     *bufi;


    memset(&bheader, 0, sizeof(struct binary_output_header));
    to_write = (int *) malloc(sizeof(int) * num_halos);

    ids = (int64_t *) malloc(sizeof(int64_t) * num_halos);
    p_start = (int64_t *) malloc(sizeof(int64_t) * num_halos);


    // Output Halos
    if (num_halos) {
        memcpy(min, halos[0].pos, sizeof(float) * 3);
        memcpy(max, halos[0].pos, sizeof(float) * 3);
    }
    for (i = 0; i < num_halos; i++) {
        if (!_should_print(halos + i, bounds)) {
            to_write[i] = 0;
            continue;
        }

        to_write[i] = 1;

        for (j = 0; j < 3; j++) {
            if (min[j] > halos[i].pos[j])
                min[j] = halos[i].pos[j];
            if (max[j] < halos[i].pos[j])
                max[j] = halos[i].pos[j];
        }

        ids[id] = id + id_offset;
        p_start[id] = bheader.num_particles;
        bheader.num_halos++;
        bheader.num_particles += tmp.num_p;
        id++;
    }

    num_write = bheader.num_halos;

    get_output_filename(filename, 1024, snap, chunk, "hdf5");
    file_id = check_H5Fcreate(filename, H5F_ACC_TRUNC);

    dims1[0] = num_write;

    // ID
    dataspace_id = check_H5Screate_simple(1, dims1, NULL); 
    datatype_id = check_H5Tcopy(H5T_NATIVE_LLONG);
    dataset_id = check_H5Dcreate(file_id, "/ID", datatype_id, dataspace_id);
    check_H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
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

/*
    // pos
    buff = (float *) malloc(sizeof(float) * 6 * num_write);

    dims2[0] = num_write;
    dims2[1] = 6;

    for (i = 0; i < num_halos; i++) {
        if (!to_write[i])
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
    buff = (float *) malloc(sizeof(float) * 3 * num_write);

    dims2[0] = num_write;
    dims2[1] = 3;

    for (i = 0; i < num_halos; i++) {
        if (!to_write[i])
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
    buff = (float *) malloc(sizeof(float) * 3 * num_write);

    dims2[0] = num_write;
    dims2[1] = 3;

    for (i = 0; i < num_halos; i++) {
        if (!to_write[i])
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

    // Mvir
    buff = (float *) malloc(sizeof(float) * num_halos);

    for (i = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buff[i] = halos[i].m;
    }

    dataspace_id = H5Screate_simple(1, dims1, NULL); 
    datatype_id = H5Tcopy(H5T_NATIVE_FLOAT);
    dataset_id = H5Dcreate2(file_id, "/Mvir", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buff);
    status = H5Dclose(dataset_id);
    status = H5Tclose(datatype_id);
    status = H5Sclose(dataspace_id);
    free(buff);
*/


    // Output header
    fill_binary_header(&bheader, snap, chunk);
    bheader.particle_type = PARTICLE_TYPE_IDS;
    if (bounds)
        memcpy(bheader.bounds, bounds, sizeof(float) * 6);
    else {
        memcpy(bheader.bounds, min, sizeof(float) * 3);
        memcpy(&(bheader.bounds[3]), max, sizeof(float) * 3);
    }


    check_H5Fclose(file_id);
}

void load_hdf5_header(int64_t snap, int64_t chunk,
                      struct binary_output_header *bheader) {
    /*
    char  buffer[1024];
    FILE *input;
    get_output_filename(buffer, 1024, snap, chunk, "bin");
    input = check_fopen(buffer, "rb");
    check_fread(bheader, sizeof(struct binary_output_header), 1, input);
    assert(bheader->magic == ROCKSTAR_MAGIC);
    fclose(input);
    */
}

void load_hdf5_halos(int64_t snap, int64_t chunk,
                     struct binary_output_header *bheader,
                     struct halo **halos, int64_t **part_ids,
                     int64_t coalesced) {
    /*
    char    buffer[1024];
    FILE   *input;
    int64_t i, j;

    if (!coalesced)
        get_output_filename(buffer, 1024, snap, chunk, "bin");
    else
        get_output_filename(buffer, 1024, snap, chunk, "coalesced.bin");
    input = check_fopen(buffer, "rb");
    check_fread(bheader, sizeof(struct binary_output_header), 1, input);
    assert(bheader->magic == ROCKSTAR_MAGIC);
    assert(bheader->num_halos >= 0);
    assert(bheader->num_particles >= 0);
    check_realloc_s(*halos, sizeof(struct halo), bheader->num_halos);
    check_realloc_s(*part_ids, sizeof(int64_t), bheader->num_particles);
    i = check_fread(*halos, sizeof(struct halo), bheader->num_halos, input);
    j = check_fread(*part_ids, sizeof(int64_t), bheader->num_particles, input);
    if (i != bheader->num_halos || j != bheader->num_particles) {
        fprintf(stderr, "[Error] Truncated input file %s!\n", buffer);
        exit(1);
    }
    fclose(input);

    read_binary_header_config(bheader);
    */
}


#endif /* ENABLE_HDF5 */
