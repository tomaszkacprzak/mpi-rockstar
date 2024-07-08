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


void write_hdf5_dataset(hid_t HDF_FileID, char *dataid, hid_t type, hsize_t rank, hsize_t *dims, void *data) {
    hid_t dataspace_id = H5Screate_simple(rank, dims, NULL); 
    hid_t datatype_id = H5Tcopy(type);
    hid_t dataset_id = H5Dcreate2(HDF_FileID, dataid, datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check_H5Dwrite(dataset_id, type, data);
    check_H5Dclose(dataset_id);
    check_H5Tclose(datatype_id);
    check_H5Sclose(dataspace_id);
}

void write_hdf5_header(struct binary_output_header bheader) {
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
    char                        filename[1024];
    int64_t                     i, j, num_write, id = 0;
    int64_t                    *ids, *p_start;
    int                        *to_write; 
    struct binary_output_header bheader;
    hid_t                       file_id;
    hsize_t                     dims1[1], dims2[2], dims3[2], dims4[2];
    float                      *buffer_float;
    int64_t                    *buffer_int;


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
        bheader.num_particles += halos[i].num_p;
        id++;
    }

    num_write = bheader.num_halos;

    // set data dimensions and buffers
    dims1[0] = num_write;

    dims2[0] = num_write;
    dims2[1] = 3;

    dims3[0] = num_write;
    dims3[1] = 4;

    dims4[0] = num_write;
    dims4[1] = 6;

    buffer_float = (float *) malloc(sizeof(float) * 6 * num_write);
    buffer_int = (int64_t *) malloc(sizeof(int64_t) * num_write);


    // Create file
    get_output_filename(filename, 1024, snap, chunk, "hdf5");
    file_id = check_H5Fcreate(filename, H5F_ACC_TRUNC);

    // ID
    write_hdf5_dataset(file_id, "/ID", H5T_NATIVE_LLONG, 1, dims1, ids);
    free(ids);

    // ParticleStart
    write_hdf5_dataset(file_id, "/ParticleStart", H5T_NATIVE_LLONG, 1, dims1, p_start);
    free(p_start);

    // Coordinates
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++){
            buffer_float[id * 3 + j] = halos[i].pos[j];
        }
        id++;
    }
    write_hdf5_dataset(file_id, "/Coordinates", H5T_NATIVE_FLOAT, 2, dims2, buffer_float);

    // Velocities
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++){
            buffer_float[id * 3 + j] = halos[i].pos[j+3];
        }
        id++;
    }
    write_hdf5_dataset(file_id, "/Velocities", H5T_NATIVE_FLOAT, 2, dims2, buffer_float);

    // CoreVelocities
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++){
            buffer_float[id * 3 + j] = halos[i].corevel[j+3];
        }
        id++;
    }
    write_hdf5_dataset(file_id, "/CoreVelocities", H5T_NATIVE_FLOAT, 2, dims2, buffer_float);

    // BulkVelocities
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++){
            buffer_float[id * 3 + j] = halos[i].bulkvel[j+3];
        }
        id++;
    }
    write_hdf5_dataset(file_id, "/BulkVelocities", H5T_NATIVE_FLOAT, 2, dims2, buffer_float);

    // Mvir
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].m;
        id++;
    }
    write_hdf5_dataset(file_id, "/Mvir", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // NumberParticles
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_int[id] = halos[i].num_p;
        id++;
    }
    write_hdf5_dataset(file_id, "/NumberParticles", H5T_NATIVE_LLONG, 1, dims1, buffer_int);

    free(buffer_float);
    free(buffer_int);

    /*
    // Output Particles
    if (output_particles) {
        for (i = 0; i < num_halos; i++) {
            if (!_should_print(halos + i, bounds))
                continue;
            for (j = 0; j < halos[i].num_p; j++)
                _append_to_buffer(&(p[halos[i].p_start + j].id),
                                  sizeof(int64_t), output);
        }
    }
    _clear_buffer(output);
    */

    free(to_write);

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
