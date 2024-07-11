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


void write_hdf5_dataset(hid_t HDF_FileID, char *dataid, hid_t type,
                        hsize_t rank, hsize_t *dims, void *data) {
    hid_t dataspace_id = check_H5Screate_simple(rank, dims, NULL); 
    hid_t datatype_id = check_H5Tcopy(type);
    hid_t dataset_id = check_H5Dcreate(HDF_FileID, dataid, datatype_id, dataspace_id);
    check_H5Dwrite(dataset_id, type, data);
    check_H5Dclose(dataset_id);
    check_H5Tclose(datatype_id);
    check_H5Sclose(dataspace_id);
}

void write_hdf5_header(hid_t HDF_FileID, struct binary_output_header *bh) {
    hid_t attr_id, attrspace_id, attrtype_id;
    hsize_t dim[1];

    hid_t group_id = check_H5Gcreate(HDF_FileID, "/Header");

    // Magic
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "Magic", H5T_NATIVE_ULLONG, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_ULLONG, &bh->magic);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // Snap
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "Snap", H5T_NATIVE_LLONG, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_LLONG, &bh->snap);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // Chunk
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "Chunk", H5T_NATIVE_LLONG, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_LLONG, &bh->chunk);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // ScaleFactor
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "ScaleFactor", H5T_NATIVE_FLOAT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_FLOAT, &bh->scale);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // Omega_m
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "Omega_m", H5T_NATIVE_FLOAT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_FLOAT, &bh->Om);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // Omega_L
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "Omega_L", H5T_NATIVE_FLOAT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_FLOAT, &bh->Ol);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // H0
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "H0", H5T_NATIVE_FLOAT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_FLOAT, &bh->h0);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // Bounds
    attrspace_id = check_H5Screate(H5S_SIMPLE);
    dim[0] = 6;
    check_H5Sset_extent_simple(attrspace_id, 1, dim, NULL);
    attr_id = check_H5Acreate(group_id, "Bounds", H5T_NATIVE_FLOAT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_FLOAT, &bh->bounds[0]);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // NumberHalos
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "NumberHalos", H5T_NATIVE_LLONG, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_LLONG, &bh->num_halos);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // NumberParticles
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "NumberParticles", H5T_NATIVE_LLONG, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_LLONG, &bh->num_particles);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // BoxSize
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "BoxSize", H5T_NATIVE_FLOAT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_FLOAT, &bh->box_size);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // ParticleMass
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "ParticleMass", H5T_NATIVE_FLOAT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_FLOAT, &bh->particle_mass);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // ParticleType
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "ParticleType", H5T_NATIVE_LLONG, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_LLONG, &bh->particle_type);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // FormatRevision
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "FormatRevision", H5T_NATIVE_INT, attrspace_id);
    check_H5Awrite(attr_id, H5T_NATIVE_INT, &bh->format_revision);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);

    // Version
    attrtype_id = check_H5Tcopy(H5T_C_S1);
    check_H5Tset_size(attrtype_id, VERSION_MAX_SIZE);
    attrspace_id = check_H5Screate(H5S_SCALAR);
    attr_id = check_H5Acreate(group_id, "Version", attrtype_id, attrspace_id);
    check_H5Awrite(attr_id, attrtype_id, &bh->rockstar_version[0]);
    check_H5Aclose(attr_id);
    check_H5Sclose(attrspace_id);
    check_H5Tclose(attrtype_id);

    check_H5Gclose(group_id);
}

void output_hdf5(int64_t id_offset, int64_t snap, int64_t chunk,
                 float *bounds, int64_t output_particles) {
    float                       max[3] = {0}, min[3] = {0};
    char                        filename[1024];
    int64_t                     i, j, num_write, num_particles;
    int64_t                     offset, id;
    int64_t                    *ids, *p_start;
    int                        *to_write; 
    struct binary_output_header bheader;
    hid_t                       HDF_FileID, HDF_FileID_Part;
    hsize_t                     dims1[1], dims3[2], dims6[2], dims_part[1];
    float                      *buffer_float;
    int64_t                    *buffer_int, *buffer_id;


    memset(&bheader, 0, sizeof(struct binary_output_header));
    to_write = (int *) malloc(sizeof(int) * num_halos);
    ids = (int64_t *) malloc(sizeof(int64_t) * num_halos);
    p_start = (int64_t *) malloc(sizeof(int64_t) * num_halos);


    // Output Halos
    if (num_halos) {
        memcpy(min, halos[0].pos, sizeof(float) * 3);
        memcpy(max, halos[0].pos, sizeof(float) * 3);
    }
    for (i = 0, id = 0; i < num_halos; i++) {
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
    num_particles = bheader.num_particles;

    // set data dimensions and buffers
    dims1[0] = num_write;

    dims3[0] = num_write;
    dims3[1] = 3;

    dims6[0] = num_write;
    dims6[1] = 6;

    buffer_float = (float *) malloc(sizeof(float) * 6 * num_write);
    buffer_int = (int64_t *) malloc(sizeof(int64_t) * num_write);


    // Create file
    get_output_filename(filename, 1024, snap, chunk, "hdf5");
    HDF_FileID = check_H5Fcreate(filename, H5F_ACC_TRUNC);

    // ID
    write_hdf5_dataset(HDF_FileID, "/ID", H5T_NATIVE_LLONG, 1, dims1, ids);
    free(ids);

    // ParticleStart
    write_hdf5_dataset(HDF_FileID, "/ParticleStart", H5T_NATIVE_LLONG, 1, dims1, p_start);
    free(p_start);

    // Coordinates
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].pos[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Coordinates", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

    // Velocities
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].pos[j+3];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Velocities", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

    // CoreVelocities
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].corevel[j+3];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/CoreVelocities", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

    // BulkVelocities
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].bulkvel[j+3];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/BulkVelocities", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

    // Mvir
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].m;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Mvir", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // Rvir
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].r;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Rvir", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // ChildRadius
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].child_r;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/ChildRadius", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // VmaxRadius
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].vmax_r;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/VmaxRadius", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // Mgrav
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].mgrav;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Mgrav", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // Vmax
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].vmax;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Vmax", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

#ifdef OUTPUT_RVMAX
    // RVmax
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].rvmax;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/RVmax", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);
#endif

    // ScaleRadius
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].rs;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/ScaleRadius", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // ScaleRadiusKlypin
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].klypin_rs;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/ScaleRadiusKlypin", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // Vrms
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].vrms;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Vrms", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // AngularMometnum
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].J[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/AngularMomentum", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

    // Energy
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].energy;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Energy", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // Spin
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].spin;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Spin", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // M200b
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].alt_m[0];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/M200b", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // M200c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].alt_m[1];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/M200c", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // M500c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].alt_m[2];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/M500c", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // M2500c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].alt_m[3];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/M2500c", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // Xoff
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].Xoff;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Xoff", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // Voff
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].Voff;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Voff", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // IntermediateToMajorAxisRatio
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].b_to_a;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/IntermediateToMajorAxisRatio", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // MinorToMajorAxisRatio
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].c_to_a;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/MinorToMajorAxisRatio", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // IntermediateToMajorAxisRatio500c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].b_to_a2;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/IntermediateToMajorAxisRatio500c", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // MinorToMajorAxisRatio500c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].c_to_a2;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/MinorToMajorAxisRatio500c", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // MajorAxis
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].A[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/MajorAxis", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

    // MajorAxis500c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].A2[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/MajorAxis500c", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

#ifdef OUTPUT_INTERMEDIATE_AXIS
    // IntermediateAxis
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].A_I[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/IntermediateAxis", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);

    // IntermediateAxis500c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 3; j++)
            buffer_float[id * 3 + j] = halos[i].A_I2[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/IntermediateAxis500c", H5T_NATIVE_FLOAT, 2, dims3, buffer_float);
#endif

    // BullockSpin
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].bullock_spin;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/BullockSpin", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // KineticToPotentialRatio
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].kin_to_pot;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/KineticToPotentialRatio", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // KineticToPotentialRatio
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].kin_to_pot;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/KineticToPotentialRatio", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // M_pe_Behroozi
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].m_pe_b;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/M_pe_Behroozi", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // M_pe_Diemer
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].m_pe_d;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/M_pe_Diemer", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // HalfMassRadius
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].halfmass_radius;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/HalfMassRadius", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

#ifdef OUTPUT_INERTIA_TENSOR
    // InertiaTensor
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 6; j++)
            buffer_float[id * 6 + j] = halos[i].inertia_tensor[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/InertiaTensor", H5T_NATIVE_FLOAT, 2, dims6, buffer_float);

    // InertiaTensor500c
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        for (j = 0; j < 6; j++)
            buffer_float[id * 6 + j] = halos[i].inertia_tensor2[j];
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/InertiaTensor500c", H5T_NATIVE_FLOAT, 2, dims6, buffer_float);
#endif

    // NumberParticles
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_int[id] = halos[i].num_p;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/NumberParticles", H5T_NATIVE_LLONG, 1, dims1, buffer_int);

    // NumberChildParticles
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_int[id] = halos[i].num_child_particles;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/NumberChildParticles", H5T_NATIVE_LLONG, 1, dims1, buffer_int);

    // DescendantID
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_int[id] = halos[i].desc;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/DescendantID", H5T_NATIVE_LLONG, 1, dims1, buffer_int);

    // Flags
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_int[id] = halos[i].flags;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/Flags", H5T_NATIVE_LLONG, 1, dims1, buffer_int);

    // NumberCoreParticles
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_int[id] = halos[i].n_core;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/NumberCoreParticles", H5T_NATIVE_LLONG, 1, dims1, buffer_int);

    // MinPosError
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].min_pos_err;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/MinPosError", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // MinVelError
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].min_vel_err;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/MinVelError", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    // MinBulkVelError
    for (i = 0, id = 0; i < num_halos; i++) {
        if (!to_write[i])
            continue;
        buffer_float[id] = halos[i].min_bulkvel_err;
        id++;
    }
    write_hdf5_dataset(HDF_FileID, "/MinBulkVelError", H5T_NATIVE_FLOAT, 1, dims1, buffer_float);

    free(buffer_float);
    free(buffer_int);

    // Output header
    fill_binary_header(&bheader, snap, chunk);
    bheader.particle_type = PARTICLE_TYPE_IDS;
    if (bounds)
        memcpy(bheader.bounds, bounds, sizeof(float) * 6);
    else {
        memcpy(bheader.bounds, min, sizeof(float) * 3);
        memcpy(&(bheader.bounds[3]), max, sizeof(float) * 3);
    }
    write_hdf5_header(HDF_FileID, &bheader);

    // Output Particles
    if (output_particles) {
        get_output_filename(filename, 1024, snap, chunk, "id.hdf5");
        HDF_FileID_Part = check_H5Fcreate(filename, H5F_ACC_TRUNC);

        buffer_id = (int64_t *) malloc(sizeof(int64_t) * num_particles);
        dims_part[0] = num_particles;
        for (i = 0, offset = 0; i < num_halos; i++) {
            if (!to_write[i])
                continue;
            for (j = 0; j < halos[i].num_p; j++)
                buffer_id[offset+j] = p[halos[i].p_start + j].id;
            offset += halos[i].num_p;
        }
        write_hdf5_dataset(HDF_FileID_Part, "/ParticleID",
                           H5T_NATIVE_LLONG, 1, dims1, buffer_id);
        free(buffer_id);

        write_hdf5_header(HDF_FileID_Part, &bheader);

        check_H5Fclose(HDF_FileID_Part);
    }


    free(to_write);
    check_H5Fclose(HDF_FileID);
}

void load_hdf5_header(int64_t snap, int64_t chunk,
                      struct binary_output_header *bheader) {
    char  filename[1024];
    hid_t HDF_FileID;

    get_output_filename(filename, 1024, snap, chunk, "hdf5");
    HDF_FileID = check_H5Fopen(filename, H5F_ACC_RDONLY);

    check_H5Fclose(HDF_FileID);
    /*
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
    char  filename[1024];
    hid_t HDF_FileID;

    if (!coalesced)
        get_output_filename(filename, 1024, snap, chunk, "hdf5");
    else
        get_output_filename(filename, 1024, snap, chunk, "coalesced.hdf5");

    HDF_FileID = check_H5Fopen(filename, H5F_ACC_RDONLY);

    check_H5Fclose(HDF_FileID);
    /*
    int64_t i, j;

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
