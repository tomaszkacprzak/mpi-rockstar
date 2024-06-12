#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "meta_io.h"
#include "../config_vars.h"
#include "../rockstar.h"
#include "../groupies.h"
#include "../check_syscalls.h"
#include "../version.h"
#include "../halo.h"


void output_hdf5(int64_t id_offset, int64_t snap, int64_t chunk,
                 float *bounds, int64_t output_particles) {
    struct halo                 tmp;
    char                        filename[1024];
    int64_t                     i, j, id = 0;
    struct binary_output_header bheader;


    get_output_filename(filename, 1024, snap, chunk, "hdf5");

    hid_t HDF_FileID = check_H5Fopen(filename, H5F_ACC);

    // Output Halos
    for (i = 0; i < num_halos; i++) {
        tmp         = halos[i];
        tmp.id      = id + id_offset;
        tmp.p_start = bheader.num_particles;
        _append_to_buffer(&tmp, sizeof(struct halo), output);
        bheader.num_halos++;
        bheader.num_particles += tmp.num_p;
        id++;
    }


    // Output header
    fill_binary_header(&bheader, snap, chunk);
    bheader.particle_type = PARTICLE_TYPE_IDS;
    if (bounds)
        memcpy(bheader.bounds, bounds, sizeof(float) * 6);
    else {
        memcpy(bheader.bounds, min, sizeof(float) * 3);
        memcpy(&(bheader.bounds[3]), max, sizeof(float) * 3);
    }
    rewind(output);
    check_fwrite(&bheader, sizeof(struct binary_output_header), 1, output);
    fclose(output);
}
