#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "io_pkdgrav3lcp.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../universal_constants.h"

#define PKDGRAV3LCP_PARTICLE_BYTES 40

void load_particles_pkdgrav3lcp(char *filename, struct particle **p, int64_t *num_p) {
    FILE   *input;
    int64_t file_particles, i, j;
    int64_t id;
    float   pos[3], vel[3];
    char    trash[8];
    double  vel_unit;

    input = check_fopen(filename, "rb");
    fseek(input, 0, SEEK_END);
    file_particles = ftell(input) / PKDGRAV3LCP_PARTICLE_BYTES;
    if (file_particles <= 0) {
        fclose(input);
        return;
    }
    fseek(input, 0, SEEK_SET);

    *p = check_realloc(*p, ((*num_p) + file_particles) * sizeof(struct particle),
                       "Adding pkdgrav3lcp particles.");

    vel_unit = ((300.0 * ROCKSTAR_BOX_SIZE) / sqrt(8.0 * M_PI)) * PKDGRAV3LCP_VEL_SCALE;

    for (i = 0; i < file_particles; i++) {
        check_fread(&id, sizeof(int64_t), 1, input);
        check_fread(pos, sizeof(float), 3, input);
        check_fread(vel, sizeof(float), 3, input);
        check_fread(trash, 1, 8, input); /* potential + padding */

        (*p)[(*num_p) + i].id = id;
        for (j = 0; j < 3; j++) {
            (*p)[(*num_p) + i].pos[j] = pos[j] * ROCKSTAR_BOX_SIZE * PKDGRAV3LCP_POS_SCALE + PKDGRAV3LCP_POS_SHIFT;
            (*p)[(*num_p) + i].pos[j + 3] = vel[j] * vel_unit;
        }
    }

    *num_p += file_particles;
    fclose(input);

    if (!ROCKSTAR_PARTICLE_MASS || ROCKSTAR_RESCALE_PARTICLE_MASS) {
        ROCKSTAR_PARTICLE_MASS = pow(ROCKSTAR_BOX_SIZE, 3) * 2.77536627208e11;
    }
    ROCKSTAR_AVG_PARTICLE_SPACING =
        cbrt(ROCKSTAR_PARTICLE_MASS / (ROCKSTAR_Om * CRITICAL_DENSITY));
    ROCKSTAR_SCALE_NOW = 1.0;
}

