#include <math.h>
#include <assert.h>
#include "config_vars.h"

int _check_bounds(float *pos_i, float *pos_f, float *bounds) {
    int64_t i;
    for (i = 0; i < 3; i++) {
        if (pos_i[i] > bounds[i + 3]) {
            if (bounds[i] < 0 && ROCKSTAR_PERIODIC) {
                pos_f[i] = pos_i[i] - ROCKSTAR_BOX_SIZE;
                if (pos_f[i] > bounds[i + 3] || pos_f[i] < bounds[i])
                    return 0;
            } else
                return 0;
        } else if (pos_i[i] < bounds[i]) {
            if (bounds[i + 3] > ROCKSTAR_BOX_SIZE && ROCKSTAR_PERIODIC) {
                pos_f[i] = pos_i[i] + ROCKSTAR_BOX_SIZE;
                if (pos_f[i] > bounds[i + 3] || pos_f[i] < bounds[i])
                    return 0;
            } else
                return 0;
        } else
            pos_f[i] = pos_i[i];
    }
    return 1;
}

void wrap_into_box(float *pos) {
    int64_t i;
    if (!ROCKSTAR_PERIODIC || !ROCKSTAR_BOX_SIZE)
        return;
    for (i = 0; i < 3; i++) {
        if (pos[i] > ROCKSTAR_BOX_SIZE)
            pos[i] -= ROCKSTAR_BOX_SIZE;
        else if (pos[i] < 0)
            pos[i] += ROCKSTAR_BOX_SIZE;
    }
}

int _check_bounds_raw(float *pos_i, float *bounds) {
    int64_t i;
    for (i = 0; i < 3; i++)
        if (pos_i[i] >= bounds[i + 3] || pos_i[i] < bounds[i])
            return 0;
    return 1;
}

// Note that b2 is extended by "overlap" in all directions
int bounds_overlap(float *b1, float *b2, float *b3, double overlap) {
    int64_t i;
    int     wrap = 0, max_wrap;
    float   min, max;
    for (i = 0; i < 3; i++) {
        b3[i] = min = b2[i] - overlap;
        b3[i + 3] = max = b2[i + 3] + overlap;
        wrap            = ((min < 0) && ROCKSTAR_PERIODIC) ? -1 : 0;
        max_wrap        = ((max > ROCKSTAR_BOX_SIZE) && ROCKSTAR_PERIODIC) ? 2 : 1;
        for (; wrap < max_wrap; wrap++) {
            if (((b1[i] + wrap * ROCKSTAR_BOX_SIZE) < max) &&
                ((b1[i + 3] + wrap * ROCKSTAR_BOX_SIZE) > min))
                break;
        }
        if (wrap == max_wrap)
            return 0;
    }
    return 1;
}

void bounds_union(float *b1, float *b2, float *result) {
    int64_t i = 0;
    for (i = 0; i < 3; i++)
        result[i] = (b1[i] < b2[i]) ? b1[i] : b2[i];
    for (; i < 6; i++)
        result[i] = (b1[i] > b2[i]) ? b1[i] : b2[i];
}
