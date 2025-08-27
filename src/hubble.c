/* Courtesy of Matt Becker. */
#include <math.h>
#include <assert.h>
#include "config_vars.h"
#include "hubble.h"

/* See http://arxiv.org/pdf/astro-ph/0208512v1.pdf */
double _weff(double a) {
    if (a != 1.0)
        return ROCKSTAR_W0 + ROCKSTAR_WA - ROCKSTAR_WA * (a - 1.0) / log(a);
    else
        return ROCKSTAR_W0;
}

double hubble_scaling(double z) {
    double z1 = 1.0 + z;
    double a  = 1.0 / z1;
    return sqrt(ROCKSTAR_Om * (z1 * z1 * z1) + ROCKSTAR_Ol * pow(a, -3.0 * (1.0 + _weff(a))));
}
