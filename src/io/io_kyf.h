#ifndef _IO_KYF_H
#define _IO_KYF_H

#include <stdint.h>
#include "../particle.h"

void load_particles_kyf(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_KYF_H */
