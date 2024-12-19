#ifndef FOF_H
#define FOF_H

#include <stdint.h>
#include "particle.h"

struct fof {
    int64_t          num_p;
    struct particle *particles;
};

struct smallfof {
    int64_t root;
};

struct FOFInfo {
    struct fof *fofs;
    int64_t     num_fofs, num_alloced_fofs;
    int64_t     num_boundary_fofs;

    struct smallfof *smallfofs;
    int64_t          num_smallfofs, num_alloced_smallfofs;

    struct particle *root_p;
    int64_t         *particle_smallfofs;
    int64_t          num_particles, num_alloced_particles;
};

void init_fofinfo(struct FOFInfo *fofinfo);
void free_fofinfo(struct FOFInfo *fofinfo);

void    init_particle_smallfofs(int64_t num_p, struct particle *particles,
                                struct FOFInfo *fofinfo);
void    link_particle_to_fof(struct particle *p, int64_t n,
                             struct particle **links, struct FOFInfo *fofinfo);
void    link_fof_to_fof(struct particle *p, int64_t n, struct particle **links,
                        const struct FOFInfo *fofinfo);
int64_t tag_boundary_particle(struct particle *p, struct FOFInfo *fofinfo);

void        build_fullfofs(struct FOFInfo *fofinfo);
struct fof *return_fullfofs(int64_t *num_f, int64_t *num_bf,
                            struct FOFInfo *fofinfo);
void copy_fullfofs(struct fof **base, int64_t *num_f, int64_t *num_alloced_f,
                   struct FOFInfo *fofinfo);

void partition_sort_particles(int64_t min, int64_t max,
                              struct particle *particles, int64_t *assignments);

#endif /* FOF_H */
