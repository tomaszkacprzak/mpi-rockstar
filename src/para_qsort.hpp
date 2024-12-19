#ifndef _PARA_QSORT_H_
#define _PARA_QSORT_H_

#include <algorithm>
#include <stdint.h>
#include <mpi.h>
#include <omp.h>

MPI_Datatype create_mpi_particle_type();
void         timed_output(const char *, ...);

template <typename T> T      *reallocate(T *, size_t);
template <typename T> T      *allocate(size_t);
template <typename Size> void clear_counts(Size[], size_t);
template <typename T>
int64_t exchange_data(const T[], int[], T *&, const int, MPI_Datatype);

#define para_qsort_THRESHOLD (5000)

template <typename T> static inline T median(const T a, const T b, const T c) {
    return std::max(std::min(a, b), std::min(std::max(a, b), c));
}

/* ### qsort_partitioning ### */
template <typename T, typename U>
static inline void qsort_partitioning_key(U a[], T key[], int &left, int &right,
                                          T pivot) {
    int ii = left;
    int jj = right;

    /* partitioning */
    while (1) {
        while (key[ii] < pivot)
            ii++;
        while (pivot < key[jj])
            jj--;
        if (ii >= jj)
            break;
        std::swap(a[ii], a[jj]);
        std::swap(key[ii], key[jj]);
        ii++;
        jj--;
    }

    left  = ii;
    right = jj;
}

/* ### single_qsort ### */
template <typename T, typename U>
static void single_qsort_key(U a[], T key[], int left, int right) {
    if (left < right) {
        int     ii = left, jj = right;
        const T pivot = median(key[ii], key[jj], key[(ii + jj) / 2]);

        qsort_partitioning_key(a, key, ii, jj, pivot);
        single_qsort_key(a, key, left, ii - 1);
        single_qsort_key(a, key, jj + 1, right);
    }
}

/* ### para_qsort_internal ### */
template <typename T, typename U>
static void para_qsort_internal_key(U a[], T key[], int left, int right) {
    int length = right - left;
    if (length < para_qsort_THRESHOLD) {
        single_qsort_key(a, key, left, right);
        return;
    }

    int     ii = left, jj = right;
    const T pivot = median(key[ii], key[jj], key[(ii + jj) / 2]);

    qsort_partitioning_key(a, key, ii, jj, pivot);

#pragma omp task
    para_qsort_internal_key(a, key, left, jj);
#pragma omp task
    para_qsort_internal_key(a, key, ii, right);
}

/* ### para_qsort ### */
template <typename T, typename U>
static void para_qsort_key(U a[], T key[], int left, int right) {
    if (omp_in_parallel() != 0) {
        single_qsort_key(a, key, left, right);
        return;
    }

#pragma omp parallel
    {
#pragma omp single nowait
        {
            para_qsort_internal_key(a, key, left, right);
        }
    }
}

/* ### para_qsort_axis ### */
template <typename T, typename U>
static void para_qsort_axis(U p[], int left, int right, const int axis) {
    if (omp_in_parallel() != 0) {
        single_qsort_axis(p, left, right, axis);
        return;
    }

#pragma omp parallel
    {
#pragma omp single nowait
        {
            para_qsort_internal_axis(p, left, right, axis);
        }
    }
}

#endif /* _PARA_QSORT_H_ */
