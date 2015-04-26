#ifndef RNG_GPU_H
#define RNG_GPU_H
#include <curand.h>
#include <mpi.h>

#include "rng_vectorized.h"

void setDevice(int commRank, int commSize);

class RngGpu: public RngVectorized
{
public:
    RngGpu(int n_urandom_numbers, int n_grandom_numbers);
    ~RngGpu();
    double next_grand(const double mean, const double dev);
    double next_urand();

private:
    double *urandom_numbers;
    int n_urandom_numbers;
    int index_urandom_numbers;
    double *grandom_numbers;
    int n_grandom_numbers;
    int index_grandom_numbers;

    double *dev_urandom_numbers;
    double *dev_grandom_numbers;
    curandGenerator_t gen;
};

#endif
