#include <iostream>
#include <map>
#include <vector>
#include "mkl_vsl.h"

#include "rng_vsl.h"

#define BRNG VSL_BRNG_MCG31
#define METHOD VSL_RNG_METHOD_GAUSSIAN_ICDF
#define SEED 0

VSLStreamStatePtr stream;

using namespace std;


void vsl_init() {
    vslNewStream(&stream, BRNG, SEED);
}

void vsl_cache_init(double *urandom_numbers, int n_urandom_numbers,
                    double *grandom_numbers, int n_grandom_numbers) {
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_urandom_numbers, urandom_numbers, 0.0, 1.0);
    vdRngGaussian(METHOD, stream, n_grandom_numbers, grandom_numbers, 0.0, 1.0);
}

void vsl_close() {
    vslDeleteStream(&stream);
}
