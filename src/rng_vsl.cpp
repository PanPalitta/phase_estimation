#include <iostream>
#include <map>
#include <vector>
#include "mkl_vsl.h"

#include "rng.h"

#define BRNG VSL_BRNG_MCG31
#define METHOD VSL_RNG_METHOD_GAUSSIAN_ICDF

VSLStreamStatePtr stream;

using namespace std;

RngVsl::RngVsl(bool _gaussian, int _n_random_numbers, int seed, int rank):
    n_random_numbers(_n_random_numbers), RngVectorized(_gaussian) {
    vslNewStream(&stream, BRNG, seed + rank);
    random_numbers = new double[n_random_numbers];
    index_random_numbers = 0;
    if (gaussian) {
        vdRngGaussian(METHOD, stream, n_random_numbers, random_numbers, 0.0, 1.0);
        }
    else {
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_random_numbers, random_numbers, 0.0, 1.0);
        }
    }

double RngVsl::next_rand(const double mean, const double dev) {
    if (index_random_numbers >= n_random_numbers) {
        index_random_numbers = 0;
        if (gaussian) {
            vdRngGaussian(METHOD, stream, n_random_numbers, random_numbers, 0.0, 1.0);
            }
        else {
            vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_random_numbers, random_numbers, 0.0, 1.0);
            }
        }
    return random_numbers[index_random_numbers++] * dev + mean;
    }

RngVsl::~RngVsl() {
    delete[] random_numbers;
    vslDeleteStream(&stream);
    }
