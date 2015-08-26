#ifndef RNG_VSL_H
#define RNG_VSL_H

#include "rng_vectorized.h"

using namespace std;

class RngVsl: public RngVectorized {
    public:
        RngVsl(int n_urandom_numbers, int n_grandom_numbers);
        ~RngVsl();
        double next_grand(const double mean, const double dev);
        double next_urand();

    private:
        double *urandom_numbers;
        int n_urandom_numbers;
        int index_urandom_numbers;
        double *grandom_numbers;
        int n_grandom_numbers;
        int index_grandom_numbers;
    }

#endif
