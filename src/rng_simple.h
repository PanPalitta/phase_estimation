#ifndef RNG_SIMPLE_H
#define RNG_SIMPLE_H

#include "rng.h"

using namespace std;

class RngSimple: public Rng {
    public:
        RngSimple(int n_urandom_numbers, int n_grandom_numbers);
        ~RngSimple();
        double next_grand(const double mean, const double dev);
        double next_urand();
    };

#endif // RNG_SIMPLE_H
