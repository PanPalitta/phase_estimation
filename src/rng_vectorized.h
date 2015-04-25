#ifndef RNG_VECTORIZED_H
#define RNG_VECTORIZED_H

#include "rng.h"

using namespace std;

class RngVectorized: public Rng
{
public:
    RngVectorized(int n_urandom_numbers, int n_grandom_numbers);
    ~RngVectorized();
    double next_grand(const double mean, const double dev);
    double next_urand();

private:
    double *urandom_numbers;
    int n_urandom_numbers;
    int index_urandom_numbers;
    double *grandom_numbers;
    int n_grandom_numbers;
    int index_grandom_numbers;
};

#endif // RNG_VECTORIZED_H
