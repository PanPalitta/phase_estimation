#include <iostream>
#include <cstdlib>
#include <cstring>

#include "rng.h"

RngSimple::RngSimple(int n_urandom_numbers, int n_grandom_numbers, int seed, 
                     int rank) {
    srand(seed + rank);
    }

RngSimple::~RngSimple() {
    }

double RngSimple::next_grand(const double mean, const double dev) {
    return (double(rand()) / RAND_MAX - 0.5) * 2 * dev + mean;
    }

double RngSimple::next_urand() {
    return double(rand()) / RAND_MAX;
    }
