#include <iostream>
#include <cstdlib>
#include <cstring>

#include "rng.h"

RngSimple::RngSimple(bool _gaussian, int n_random_numbers, int seed, int rank):
    RngBase(_gaussian) {
    srand(seed + rank);
    }

RngSimple::~RngSimple() {
    }

double RngSimple::next_rand(const double mean, const double dev) {
    if (gaussian) {
        return (double(rand()) / RAND_MAX - 0.5) * 2 * dev + mean;
        }
    else {
        return double(rand()) / RAND_MAX;
        }
    }
