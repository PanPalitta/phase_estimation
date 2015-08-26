#include <iostream>
#include <cstdlib>
#include <cstring>

#include "rng_simple.h"

RngSimple::RngSimple(int n_urandom_numbers, int n_grandom_numbers) {
    }

RngSimple::~RngSimple() {
    }

double RngSimple::next_grand(const double mean, const double dev) {
    return (double(rand()) / RAND_MAX - 0.5) * 2 * dev + mean;
    }

double RngSimple::next_urand() {
    return double(rand()) / RAND_MAX;
    }
