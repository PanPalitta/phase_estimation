#include <cmath>

#include "problem.h"

void Problem::modulo(double *can1) {
    for(int i = 0; i < num; ++i) {
        if(can1[i] < lower_bound[i]) {
            can1[i] = lower_bound[i];
        } else if (can1[i] > upper_bound[i]) {
            can1[i] = upper_bound[i];
        }
    }
}

void Problem::normalize(double *can1) {
    double norm = 0;
    for (int i = 0; i < num; ++i) {
        norm += abs(can1[i] * can1[i]);
    }
    norm = sqrt(norm);
    for (int i = 0; i < num; ++i) {
        can1[i] = can1[i] / norm;
    }
}
