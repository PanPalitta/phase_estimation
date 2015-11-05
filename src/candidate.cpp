#include "candidate.h"

Candidate::~Candidate() {
    if (is_candidate_initialized) {
        delete[] can_best;
        delete[] contender;
        delete[] global_best;
       
        delete[] best_fit;
        delete[] cont_fit;
        delete[] global_fit;
    }
  
    if (is_velocity_initialized) {
        delete[] velocity;
    }
}

void Candidate::init_can(int numvar, int fit_size) {
    if(numvar <= 0) {
        throw out_of_range("numvar must be positive.");
        }
    if(fit_size <= 0) {
        throw invalid_argument("num_fit must be positive.");
        }
    num = numvar;
    can_best = new double[num];
    contender = new double[num];
    global_best = new double[num];

    num_fit = fit_size;
    best_fit = new double[num_fit];
    cont_fit = new double[num_fit];
    global_fit = new double[num_fit];
    is_candidate_initialized = true;
}

void Candidate::init_velocity() { //This function can only be called after init_can. What can we do make sure it is safe to use?
    velocity = new double[num];
    is_velocity_initialized = true;
    }

void Candidate::update_cont(double *input) {
    memcpy(contender, input, num * sizeof(double));
    }

void Candidate::update_vel(double *input) {
    memcpy(velocity, input, num * sizeof(double));
    }

void Candidate::update_best() {
    memcpy(can_best, contender, num * sizeof(double));
    memcpy(best_fit, cont_fit, num_fit * sizeof(double));
    best_times = times;
    }

void Candidate::update_global(double *input) {
    memcpy(global_best, input, num * sizeof(double));
    }

void Candidate::put_to_global() {
    memcpy(global_best, can_best, num * sizeof(double));
    memcpy(global_fit, best_fit, num_fit * sizeof(double));
    global_times = best_times;
    }

void Candidate::read_cont(double *output) {
    memcpy(output, contender, num * sizeof(double));
    }

void Candidate::read_vel(double *output) {
    memcpy(output, velocity, num * sizeof(double));
    }

void Candidate::read_best(double *output) {
    memcpy(output, can_best, num * sizeof(double));
    }

void Candidate::read_global(double *output) {
    memcpy(output, global_best, num * sizeof(double));
    }

void Candidate::write_contfit(double *fit, int tt) {
    memcpy(cont_fit, fit, num_fit * sizeof(double));
    times = tt;
    }

void Candidate::write_bestfit(double *fit) {
    for(int i = 0; i < num_fit; i++) {
        best_fit[i] = best_fit[i] * double(best_times) + fit[i];
        }
    best_times += 1;
    for(int i = 0; i < num_fit; i++) {
        best_fit[i] = best_fit[i] / double(best_times);
        }
    }

void Candidate::write_globalfit(double *fit) {
    memcpy(global_fit, fit, num_fit * sizeof(double));
    }
