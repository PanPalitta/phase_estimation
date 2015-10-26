#ifndef CANDIDATE_H
#define CANDIDATE_H
#include <stdexcept>

using namespace std;

template<typename typeT>
class Candidate {
    public:
        Candidate() {};
        ~Candidate();

        void init_can(int numvar, int fit_size);
        void init_velocity();

        void update_cont(typeT *input);
        void update_vel(typeT *input);
        void update_best();
        void update_global(typeT *input);
        void put_to_global();

        void read_cont(typeT *output);
        void read_vel(typeT *output);
        void read_best(typeT *output);
        void read_global(typeT *output);

        void write_contfit(double *fit, int tt);
        void write_bestfit(double *fit);
        void write_globalfit(double *fit);

        double read_contfit(int i) {
            return cont_fit[i];
            }
        double read_bestfit(int i) {
            return best_fit[i];
            }
        double read_globalfit(int i) {
            return global_fit[i];
            }
        int read_bestt() {
            return best_times;
            }

    private:

        int num, num_fit;
        double *best_fit, *cont_fit, *global_fit;
        int times, best_times, global_times;//number of samples used to calculate average best_fit

        //memory arrays
        typeT *can_best;
        typeT *contender;
        typeT *velocity;
        typeT *global_best;

    };

template<typename typeT>
Candidate<typeT>::~Candidate() {
    //delete[] must be commented out for unit testing in order for the framework to operate properly.
    delete[] can_best;
    delete[] contender;
    delete[] velocity;
    delete[] global_best;

    delete[] best_fit;
    delete[] cont_fit;
    delete[] global_fit;
    }

template<typename typeT>
void Candidate<typeT>::init_can(int numvar, int fit_size) {
    if(numvar <= 0) {
        throw out_of_range("numvar must be positive.");
        }
    if(fit_size <= 0) {
        throw invalid_argument("num_fit must be positive.");
        }
    num = numvar;
    can_best = new typeT[num];
    contender = new typeT[num];
    global_best = new typeT[num];

    num_fit = fit_size;
    best_fit = new double[num_fit];
    cont_fit = new double[num_fit];
    global_fit = new double[num_fit];
    }

template<typename typeT>
void Candidate<typeT>::init_velocity() { //This function can only be called after init_can. What can we do make sure it is safe to use?
    velocity = new typeT[num];
    }

template<typename typeT>
void Candidate<typeT>::update_cont(typeT *input) {
    memcpy(contender, input, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::update_vel(typeT *input) {
    memcpy(velocity, input, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::update_best() {
    memcpy(can_best, contender, num * sizeof(typeT));
    memcpy(best_fit, cont_fit, num_fit * sizeof(double));
    best_times = times;
    }

template<typename typeT>
void Candidate<typeT>::update_global(typeT *input) {
    memcpy(global_best, input, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::put_to_global() {
    memcpy(global_best, can_best, num * sizeof(typeT));
    memcpy(global_fit, best_fit, num_fit * sizeof(double));
    global_times = best_times;
    }

template<typename typeT>
void Candidate<typeT>::read_cont(typeT *output) {
    memcpy(output, contender, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::read_vel(typeT *output) {
    memcpy(output, velocity, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::read_best(typeT *output) {
    memcpy(output, can_best, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::read_global(typeT *output) {
    memcpy(output, global_best, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::write_contfit(double *fit, int tt) {
    for(int i = 0; i < num_fit; i++) {
        cont_fit[i] = fit[i] / double(tt);
        }
    times = tt;
    }

template<typename typeT>
void Candidate<typeT>::write_bestfit(double *fit) {
    for(int i = 0; i < num_fit; i++) {
        best_fit[i] = best_fit[i] * double(best_times) + fit[i];
        }
    best_times += 1;
    for(int i = 0; i < num_fit; i++) {
        best_fit[i] = best_fit[i] / double(best_times);
        }
    }

template<typename typeT>
void Candidate<typeT>::write_globalfit(double *fit) {
    memcpy(global_fit, fit, num_fit * sizeof(double));
    }

#endif // CANDIDATE_H
