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

        void write_contfit(double *fit, int tt) {
            int i;
            for(i = 0; i < num_fit; i++) {
                cont_fit[i] = fit[i] / double(tt);
                //cout<<cont_fit[i]<<",";
                }
            times = tt;
            //cout<<endl;
            }
        void write_bestfit(double *fit) {
            int i;
            for(i = 0; i < num_fit; i++) {
                best_fit[i] = best_fit[i] * double(best_times) + fit[i];
                }
            best_times += 1;
            for(i = 0; i < num_fit; i++) {
                best_fit[i] = best_fit[i] / double(best_times);
                }
            }
        void write_globalfit(double *fit) {
            int i;
            for(i = 0; i < num_fit; i++) {
                global_fit[i] = fit[i];
                }
            }

//    private:
        int num;
        int num_fit;
        typeT *can_best;
        typeT *contender;
        typeT *velocity;
        typeT *global_best;
        double *best_fit, *cont_fit, *global_fit; //should this be double or another type?
        int times, best_times, global_times; //number of samples used to calculate average best_fit

    };

template<typename typeT>
Candidate<typeT>::~Candidate() {
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
        throw out_of_range("numvar is not positive.");
        }
    if(fit_size <= 0) {
        throw invalid_argument("fit size is not positive");
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
void Candidate<typeT>::init_velocity() {
    velocity = new typeT[num];
    }


template<typename typeT>
void Candidate<typeT>::update_cont(typeT *input) {
    memcpy(contender, input, num * sizeof(typeT));
    }

template<typename typeT>
void Candidate<typeT>::update_vel(typeT *input) {
    int i;
    for (i = 0; i < num; i++) {
        velocity[i] = input[i];
        }
    }

template<typename typeT>
void Candidate<typeT>::update_best() {
    int i;
    for (i = 0; i < num; i++) {
        can_best[i] = contender[i];
        }
    for (i = 0; i < num_fit; i++) {
        best_fit[i] = cont_fit[i];
        }
    best_times = times;
    }

template<typename typeT>
void Candidate<typeT>::update_global(typeT *input) {
    int i;
    for (i = 0; i < num; i++) {
        global_best[i] = input[i];
        }
    }

template<typename typeT>
void Candidate<typeT>::put_to_global() {
    int i;
    for (i = 0; i < num; i++) {
        global_best[i] = can_best[i];
        }
    for (i = 0; i < num_fit; i++) {
        global_fit[i] = best_fit[i];
        }
    global_times = best_times;
    }

template<typename typeT>
void Candidate<typeT>::read_cont(typeT *output) {
    int i;
    for (i = 0; i < num; i++) {
        output[i] = contender[i];
        }
    }

template<typename typeT>
void Candidate<typeT>::read_vel(typeT *output) {
    int i;
    for (i = 0; i < num; i++) {
        output[i] = velocity[i];
        }
    }

template<typename typeT>
void Candidate<typeT>::read_best(typeT *output) {
    int i;
    for (i = 0; i < num; i++) {
        output[i] = can_best[i];
        }
    }

template<typename typeT>
void Candidate<typeT>::read_global(typeT *output) {
    int i;
    for (i = 0; i < num; i++) {
        output[i] = global_best[i];
        }
    }


#endif // CANDIDATE_H
