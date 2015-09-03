#ifndef PSO_H
#define PSO_H

#include "mpi_optalg.h"

class PSO : public OptAlg {
    public:
        PSO() {};
        PSO(Problem *problem_ptr): w(0.8), phi1(0.6), phi2(1.0), v_max(0.2) {
            this->prob = problem_ptr;
            this->num = this->prob->num;
            }
        ~PSO() {};

        void put_to_best(int my_rank, int total_pop, int nb_proc);
        void combination(int my_rank, int total_pop, int nb_proc);
        void selection(int my_rank, int total_pop, int nb_proc);
        void write_param(double *param_array);
        void fit_to_global() {};
        void find_global(int my_rank, int total_pop, int nb_proc);

    private:
        double w, phi1, phi2, v_max;
    };
#endif
