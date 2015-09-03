#ifndef DE_H
#define DE_H

#include "mpi_optalg.h"

class DE : public OptAlg {
    public:
        DE() {};
        DE(Problem *problem_ptr): F(0.1), Cr(0.6) {
            this->prob = problem_ptr;
            this->num = this->prob->num;
            }
        ~DE() {};

        void put_to_best(int my_rank, int total_pop, int nb_proc);
        void combination(int my_rank, int total_pop, int nb_proc);
        void selection(int my_rank, int total_pop, int nb_proc);
        void write_param(double *param_array);
        void fit_to_global();
        void find_global(int my_rank, int total_pop, int nb_proc) {};

    private:
        double F, Cr;
    };
#endif // DE_H
