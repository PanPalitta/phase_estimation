#ifndef MPIOPTALG_H
#define MPIOPTALG_H
#include <mpi.h>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

#include "problem.h"
#include "candidate.h"
#include "rng.h"
//#include "aux_functions.h"

class OptAlg {
    public:
        OptAlg() {};
        OptAlg(Problem *problem_ptr, Rng *gaussian_rng): gaussian_rng(gaussian_rng) {
            prob = problem_ptr;
            num = prob->num;
            num_fit = prob->num_fit;
            }
        virtual ~OptAlg() {
            delete[] pop;
            };

        virtual void put_to_best(int my_rank, int total_pop, int nb_proc) {};
        virtual void combination(int my_rank, int total_pop, int nb_proc) {};
        virtual void selection(int my_rank, int total_pop, int nb_proc) {};
        virtual void write_param(double *param_array) {};
        virtual void read_param(double *param_array) {};
        virtual void fit_to_global() {};
        virtual void find_global(int my_rank, int total_pop, int nb_proc) {};


        Problem* prob;

        //functions and variables that are not algorithm specific
        void Init_population(int psize);
        void Init_previous(double prev_dev, double new_dev, int psize, double *prev_soln);

        void Cont_fitness(int p);
        void Best_fitness(int p);
        void update_popfit();

        void set_success(int iter, bool goal);
        bool check_success(int t, int D, double fit, double slope, double intercept);

        //Selecting solution
        double Final_select(int my_rank, int total_pop, int nb_proc, double *fit, double *solution, double *fitarray);
        double avg_Final_select(double* solution, int repeat, int my_rank, int total_pop, int nb_proc, double *soln_fit);

        void dev_gen(double *dev_array, double prev_dev, double new_dev, int cut_off);
        int find_max(double *fit, int total_pop);

        bool success, policy_type;
        int num, num_fit;

    protected:
        Rng *gaussian_rng;
        int pop_size, T, t;
        Candidate *pop;
        bool goal;
    };

class DE : public OptAlg {
    public:
        DE() {};
        DE(Problem *problem_ptr, Rng *gaussian_rng): OptAlg(problem_ptr, gaussian_rng), F(0.1), Cr(0.6) {};
        ~DE() {};

        void put_to_best(int my_rank, int total_pop, int nb_proc);
        void combination(int my_rank, int total_pop, int nb_proc);
        void selection(int my_rank, int total_pop, int nb_proc);
        void fit_to_global();
        void find_global(int my_rank, int total_pop, int nb_proc) {};
        void write_param(double *param_array);
        void read_param(double *param_array);

    private:
        double F, Cr;

        void family_gen(int* fam, int p, int fam_size, int total_pop);

    };

class PSO : public OptAlg {
    public:
        PSO() {};
        PSO(Problem *problem_ptr, Rng *gaussian_rng): OptAlg(problem_ptr, gaussian_rng), w(0.8), phi1(0.6), phi2(1.0), v_max(0.2) {};
        ~PSO() {};

        void put_to_best(int my_rank, int total_pop, int nb_proc);
        void combination(int my_rank, int total_pop, int nb_proc);
        void selection(int my_rank, int total_pop, int nb_proc);
        void write_param(double *param_array);
        void read_param(double *param_array);
        void fit_to_global() {};
        void find_global(int my_rank, int total_pop, int nb_proc);

    private:
        inline void find_index(int *prev, int *forw, int p, int total_pop);
        inline int find_fitness(int prev, double prev_fit, int forw, double forw_fit, int p, double fit);
        double w, phi1, phi2, v_max;
    };

#endif // OPTALG_H
