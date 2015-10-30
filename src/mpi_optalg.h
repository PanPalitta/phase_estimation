#ifndef MPIOPTALG_H
#define MPIOPTALG_H
#include <mpi.h>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

#include "problem.h"
#include "candidate.h"

class OptAlg {
    public:
        OptAlg() {};
        OptAlg(Problem *problem_ptr) {
            prob = problem_ptr;
            num = prob->num;
            num_fit = prob->num_fit;
            }
        virtual ~OptAlg() {};

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

        double rand_Gaussian(double mean, double dev);
        void dev_gen(double *dev_array, double prev_dev, double new_dev, int cut_off);

        //Selecting solution
        double Final_select(int my_rank, int total_pop, int nb_proc, double *fit, double *solution, double *fitarray);
        double avg_Final_select(double* solution, int repeat, int my_rank, int total_pop, int nb_proc, double *soln_fit);

        int find_max(double *fit, int total_pop);

        //function for accept-reject criteria: should be moved to problem class at some point
        bool check_policy(double error, double sharp);//specific to phase estimation
        //calculating linear regression
        void linear_fit(int data_size, double *x, double *y, double *slope, double *intercept, double *mean_x);
        double error_interval(double *x, double *y, double mean_x, int data_size, double *SSres, double slope, double intercept);
        double error_update(int old_size, double *SSres, double *mean_x, double slope, double intercept, double *y, double *x);
        //calculating quantile
        double quantile(double p);
        inline double inv_erf(double x);
        inline int sgn(double x);

        bool success, policy_type;
        int num, num_fit;

    protected:
        int pop_size, T, t;
        Candidate *pop;
        bool goal;
    };

class DE : public OptAlg {
    public:
        DE() {};
        DE(Problem *problem_ptr): F(0.1), Cr(0.6) {
            this->prob = problem_ptr;
            this->num = this->prob->num;
            this->num_fit = prob->num_fit;
            }
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
        PSO(Problem *problem_ptr): w(0.8), phi1(0.6), phi2(1.0), v_max(0.2) {
            this->prob = problem_ptr;
            this->num = this->prob->num;
            this->num_fit = prob->num_fit;
            }
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
