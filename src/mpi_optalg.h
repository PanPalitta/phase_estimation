#ifndef MPIOPTALG_H
#define MPIOPTALG_H
#include <mpi.h>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

#include "problem.h"
#include "candidate.h"
#include "rng.h"
#include "aux_functions.h"

/*! \brief OptAlg class contains the functions that can be used to contruct an optimization algorithm in the main function.
    The object needs the pointer to the problem object and the RNG object in order to instantiate.
*/

class OptAlg {
    public:
        OptAlg() {};
        OptAlg(Problem *problem_ptr, Rng *gaussian_rng, int pop_size): gaussian_rng(gaussian_rng) {
            prob = problem_ptr; //Store the pointer to the problem object.
            num = prob->num; //Store the number of variables.
            num_fit = prob->num_fit; //Store the number of fitness values.

            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
            total_pop = pop_size; //Store the population size.
            }
        virtual ~OptAlg() {
            delete[] pop;
            };

        /*List of functions which are specific to the optimization functions*/
        virtual void put_to_best() {}; /*!< This function main purpose is to prepare the population for optimization
                                        *after the population is initialized and the first calculation of the fitness function.*/
        virtual void combination() {}; /*!< This function generates new candidates that competes with the existing population.
                                        * This is one potential bottleneck in the communication between processors.*/
        virtual void selection() {}; /*!< This function performs selects the candidates for the next round of iteration.*/
        virtual void write_param(double *param_array) {}; /*!< This function let the user specify new values to the search parameters.*/
        virtual void read_param(double *param_array) {}; /*!< This function reads the search parameters that are currently in use.*/
        virtual void fit_to_global() {}; /*!< This function copies memories on best array to global array.*/
        virtual void find_global() {}; /*!< This function finds global best solution within a subset of the population and copy it to the global array.*/


        Problem* prob;

        //functions and variables that are not algorithm specific
        void Init_population(int psize); /*!< A function to initialize the population using a uniform probability distribution.*/
        void Init_previous(double prev_dev, double new_dev, int psize, double *prev_soln); /*!< A function to initialize the population using the solution from N-1 variables.*/

        void Cont_fitness(int p); /*!< A function to compute the fitness value of the solution in contender array.*/
        void Best_fitness(int p); /*!< A function to compute the fitness value of the solution in best array.*/
        void update_popfit(); /*!< A function to compute the fitness value of the solution in best array and update the mean fitness value.*/

        void set_success(int iter, bool goal); /*!<A function to set the condition for which the optimization algorithm would accept the solution and stop.*/
        bool check_success(int t, double *current_fitarray, double *memory_fitarray, int data_size, double t_goal, bool *mem_ptype, int *numvar, int N_cut, double *memory_forT); /*!< A function for checking whether the current solution/iteration satisfied the accept condition.*/

        //Selecting solution
        double Final_select(double *fit, double *solution, double *fitarray); /*!< A function to find the candidate with the highest fitness value in the population.*/
        double avg_Final_select(double* solution, int repeat, double *soln_fit, double *fitarray); /*!< A function to compute the fitness function for 10 times and find the candidate with the highest fitness value in the population*/

        void dev_gen(double *dev_array, double prev_dev, double new_dev, int cut_off); /*!< A function for generating the array that indicates the width of the Gaussian distribution for each of the variable.
                                                                                        * This array is used by Init_previous(). */
        int find_max(double *fit); /*!< A function to find the candidate that has the highet fitness value
                                    * (for only one of the many fitness values that a candidate posseses).*/

        bool success; /*!< Variable for indicating whether a solution is accepted.*/
        bool policy_type; /*!< Variable for indicating whether the policy can be used. */

    protected:
        int num; /*!< Number of variables*/
        int num_fit; /*!< Number of fitness values*/
        Rng *gaussian_rng; /*!< RNG object pointer*/
        int pop_size; /*!< The number of candidates on a processor.*/
        int T; /*!< The goal in number of iteration.*/
        int t; /*!<Time variables*/
        Candidate *pop; /*!< Pointer to population.*/
        bool goal; /*!< The variables indicating whether an accept-reject criteria is used.*/

        int total_pop; /*!< The total number of candidates in the population.*/
        int my_rank; /*!< Variable for MPI*/
        int nb_proc;/*!< Variables for MPI*/

    };

/*! \brief DE class contains the functions that are specific to this particular optimization algorithm.*/

class DE : public OptAlg {
    public:
        DE() {};
        DE(Problem *problem_ptr, Rng *gaussian_rng, int pop_size): OptAlg(problem_ptr, gaussian_rng, pop_size), F(0.1), Cr(0.6) {};
        ~DE() {};

        void put_to_best();
        void combination();
        void selection();
        void fit_to_global();
        void find_global() {};
        void write_param(double *param_array);
        void read_param(double *param_array);

    private:
        double F; /*!< The search variable used in combining the base with the differential part.*/
        double Cr; /*!< The crossover rate.*/

        void family_gen(int* fam, int p, int fam_size); /*!< This function generates three random numbers indicating the candidates
                                                        * that are used by combination() to generate the next generation of candidate.*/

    };

/*! \brief PSO class contains the functions that are specific to this particular optimization algorithm.*/

class PSO : public OptAlg {
    public:
        PSO() {};
        PSO(Problem *problem_ptr, Rng *gaussian_rng, int pop_size): OptAlg(problem_ptr, gaussian_rng, pop_size), w(0.8), phi1(0.6), phi2(1.0), v_max(0.2) {};
        ~PSO() {};

        void put_to_best();
        void combination();
        void selection();
        void write_param(double *param_array);
        void read_param(double *param_array);
        void fit_to_global() {};
        void find_global();

    private:
        inline void find_index(int *prev, int *forw, int p); /*!< A function for finding the index of neighbors of a candidate.*/
        inline int find_fitness(int prev, double prev_fit, int forw, double forw_fit, int p, double fit); /*!< A functoin for getting the index of the candidate in the neighborhood that has the highest fitness value.*/
        double w; /*!< The inertia variable to be multipled to the velocity.*/
        double phi1; /*!< The local search variable.*/
        double phi2; /*!< The global search variable.*/
        double v_max; /*!< The maximum speed for every variable.*/
    };

#endif // OPTALG_H
