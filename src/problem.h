#ifndef PROBLEM_H
#define PROBLEM_H

using namespace std;

/*! \brief Problem class contains the prototype of the functions in the optimization problem that OptAlg class needs.
*/

class Problem {
    public:
        Problem() {};
        virtual ~Problem() {
            delete[] upper_bound;
            delete[] lower_bound;
            }

        virtual void fitness(double *soln, double *fitarray) {
            /*! A function intend to be a wrapper for changing conditions in which the fitness function is evaluated.*/
            }
        virtual void avg_fitness(double *soln, int K, double *fitarray) {
            /*! A function for calculating the fitness value.
            It allows a number of sample K to be passed into the function in case the fitness function is a statistical 'quantity'(?)*/
            }
        virtual bool T_condition(double *fitarray, int *numvar, int N_cut, bool *mem_ptype, double *memory_forT) {
            /*! A function for calculating additional conditions for when the optimization algorithm is set to accept solution after time T.*/
            return 0;
            }
        virtual bool error_condition(double *current_fitarray, double *memory_fitarray, int data_size, double goal) {
            /*! A function for calculating additional conditions for when optimization algorithm is set to accept solution from error bound.*/
            return 0;
            }
        virtual void boundary(double *can1) {
            /*! This function is used to keep the solution candidate within the boundary of the search space.*/
            }

        double *lower_bound;/*!< Pointer to array storing the lower bound of the variables*/
        double *upper_bound;/*!< Pointer to array storing the upper bound of the variables*/

        int num; /*!<number of variables in the problem*/
        int num_repeat; /*!<number of repeats*/
        int num_fit; /*<number of fitnesses*/
    };

#endif // PROBLEM_H
