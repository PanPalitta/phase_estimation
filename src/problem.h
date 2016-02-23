#ifndef PROBLEM_H
#define PROBLEM_H

using namespace std;

/*! \brief Problem class contains the prototype of the functions in the optimization problem that OptAlg class needs to access.
*/

class Problem {
    public:
        Problem() {};
        virtual ~Problem() {
            delete[] upper_bound;
            delete[] lower_bound;
            }

        virtual void fitness(double *soln, double *fitarray) {
            //return 0;//Do we still need this here?
            }
        virtual void avg_fitness(double *soln, int K, double *fitarray) {
            /*! A function for calculating the fitness value.
            It allows a number of sample K to be passed into the function in case the fitness function is a statistical 'quantity'(?)*/
            }
        virtual void T_condition(double *fitarray, int *numvar, int N_cut, bool *mem_ptype) {
            /*! A function for calculating additional conditions for when the optimization algorithm is set to accept solution after time T.*/
            }
        virtual bool error_condition(double *memory_fitarray, int data_size, double t_goal) {
            /*! A function for calculating additional conditions for when optimization algorithm is set to accept solution from error bound.*/
            return 0;
            }
        virtual void boundary(double *can1) {
            }

        double *lower_bound;
        double *upper_bound;

        int num; //number of variables in the problem
        int num_repeat; //number of repeats (calculated from num of var)
        int num_fit; //number of fitness
    };

#endif // PROBLEM_H
