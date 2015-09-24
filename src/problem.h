#ifndef PROBLEM_H
#define PROBLEM_H

using namespace std;

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
            //return 0;   //K is the number of samples to calculate the average
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
