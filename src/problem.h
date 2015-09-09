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
            //return 0;
            }
        virtual void avg_fitness(double *soln, int K, double *fitarray) {
            //return 0;   //K is the number of samples to calculate the average
            }

        double *lower_bound;
        double *upper_bound;

        int num; //number of variables in the problem
        int num_repeat; //number of repeats (calculated from num of var)
        int num_fit; //number of fitness

        /*functions that don't need to be linked to specific problems*/
        void modulo(double *can1);
        void normalize(double *can1);

    };

#endif // PROBLEM_H
