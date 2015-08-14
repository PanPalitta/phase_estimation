#ifndef PROBLEM_H
#define PROBLEM_H

#include <complex>

using namespace std;
typedef complex<double> dcmplx;

template<class typeT>
class Problem
{
public:
    Problem() {};
    virtual ~Problem() {
        delete[] upper_bound;
        delete[] lower_bound;
    }

    virtual void fitness(typeT *soln, double *fitarray) {
        //return 0;
    }
    virtual void avg_fitness(typeT *soln, int K, double *fitarray) {
        //return 0;   //K is the number of samples to calculate the average
    }

    typeT *lower_bound;
    typeT *upper_bound;

    int num; //number of variables in the problem
    int num_repeat; //number of repeats (calculated from num of var)
    int num_fit; //number of fitness

    /*functions that don't need to be linked to specific problems*/
    void modulo(typeT *can1);
    void normalize(typeT *can1);

};

template<typename typeT>
void Problem<typeT>::modulo(typeT *can1) {
    int i;

    for(i=0; i<num; ++i) {
        if(dcmplx(can1[i]).real()<dcmplx(lower_bound[i]).real()) {
            dcmplx(can1[i]).real()=dcmplx(lower_bound[i]).real();
        }
        else if (dcmplx(can1[i]).real()>dcmplx(upper_bound[i]).real()) {
            dcmplx(can1[i]).real()=dcmplx(upper_bound[i]).real();
        }
        else {
            /*doesn't change anything*/
        }
        if(dcmplx(can1[i]).imag()!=0) {
            if (dcmplx(can1[i]).imag()<dcmplx(lower_bound[i]).imag()) {
                dcmplx(can1[i]).imag()=dcmplx(lower_bound[i]).imag();
            }
            else if (dcmplx(can1[i]).imag()>dcmplx(upper_bound[i]).imag()) {
                dcmplx(can1[i]).imag()=dcmplx(upper_bound[i]).imag();
            }
            else {
                /*doesn't change anything*/
            }
        }

    }

}

template<typename typeT>
void Problem<typeT>::normalize(typeT *can1) {
    int i;
    double norm=0;
    for(i=0; i<num; ++i) {
        norm+=abs(can1[i]*can1[i]);
    }
    norm=sqrt(norm);
    for(i=0; i<num; ++i) {
        can1[i]=can1[i]/norm;
    }
}

#endif // PROBLEM_H
