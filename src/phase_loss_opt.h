#ifndef PHASE_LOSS_H
#define PHASE_LOSS_H

#include<complex>
#include<cmath>
#include<iostream>
#include<cstdlib>
#include<ctime>

/*NOTE on using fitness functions for phase estimation problem with loss*/
/* Both avg_fitness() and fitness() contain the same code.
 * The policies are learned without loss (loss in avg_fitness() set to zero).
 * Then the policies are selected based on its mean fitness value for lossy 
 * interferometer (loss in fitness() set to other than zero) which is called 
 * through avg_Final_select() in OptAlg class.
 * */


#include "problem.h"
class Phase: public Problem<double>
{
public:
	Phase(int numvar);
	~Phase();

	double fitness(double *soln);
	double avg_fitness(double *soln, int K);

private:
	double lower;
	double upper;

	//variables for WK state generation
	dcmplx *input_state;
	double *sqrtfac_mat; //matrix to keep values of square roots of factorials
	double *overfac_mat; //matrix to keep values of one over factorials
	double cosN;
	double tan_beta;
	//state for use with measurement
	dcmplx *state;
	dcmplx *update0;
	dcmplx *update1;
	//functions to generate WK_state
	void sqrtfac(double *fac_mat);
	void one_over_fac(double *over_mat);
	double cal_spart(int n, int k, int N);//N is is the same as total number of photon 'num', but I'll leave it like this.
	void WK_state();
	//Measurement function
	bool outcome(double phi, double PHI, int N);//N is the number of photons currently available, not equal to 'num'
	bool noise_outcome(double phi, double PHI, int N);
	void state_loss(int N);
	double rand_Gaussian(double mean,double dev);
	double mod_2PI(double PHI);
};




#endif // PHASE_H
