#ifndef PHASE_LOSS_H
#define PHASE_LOSS_H

#include <complex>
#if HAVE_CONFIG_H
#include <config.h>
#endif

#define DEV_N 0.0 //level of noise in operator
//Normal and Hat noise parameters
#define THETA_DEV 1.41421356237//M_PI;//phase noise level
//RTN parameter
#define Ps 0.5 //probability of telegraph noise
//SND parameter
#define RATIO 1.176959769 //ratio between alpha and sigma
//Lognormal parameters
#define MU 2.2214//variance
#define THETA 0.2715 //skewness, variance


#include "problem.h"
#include "rng.h"
#include "aux_functions.h"

typedef complex<double> dcmplx;

/*! \brief Phase class for the problem of adaptive interferometric phase estimation including noise and loss.
*
* This class can be replaced by any optimization problem of interest written by the user.
* The problem is initialized in the main function through the Problem class pointer.
*/

class Phase: public Problem {
    public:
        Phase(const int numvar, Rng *gaussian_rng, Rng *uniform_rng);
        ~Phase();

        void fitness(double *soln, double *fitarray);
        void avg_fitness(double *soln, const int K, double *fitarray);
        bool T_condition(double *fitarray, int *numvar, int N_cut, bool *mem_ptype, double *memory_forT);
        bool error_condition(double *current_fitarray, double *memory_fitarray, int data_size, double goal);
        void boundary(double *can1);

    private:
        double lower; //The lower bound of the variables
        double upper; //The upper bound of the variables
        double loss; //The photon loss rate

        Rng *gaussian_rng, *uniform_rng;

        //variables for WK state generation
        dcmplx *input_state; //The array containing the WK state
        double *sqrtfac_mat; //matrix to keep values of square roots of factorials
        double *overfac_mat; //matrix to keep values of one over factorials
        double *sqrt_cache; //array to avoid calculation of expensive sqrt calls for integers
        double tan_beta;

        dcmplx *state; //state for use with measurement
        dcmplx *update0; //variables for noise_output function
        dcmplx *update1;

        //functions to generate WK_state
        inline void sqrtfac(double *fac_mat); //A function for calculating square root of factorials for state generation.
        inline void one_over_fac(double *over_mat); //A function for calculating one over square root of factorials for state generation.
        inline double cal_spart(const int n, const int k, const int N);//N is the same as total number of photon 'num'.
        void WK_state(); //A funtion generating the WK state
        //Measurement functions
        inline bool noise_outcome(const double phi, const double PHI, const int N); //A function for simulating a photon going through a noisy Mach-Zehnder interferometer.
        inline void state_loss(const int N); //A function for simulating state change under loss of a photon.
        inline double mod_2PI(double PHI); //A function to perform modulo 2PI on phase.

        bool check_policy(double error, double sharp); //A function for checking whether the policy is resilient to loss. This is called by the T_condition().

	double rand_Gaussian(double mean, double dev);
	double rand_Hat(double PHI, double dev);
	//lognormal-distribution noise
	inline double inv_erf(double x); /*!< Function for calculating the inverse of an error function.*/
	inline int sgn(double x); /*!< Sign function.*/
	double Lognormal(double mu, double sigma, double peak);
	//Random telegraph noise
	double rand_RTN(double PHI,double ps,double dev);
	//Skewed normal noise
	double rand_skewed(double mean, double dev, double ratio);

    };
#endif // PHASE_H
