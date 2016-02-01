//#include <cstdlib>
#include <stdexcept>
#include <cmath>

using namespace std;

//function for accept-reject criteria: should be moved to problem class at some point
bool check_type(int t,int T,int *numvar,int N_cut,bool *memptype,double *fitarray);
bool check_policy(double error, double sharp);//specific to phase estimation
//calculating linear regression
void linear_fit(int data_size, double *x, double *y, double *slope, double *intercept, double *mean_x);
double error_interval(double *x, double *y, double mean_x, int data_size, double *SSres, double slope, double intercept);
double error_update(int old_size, double *SSres, double *mean_x, double slope, double intercept, double *y, double *x);
//calculating quantile
double quantile(double p);
inline double inv_erf(double x);
inline int sgn(double x);
