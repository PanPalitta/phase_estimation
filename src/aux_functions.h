#include <stdexcept>
#include <cmath>

using namespace std;

/*! \brief This library contains the functions for linear regression.
*/

//calculating linear regression
void linear_fit(int data_size, double *x, double *y, double *slope, double *intercept, double *mean_x); /*!< Function for calculating the linear regression from a set of x,y array of size data_size.*/
double error_interval(double *x, double *y, double mean_x, int data_size, double slope, double intercept);/*!< Function for calculating the error of the current y data using the previous x and y data and the linear equation calculated from linear_fit.*/
double error_update(int old_size, double *SSres, double *mean_x, double slope, double intercept, double *y, double *x); /*!< Function for updating the error interval for new value of y data. This function has to be used after error_interval.*/
//calculating quantile
double quantile(double p); /*!< Function for calculating the quantile of a normal distribution corresponding to the percentile p.*/
inline double inv_erf(double x); /*!< Function for calculating the inverse of an error function.*/
inline int sgn(double x); /*!< Sign function.*/
