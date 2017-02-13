#include "aux_functions.h"
#include <iostream>

using namespace std;

/*############################Linear Regression############################*/
void linear_fit(int data_size, double *x, double *y, double *slope, double *intercept, double *mean_x) {
    /** This function calculates the coefficients of a linear equation.
    * The function takes two arrays of x and y of size data_size, and output the slope, intercept, and the mean_x.
    */
    if(data_size <= 0) {
        throw invalid_argument("data_size must be positive.");
        }
    double sum_x = 0;
    double sum_xx = 0;
    double sum_y = 0;
    double sum_xy = 0;

    //Preparing the data for calculating mean_x, slope, and intercept.
    for(int v = 0; v < data_size; ++v) {
        sum_x = sum_x + x[v];
        sum_xx = sum_xx + x[v] * x[v];
        sum_y = sum_y + y[v];
        sum_xy = sum_xy + x[v] * y[v];
        }

    *mean_x = sum_x / double(data_size);
    *slope = (sum_xy - sum_y * sum_x / double(data_size)) / (sum_xx - sum_x * sum_x / double(data_size));
    *intercept = sum_y / double(data_size) - *slope * (*mean_x);
    }

double error_interval(double *x, double *y, double mean_x, int data_size, double slope, double intercept) {
    /** This function calculates the error of the latest data y in the array.
    * This requires all data x and y of size data_size, slope, mean_x, and intercept from linear_fit function.
    * SSres is a variable used to store the value of (x-mean_x)^2, which can be used to reduce the time for calculating error_update.
    */
    if(data_size <= 0) {
        throw invalid_argument("data_size must be positive.");
        }
    double SSx = 0;
    double SSres = 0;   
 
    for(int i = 0; i < data_size; ++i) {
        SSres = SSres + pow(y[i] - slope * x[i] - intercept, 2);
        SSx = SSx + (x[i] - mean_x) * (x[i] - mean_x);
        }
    return sqrt(SSres / double(data_size - 1) * (1 / data_size + (pow(x[data_size] - mean_x, 2) / SSx)));
    }

double error_update(int old_size, double *SSres, double *mean_x, double slope, double intercept, double *y, double *x) {
    /** This function update the error_interval when a new value of y.
    * The function requires the coefficients of the linear equations and the data array x and y.
    */
    if(old_size <= 0) {
        throw invalid_argument("data_size must be positive.");
        }
    double SSx = 0;

    *mean_x = (*mean_x * old_size + x[old_size]) / double(old_size + 1); //updating mean_x
    *SSres = *SSres + pow(y[old_size] - slope * x[old_size] - intercept, 2);

    for(int i = 0; i < old_size + 1; ++i) {
        SSx = SSx + (x[i] - *mean_x) * (x[i] - *mean_x);
        }
    return sqrt(*SSres / double(old_size - 1) * (1 / old_size + (pow(x[old_size] - *mean_x, 2) / SSx)));
    }

/*###########################Calculating Quantile############################*/
double quantile(double p) { //p is percentile
    double result;
    try {
        result = inv_erf(p);
        }
    catch(invalid_argument) {
        p = 0.999999;
        }
    return sqrt(2) * inv_erf(2 * p - 1);
    }

inline double inv_erf(double x) {
    if(x == 1) {
        throw invalid_argument("Input leads to error.");
        }
    double a = 0.140012;
    double lnx = log(1 - x * x);
    double temp = sqrt(pow(2.0 / (M_PI * a) + lnx / 2.0, 2) - lnx / a);

    return sgn(x) * sqrt(temp - 2.0 / (M_PI * a) - lnx / 2.0);
    }

inline int sgn(double x) {
    /** Sign function of x: https://en.wikipedia.org/wiki/Sign_function
    */
    if(x < 0) {
        return -1;
        }
    else if(x == 0) {
        return 0;
        }
    else {
        return 1;
        }
    }
