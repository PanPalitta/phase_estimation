#ifndef IO_H
#define IO_H

/*! \brief This file contains the functions that are involved in setting the input parameters
 and printing out the results.
*/

void output_header(char const *output_filename, char const *time_filename); /*!<This function print the headers to output files.*/

void output_result(int num, int num_fit, double *final_fit, double *solution,
                   time_t start_time, char const *output_filename,
                   char const *time_filename); /*!<This function print the sharpness, policy, and time for N.*/

void read_config_file(char const *filename, int *pop_size, int *N_begin,
                      int *N_cut, int *N_end, int *iter, int *iter_begin,
                      int *repeat, int *seed, string *output_filename,
                      string *time_filename, string *optimization,
                      double *prev_dev, double *new_dev, double *t_goal, int *data_end);/*!<This function set the parameters for the optimization algorithm.*/

#endif // IO_H
