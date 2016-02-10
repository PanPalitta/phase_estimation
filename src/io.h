#ifndef IO_H
#define IO_H

void output_header(char const *output_filename, char const *time_filename);

template<typename typeT>
void output_result(int num, double final_fit, typeT *solution,
                   time_t start_time, char const *output_filename,
                   char const *time_filename);

void read_config_file(char const *filename, int *pop_size, int *N_begin,
                      int *N_cut, int *N_end, int *iter, int *iter_begin,
                      int *repeat, int *seed, string *output_filename,
                      string *time_filename, string *optimization,
                      double *prev_dev, double *new_dev, double *t_goal, int *data_end);

#endif // IO_H
