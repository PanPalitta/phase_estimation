#ifndef RNG_VSL_H
#define RNG_VSL_H

void vsl_init();
void vsl_cache_init(double *urandom_numbers, int n_urandom_numbers,
                    double *grandom_numbers, int n_grandom_numbers);
void vsl_close();
#endif
