#include <iostream>
#include <cuda.h>
#include <curand.h>

#include "rng_gpu.h"

using namespace std;

double *dev_urandom_numbers;
double *dev_grandom_numbers;
curandGenerator_t gen;

#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
     cout << "CUDA call error at" << __FILE__<< ":" << __LINE__ << endl;\
     }} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
     cout << "CURAND call error at" << __FILE__<< ":" << __LINE__ << endl;\
     }} while(0)

void gpu_cache_alloc(int n_urandom_numbers, int n_grandom_numbers) {
    CUDA_CALL(cudaMalloc((void **)&dev_urandom_numbers,
              n_urandom_numbers*sizeof(double)));
    CUDA_CALL(cudaMalloc((void **)&dev_grandom_numbers,
              n_grandom_numbers*sizeof(double)));
    /* Create pseudo-random number generator */
    CURAND_CALL(curandCreateGenerator(&gen,
                CURAND_RNG_PSEUDO_DEFAULT));
    /* Set seed */
    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen,
                1234ULL));
}

void gpu_cache_init(double *urandom_numbers, int n_urandom_numbers, 
                    double *grandom_numbers, int n_grandom_numbers) {
  /* Generate n floats on device */
  CURAND_CALL(curandGenerateUniformDouble(gen, dev_urandom_numbers, 
                                          n_urandom_numbers));
  CUDA_CALL(cudaMemcpy(urandom_numbers, dev_urandom_numbers, 
                       n_urandom_numbers * sizeof(double),
                       cudaMemcpyDeviceToHost));
  CURAND_CALL(curandGenerateNormalDouble(gen, dev_grandom_numbers, 
                                         n_grandom_numbers, 0.0, 1.0));
  CUDA_CALL(cudaMemcpy(grandom_numbers, dev_grandom_numbers, 
                       n_grandom_numbers * sizeof(double),
                       cudaMemcpyDeviceToHost));
}

void gpu_cache_free() {
    CURAND_CALL(curandDestroyGenerator(gen));
    CUDA_CALL(cudaFree(dev_urandom_numbers));
    CUDA_CALL(cudaFree(dev_grandom_numbers));
}
