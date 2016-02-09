#include <iostream>
#include <map>
#include <vector>
#include <mpi.h>
#include <cuda.h>
#include <curand.h>

#include "rng.h"

using namespace std;

#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
            cout << "CUDA call error at" << __FILE__<< ":" << __LINE__ << endl;\
            }} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
            cout << "CURAND call error at" << __FILE__<< ":" << __LINE__ << endl;\
            exit(-1);\
            }} while(0)

/// Note that this function was lifted from http://code.google.com/p/gpmr/
void setDevice(int commRank, int commSize) {
    int devCount;
    int deviceNum = 0;
    CUDA_CALL(cudaGetDeviceCount(&devCount));
    FILE * fp = popen("/bin/hostname", "r");
    char buf[1024];
    if (fgets(buf, 1023, fp) == NULL) strcpy(buf, "localhost");
    pclose(fp);
    string host = buf;
    host = host.substr(0, host.size() - 1);
    strcpy(buf, host.c_str());
    if (commRank == 0) {
        map<string, vector<int> > hosts;
        map<string, int> devCounts;
        hosts[buf].push_back(0);
        devCounts[buf] = devCount;

        MPI_Status stat;
        MPI_Request req;
        for (int i = 1; i < commSize; ++i) {
            MPI_Recv(buf, 1024, MPI_CHAR, i, 0, MPI_COMM_WORLD, &stat);
            MPI_Recv(&devCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &stat);

            // check to make sure each process on each node reports the same number of devices.
            hosts[buf].push_back(i);
            if (devCounts.find(buf) != devCounts.end()) {
                if (devCounts[buf] != devCount) {
                    printf("Error, device count mismatch %d != %d on %s\n", devCounts[buf], devCount, buf);
                    fflush(stdout);
                    }
                }
            else devCounts[buf] = devCount;
            }
        // check to make sure that we don't have more jobs on a node than we have GPUs.
        for (map<string, vector<int> >::iterator it = hosts.begin(); it != hosts.end(); ++it) {
            if (it->second.size() > static_cast<unsigned int>(devCounts[it->first])) {
                printf("Error, more jobs running on '%s' than devices - %d jobs > %d devices.\n",
                       it->first.c_str(), static_cast<int>(it->second.size()), devCounts[it->first]);
                fflush(stdout);
                exit(1);
                }
            }

        // send out the device number for each process to use.
        MPI_Irecv(&deviceNum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req);
        for (map<string, vector<int> >::iterator it = hosts.begin(); it != hosts.end(); ++it) {
            for (unsigned int i = 0; i < it->second.size(); ++i) {
                int devID = i;
                MPI_Send(&devID, 1, MPI_INT, it->second[i], 0, MPI_COMM_WORLD);
                }
            }
        MPI_Wait(&req, &stat);
        }
    else {
        // send out the hostname and device count for your local node, then get back the device number you should use.
        MPI_Status stat;
        MPI_Send(buf, strlen(buf) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&devCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&deviceNum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        }
    MPI_Barrier(MPI_COMM_WORLD);
    CUDA_CALL(cudaSetDevice(deviceNum));
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp devProp;
    cudaGetDeviceProperties(&devProp, device);
    cout <<  device << " " << devProp.name << " Compute Capability: " << devProp.major << "." << devProp.minor << "\n";
}

RngGpu::RngGpu(bool _gaussian, int _n_random_numbers, int seed, int rank):
    n_random_numbers(_n_random_numbers), RngVectorized(_gaussian) {
    int nb_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    setDevice(rank, nb_proc);
    CUDA_CALL(cudaMallocHost((void **)&random_numbers,
                             n_random_numbers * sizeof(double)));
    index_random_numbers = n_random_numbers;
    CUDA_CALL(cudaMalloc((void **)&dev_random_numbers,
                         n_random_numbers * sizeof(double)));
    /* Create pseudo-random number generator */
    CURAND_CALL(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
    /* Set seed */
    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, seed));
    if (gaussian) {
        CURAND_CALL(curandGenerateNormalDouble(gen, dev_random_numbers,
                                               n_random_numbers, 0.0, 1.0));
    } else {
        CURAND_CALL(curandGenerateUniformDouble(gen, dev_random_numbers,
                                                n_random_numbers));
    }
    CUDA_CALL(cudaMemcpyAsync(random_numbers, dev_random_numbers,
                              n_random_numbers * sizeof(double),
                              cudaMemcpyDeviceToHost));
    }


double RngGpu::next_rand(const double mean, const double dev) {
    if (index_random_numbers >= n_random_numbers) {
        index_random_numbers = 0;
        cudaDeviceSynchronize();
        if (gaussian) {
            CURAND_CALL(curandGenerateNormalDouble(gen, dev_random_numbers,
                                                   n_random_numbers, 0.0, 1.0));
        } else {
            CURAND_CALL(curandGenerateUniformDouble(gen, dev_random_numbers,
                                                    n_random_numbers));
        }
        CUDA_CALL(cudaMemcpyAsync(random_numbers, dev_random_numbers,
                                  n_random_numbers * sizeof(double),
                                  cudaMemcpyDeviceToHost));
        }
    return random_numbers[index_random_numbers++] * dev + mean;
    }

RngGpu::~RngGpu() {
    CURAND_CALL(curandDestroyGenerator(gen));
    CUDA_CALL(cudaFree(dev_random_numbers));
    CUDA_CALL(cudaFreeHost(random_numbers));
    }
