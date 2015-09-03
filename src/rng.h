#ifndef RNG_H
#define RNG_H
#if HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef CUDA
#include <curand.h>
#endif

using namespace std;

class RngBase {
    public:
        ~RngBase() {};
        double next_grand(const double mean, const double dev) {
            return 0.0;
            };
        double next_urand() {
            return 0.0;
            };
    };

class RngSimple: public RngBase {
    public:
        RngSimple(int n_urandom_numbers, int n_grandom_numbers, int seed, 
                  int rank);
        ~RngSimple();
        double next_grand(const double mean, const double dev);
        double next_urand();
    };

class RngVectorized: public RngBase {
    public:
        ~RngVectorized() {};
        double next_grand(const double mean, const double dev);
        double next_urand();

    private:
        double *urandom_numbers;
        int n_urandom_numbers;
        int index_urandom_numbers;
        double *grandom_numbers;
        int n_grandom_numbers;
        int index_grandom_numbers;
    };

#ifdef CUDA
class RngGpu: public RngVectorized {
    public:
        RngGpu(int n_urandom_numbers, int n_grandom_numbers, int seed, int rank);
        ~RngGpu();
        double next_grand(const double mean, const double dev);
        double next_urand();

    private:
        double *urandom_numbers;
        int n_urandom_numbers;
        int index_urandom_numbers;
        double *grandom_numbers;
        int n_grandom_numbers;
        int index_grandom_numbers;

        double *dev_urandom_numbers;
        double *dev_grandom_numbers;
        curandGenerator_t gen;
    };

#define Rng RngGpu

#elif defined(VSL)

class RngVsl: public RngVectorized {
    public:
        RngVsl(int n_urandom_numbers, int n_grandom_numbers, int seed, 
               int rank);
        ~RngVsl();
        double next_grand(const double mean, const double dev);
        double next_urand();

    private:
        double *urandom_numbers;
        int n_urandom_numbers;
        int index_urandom_numbers;
        double *grandom_numbers;
        int n_grandom_numbers;
        int index_grandom_numbers;
    }
#define Rng RngGpu

#else

#define Rng RngSimple

#endif

#endif // RNG_H
