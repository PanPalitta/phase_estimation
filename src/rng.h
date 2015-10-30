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
        RngBase(bool _gaussian): gaussian(_gaussian) {};
        ~RngBase() {};
        double next_rand(const double mean, const double dev) {
            return 0.0;
            };
    
    protected:
        bool gaussian;
    };

class RngSimple: public RngBase {
    public:
        RngSimple(bool _gaussian, int n_random_numbers, int seed, int rank);
        ~RngSimple();
        double next_rand(const double mean, const double dev);
    };

class RngVectorized: public RngBase {
    public:
        RngVectorized(bool _gaussian): RngBase(_gaussian) {};
        ~RngVectorized() {};
        double next_rand(const double mean, const double dev);
    };

#ifdef CUDA
class RngGpu: public RngVectorized {
    public:
        RngGpu(bool _gaussian, int n_random_numbers, int seed, int rank);
        ~RngGpu();
        double next_rand(const double mean, const double dev);

    private:
        double *random_numbers;
        int n_random_numbers;
        int index_random_numbers;

        double *dev_random_numbers;
        curandGenerator_t gen;
    };

#define Rng RngGpu

#elif defined(VSL)

class RngVsl: public RngVectorized {
    public:
        RngVsl(bool _gaussian, int n_random_numbers, int seed, int rank);
        ~RngVsl();
        double next_rand(const double mean, const double dev);

    private:
        double *random_numbers;
        int n_random_numbers;
        int index_random_numbers;

    }
#define Rng RngGpu

#else

#define Rng RngSimple

#endif

#endif // RNG_H
