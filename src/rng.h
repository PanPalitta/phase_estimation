#ifndef RNG_H
#define RNG_H
#if HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef CUDA
#include <curand.h>
#endif

using namespace std;

/*!  \brief Rng class store methods for generating random numbers using RngBase as its base class.
*    The class can generate uniformaly random number for approximating normally distributed numbers when 'gaussian' is specified.
*/

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

/*!   \brief RngSimple generates random numbers from C++ build-in pseudo-random number generator.
*/

class RngSimple: public RngBase {
    public:
        RngSimple(bool _gaussian, int n_random_numbers, int seed, int rank);
        ~RngSimple();
        double next_rand(const double mean, const double dev);
    };

/*! \brief RngVectorized generates random numbers into vectors in order to reduce the computational overhead.
*   There are two methods to this class: VSL and GPU.
*/

class RngVectorized: public RngBase {
    public:
        RngVectorized(bool _gaussian): RngBase(_gaussian) {};
        ~RngVectorized() {};
        double next_rand(const double mean, const double dev);
    };

#ifdef CUDA
/*! \brief RngGpu generates vectors of random number asynchronously on GPUs. One GPU is assigned per CPU.
*/

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

/*!    \brief RngVsl uses Intel's VSL library to generate a vector of random numbers. This approach reduce computational overhead.
*/

class RngVsl: public RngVectorized {
    public:
        RngVsl(bool _gaussian, int n_random_numbers, int seed, int rank);
        ~RngVsl();
        double next_rand(const double mean, const double dev);

    private:
        double *random_numbers;
        int n_random_numbers;
        int index_random_numbers;
    };
#define Rng RngVsl

#else

#define Rng RngSimple

#endif

#endif // RNG_H
