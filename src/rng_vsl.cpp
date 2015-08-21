#include <iostream>
#include <map>
#include <vector>
#include "mkl_vsl.h"

#include "rng_vsl.h"

#define BRNG VSL_BRNG_MCG31
#define METHOD VSL_RNG_METHOD_GAUSSIAN_ICDF
#define SEED 0

VSLStreamStatePtr stream;

using namespace std;

RngVsl::RngVsl(int _n_urandom_numbers, int _n_grandom_numbers):
	n_urandom_numbers(_n_urandom_numbers),
	n_grandom_numbers(_n_urandom_numbers)
{
	vslNewStream(&stream, BRNG, SEED);
	urandom_numbers = new double[n_urandom_numbers];
	index_urandom_numbers = 0;
	grandom_numbers = new double[n_grandom_numbers];
	index_grandom_numbers = 0;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_urandom_numbers, urandom_numbers, 0.0, 1.0);
	vdRngGaussian(METHOD, stream, n_grandom_numbers, grandom_numbers, 0.0, 1.0);
}

double RngVsl::next_grand(const double mean, const double dev)
{
	if (index_grandom_numbers >= n_grandom_numbers) {
		index_grandom_numbers = 0;
		vdRngGaussian(METHOD, stream, n_grandom_numbers, grandom_numbers, 0.0, 1.0);
	}
	return grandom_numbers[index_grandom_numbers++]*dev+mean;
}

double RngVsl::next_urand()
{
	if (index_urandom_numbers >= n_urandom_numbers) {
		index_urandom_numbers = 0;
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n_urandom_numbers, urandom_numbers, 0.0, 1.0);
	}
	return urandom_numbers[index_urandom_numbers++];
}

RngVsl::~RngVsl()
{
	delete[] grandom_numbers;
	delete[] urandom_numbers;
	vslDeleteStream(&stream);
}
