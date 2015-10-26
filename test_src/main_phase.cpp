#include <C:\UnitTest++-1.3\src\UnitTest++.h>
#include "rng.h"
#include "problem.h"
#include "phase_loss_opt.h"

#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <complex>

using namespace std;
typedef complex<double> dcmplx;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// To add a test, simply put the following code in the a .cpp file of your choice:
//
// =================================
// Simple Test
// =================================
//
//  TEST(YourTestName)
//  {
//  }
//
// The TEST macro contains enough machinery to turn this slightly odd-looking syntax into legal C++, and automatically register the test in a global list.
// This test list forms the basis of what is executed by RunAllTests().
//
// If you want to re-use a set of test data for more than one test, or provide setup/teardown for tests,
// you can use the TEST_FIXTURE macro instead. The macro requires that you pass it a class name that it will instantiate, so any setup and teardown code should be in its constructor and destructor.
//
//  struct SomeFixture
//  {
//    SomeFixture() { /* some setup */ }
//    ~SomeFixture() { /* some teardown */ }
//
//    int testData;
//  };
//
//  TEST_FIXTURE(SomeFixture, YourTestName)
//  {
//    int temp = testData;
//  }
//
// =================================
// Test Suites
// =================================
//
// Tests can be grouped into suites, using the SUITE macro. A suite serves as a namespace for test names, so that the same test name can be used in two difference contexts.
//
//  SUITE(YourSuiteName)
//  {
//    TEST(YourTestName)
//    {
//    }
//
//    TEST(YourOtherTestName)
//    {
//    }
//  }
//
// This will place the tests into a C++ namespace called YourSuiteName, and make the suite name available to UnitTest++.
// RunAllTests() can be called for a specific suite name, so you can use this to build named groups of tests to be run together.
// Note how members of the fixture are used as if they are a part of the test, since the macro-generated test class derives from the provided fixture class.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// run all tests
int main(int argc, char **argv) {
    return UnitTest::RunAllTests();
    }

SUITE(RNG) {
    TEST(rngen_uniform) {
        int n = 10;
        int u = 10;
        int seed = time(NULL);
        int rank = 0;
        double random1, random2;
        srand(seed + rank);
        random1 = double(rand()) / RAND_MAX;
        Rng *rng = new Rng(n, u, seed, rank);
        random2 = rng->next_urand();
        CHECK(random1 == random2);
        }

    TEST(rngen) {
        int n = 10;
        int u = 10;
        int seed = time(NULL);
        int rank = 0;
        double dev = -0.1;
        double mean = 1;
        double random1, random2;
        srand(seed + rank);
        random1 = (double(rand()) / RAND_MAX - 0.5) * 2 * dev + mean;
        Rng *rng = new Rng(n, u, seed, rank);
        random2 = rng->next_grand(mean, dev);
        CHECK(random1 == random2);
        }
    }

SUITE(phase_mainfunc) {
    TEST(constructor) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = -5;
        bool err = 0;
        Problem* problem;
        CHECK_THROW(problem = new Phase(numvar, rng), invalid_argument);
        }
    TEST(test_catch) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = -5;
        Problem* problem;
        bool err = false;
        try {
            problem = new Phase(numvar, rng);
            //This test has a problem b/c destructor calls a delete[] when the array hasn't yet been inititalized.
            }
        catch(invalid_argument) {
            err = true;
            }
        CHECK(err == true);
        }
    TEST(boundary) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 10;
        Phase* problem = new Phase(numvar, rng);
        double Phi = 4 * M_PI / (numvar + 1);
        double array[numvar];
        bool err = 0;

        for(int i = 0; i < numvar; ++i) {
            array[i] = (i - numvar / 2) * Phi;
            }
        problem->boundary(array);
        for(int i = 0; i < numvar; ++i) {
            if(array[i] > problem->upper_bound[i]) {
                err = 1;
                }
            else if(array[i] < problem->lower_bound[i]) {
                err = 1;
                }
            else {}
            }

        CHECK(err == 0);
        }
    TEST(fitness) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 10;
        Problem* problem = new Phase(numvar, rng);
        double can[numvar + 1];
        int K = numvar * numvar * 10;
        double fit[problem->num_fit];

        for(int i = 0; i < numvar; ++i) {
            can[i] = 1;
            }
        problem->avg_fitness(can, K, fit);
        CHECK(fit[0] < 1.0);
        }

    }

SUITE(phase_auxfunc) {
    TEST(mod_2PI) {
        int n = 10;
        int u = 10;
        int seed = time(NULL);
        int rank = 0;
        Rng *rng = new Rng(n, u, seed, rank);
        int numvar = 4;
        Phase* problem = new Phase(numvar, rng);
        double result;
        double phi = -3 * M_PI;
        double modifier = 0;
        result = problem->mod_2PI(phi + modifier);
        CHECK(fabs(result - (phi + 4 * M_PI)) < 0.000000000000001);
        }
    TEST(sqrtfac) {
        int n = 10;
        int u = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(n, u, seed, rank);
        int numvar = 5;
        Phase* problem = new Phase(numvar, rng);
        double array[numvar];
        double final = 1.0;
        problem->sqrtfac(array);

        for(int i = 1; i <= numvar; ++i) {
            final = final * sqrt(i);
            }
        CHECK(final == array[numvar]);
        }
    TEST(one_over_fac) {
        int n = 10;
        int u = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(n, u, seed, rank);
        int numvar = 5;
        Phase* problem = new Phase(numvar, rng);
        double array[numvar + 1];
        double final = 1.0;

        problem->one_over_fac(array);

        for(int i = 1; i <= numvar - 1; ++i) {
            final = final / i;
            }
        CHECK(final == array[numvar - 1]);
        }
    TEST(cal_spart) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 10;

        int n = 0;
        int k = 0;
        double result;
        Phase* problem = new Phase(numvar, rng);
        problem->one_over_fac(problem->overfac_mat);
        problem->tan_beta = tan(M_PI / 4);
        result = problem->cal_spart(n, k, numvar);

        for(int i = 1; i <= numvar; ++i) {
            result = result * i;
            }
        CHECK(result - 1 <= 0.000001);
        }
    }

SUITE(phase_operations) {
    TEST(state_loss) {
        int n = 10;
        int u = 10;
        int seed = time(NULL);
        int rank = 5;
        Rng *rng = new Rng(n, u, seed, rank);
        int numvar = 4;
        Phase* problem = new Phase(numvar, rng);
        double total;
        //create state
        for(int i = 0; i <= numvar; i++) {
            problem->state[i].real() = 1 / sqrt(numvar + 1);
            problem->state[i].imag() = 0;
            }
        problem->state_loss(numvar);
        total = 0.0;
        for(int i = 0; i <= numvar; i++) {
            total += problem->state[i].real() * problem->state[i].real();
            }
        CHECK(fabs(total - 1) <= 0.000000000000001);
        }
    TEST(noisy_outcome) {
        int n = 10;
        int u = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(n, u, seed, rank);
        int numvar = 4;
        Phase* problem = new Phase(numvar, rng);
        int result;
        double total = 0;

        double phi = -M_PI / 2;
        double PHI = 0;
        //state generation
        for (int i = 0; i <= numvar; ++i) {
            problem->state[i].real() = 1 / sqrt(numvar + 1);
            problem->state[i].imag() = 0;
            }
        result = problem->noise_outcome(phi, PHI, numvar);
        for (int i = 0; i <= numvar; ++i) {
            total += problem->state[i].real() * problem->state[i].real() + problem->state[i].imag() * problem->state[i].imag();
            }
        CHECK(fabs(total - 1) <= 0.000000000000001);
        }
    TEST(WK_state) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 11;
        Phase* problem = new Phase(numvar, rng);
        double total = 0;
        for(int i = 0; i <= numvar; ++i) {
            total += abs(problem->input_state[i] * conj(problem->input_state[i]));
            }
        CHECK(total - 1 <= 0.0000000000000001);
        }
    }
