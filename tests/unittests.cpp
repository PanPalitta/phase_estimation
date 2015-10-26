#include "unittest++/UnitTest++.h"
#include "rng.h"
#include "problem.h"
#include "phase_loss_opt.h"

#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <complex>

using namespace std;
typedef complex<double> dcmplx;

// run all tests
int main(int argc, char **argv) {
    return UnitTest::RunAllTests();
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
