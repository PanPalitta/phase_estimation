#include <C:\UnitTest++-1.3\src\UnitTest++.h> // change to directory for UnitTest++
#include <stdexcept>
#include "rng.h"
#include "problem.h"
#include "phase_loss_opt.h"
#include "mpi_optalg.h"

//When the only thing that we can test are things in public,
//this set of tests is no longer providing any useful information


// run all tests
int main(int argc, char **argv) {
    return UnitTest::RunAllTests();
    }

SUITE(pso_construct) {
    TEST(check_num) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new PSO(problem);
        CHECK(alg->num == numvar);
        }
    TEST(check_numfit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new PSO(problem);
        CHECK(alg->num_fit == problem->num_fit);
        }
    }
SUITE(pso_param) {
    TEST(read_param) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new PSO(problem);

        double total = 0;
        double input[4] = {0.8, 0.6, 1.0, 0.2};
        double output[4];
        alg->read_param(output);
        for(int i = 0; i < 4; ++i) {
            total += input[i] - output[i];
            }

        CHECK(total == 0);
        }

    TEST(write_param) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new PSO(problem);

        double total = 0;
        double input[4] = {1.7, 2.6, 5.6, 0.9};
        double output[4];
        alg->write_param(input);
        alg->read_param(output);
        for(int i = 0; i < 4; ++i) {
            total += input[i] - output[i];
            }
        CHECK(total == 0);
        }
    }
